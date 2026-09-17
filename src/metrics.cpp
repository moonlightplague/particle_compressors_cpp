#include "particle/pipeline.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace particle {
namespace {
Json metric(const Array &original, const Array &restored, double scale) {
  size_t n = original.size();
  double lo = std::numeric_limits<double>::infinity(), hi = -lo, sum = 0,
         squared = 0, maximum = 0;
  for (size_t i = 0; i < n; ++i) {
    double x = double(original.number(i)) / scale,
           y = double(restored.number(i)) / scale, e = y - x;
    lo = std::min(lo, x);
    hi = std::max(hi, x);
    sum += std::abs(e);
    squared += e * e;
    maximum = std::max(maximum, std::abs(e));
  }
  double mse = squared / n, rmse = std::sqrt(mse), range = hi - lo;
  double inf = std::numeric_limits<double>::infinity();
  return {{"count", n},
          {"min", lo},
          {"max", hi},
          {"range", range},
          {"max_absolute_error", maximum},
          {"mean_absolute_error", sum / n},
          {"mse", mse},
          {"mmse", mse},
          {"rmse", rmse},
          {"nrmse", range  ? rmse / range
                    : rmse ? inf
                           : 0},
          {"psnr", mse == 0 ? inf
                   : range == 0
                       ? -inf
                       : 20 * std::log10(range) - 10 * std::log10(mse)}};
}
} // namespace
Json compute_metrics(const Options &, const Json &m) {
  auto started = std::chrono::steady_clock::now();
  size_t count = m.at("count");
  auto source = read_hdf5(m.at("input_h5").get<std::string>(), count);
  auto restored =
      read_hdf5(m["artifacts"]["reconstructed_h5"].get<std::string>());
  if (restored.count() != count)
    throw std::runtime_error("Metrics count mismatch");
  std::vector<size_t> order(count);
  std::iota(order.begin(), order.end(), 0);
  std::string mapping = m.value("ordering", Json::object())
                            .value("reconstructed_rows", Json::object())
                            .value("mapping", "original_row");
  std::string alignment = "original_row";
  auto row_meta = m.value("ordering", Json::object())
                      .value("reconstructed_rows", Json::object());
  std::string artifact =
      row_meta.value("temporary_permutation_artifact", Json(nullptr))
              .is_string()
          ? row_meta["temporary_permutation_artifact"].get<std::string>()
          : "position_order";
  auto pre = m["artifacts"].value("preprocessed", Json::object());
  bool temporary = mapping != "original_row" && pre.contains(artifact) &&
                   fs::is_regular_file(pre[artifact].get<std::string>());
  if (temporary) {
    auto type = type_of(row_meta.value("temporary_permutation_dtype", "int32"));
    Array indices{type, read_bytes(pre[artifact].get<std::string>())};
    if (indices.size() != count)
      throw std::runtime_error("Metric permutation size mismatch");
    std::vector<bool> seen(count, false);
    for (size_t i = 0; i < count; ++i) {
      auto v = indices.number(i);
      if (v < 0 || v >= count || seen[size_t(v)])
        throw std::runtime_error("Invalid metric permutation");
      seen[size_t(v)] = true;
      order[i] = size_t(v);
    }
    alignment = "temporary_" + artifact;
  } else if (mapping == "id_sorted") {
    order = stable_id_order(source.data[0]);
    alignment = "particle_id";
  } else if (mapping != "original_row") {
    auto sorted = stable_id_order(source.data[0]);
    for (size_t i = 1; i < count; ++i)
      if (source.data[0].number(sorted[i - 1]) ==
          source.data[0].number(sorted[i]))
        throw std::runtime_error(
            "Cannot align duplicate particle IDs for metrics");
    for (size_t i = 0; i < count; ++i) {
      auto id = restored.data[0].number(i);
      auto it = std::lower_bound(sorted.begin(), sorted.end(), id,
                                 [&](size_t row, long double v) {
                                   return source.data[0].number(row) < v;
                                 });
      if (it == sorted.end() || source.data[0].number(*it) != id)
        throw std::runtime_error("Metrics ID mismatch");
      order[i] = *it;
    }
    alignment = "particle_id";
  }
  for (auto &a : source.data)
    a = a.gather(order);
  Json report;
  for (auto key : {"compressed_fields", "compressors", "particle_sort",
                   "runtime", "sizes", "timing"})
    report[key] = m.value(key, Json::object());
  report["order_dtype"] = m.value("order_dtype", "int64");
  report["row_comparison"] = {{"reconstructed_order", mapping},
                              {"alignment_source", alignment}};
  double scale = m["position_scale"]["value"];
  for (size_t i = 0; i < 7; ++i) {
    std::string name = fields[i];
    bool pos = i > 0 && i < 4;
    auto metric_values =
        metric(source.data[i], restored.data[i], pos ? scale : 1);
    metric_values["original_dtype"] = source.data[i].type.name;
    metric_values["reconstructed_dtype"] = restored.data[i].type.name;
    if (i == 0)
      metric_values["exact_match"] =
          source.data[i].bytes == restored.data[i].bytes;
    if (pos && !source.data[i].type.floating)
      metric_values["fixed_point_units"] =
          metric(source.data[i], restored.data[i], 1);
    report["fields"][name] = metric_values;
    auto b = m["field_error_bounds"][name];
    double requested = b["abs"], cb = b["compressor_abs"], cast = 0,
           rounding = 0;
    std::string codec = i == 0
                            ? "pcodec"
                            : m["compressors"][pos ? "positions" : "velocities"]
                                  .get<std::string>();
    if (pos) {
      cast = m["preprocess"]["positions"][name]
              ["preprocess_cast_max_abs_in_lcp_units"];
      rounding = source.data[i].type.floating ? 0 : .5 / scale;
    } else if (i >= 4 && codec == "lcp")
      cast = m["preprocess"]["velocities"][name]["preprocess_cast_max_abs"];
    if (codec == "xnyzip") {
      b = m["field_error_bounds"]
           [pos ? "positions_xnyzip" : "velocities_xnyzip"];
      requested = b["abs"];
      cb = b["compressor_abs"];
      cast = b.value("preprocess_l2_max_abs", 0.);
      rounding = 0;
    }
    double effective = b["mode"] == "relative" || codec == "xnyzip"
                           ? requested
                           : cb + cast + rounding;
    if (i == 0) {
      effective = m["error_bounds"]["id_sz3_abs"].get<double>();
      cb = effective;
    }
    double observed = metric_values["max_absolute_error"], range = b["range"];
    Json target = {
        {"mode", b["mode"]},
        {"relative_error_bound", b["relative"]},
        {"range_for_relative", range},
        {"range_units", b["range_units"]},
        {"requested_abs_bound", requested},
        {"compressor_abs_eb", cb},
        {"effective_final_abs_bound", effective},
        {"observed_relative_error",
         range           ? observed / range
         : observed == 0 ? 0
                         : std::numeric_limits<double>::infinity()},
        {"observed_max_absolute_error", observed},
        {"satisfied",
         observed <= effective + 1e-12 + 1e-6 * std::max(1., effective)}};
    if (i) {
      target["preprocess_cast_allowance"] = cast;
      if (pos)
        target["recombine_rounding_allowance"] = rounding;
    }
    if (codec == "xnyzip")
      target["norm"] = "l2";
    report["error_bound_consistency"][name] = target;
  }
  for (bool pos : {true, false}) {
    std::string group = pos ? "positions" : "velocities";
    if (m["compressors"][group] != "xnyzip")
      continue;
    double max = 0, sum = 0, squared = 0;
    for (size_t row = 0; row < count; ++row) {
      double e2 = 0;
      for (size_t k = pos ? 1 : 4; k < (pos ? 4 : 7); ++k) {
        double d = (double(restored.data[k].number(row)) -
                    double(source.data[k].number(row))) /
                   (pos ? scale : 1);
        e2 += d * d;
      }
      double e = std::sqrt(e2);
      max = std::max(max, e);
      sum += e;
      squared += e2;
    }
    report["xnyzip_l2"][group] = {{"count", count},
                                  {"norm", "l2"},
                                  {"units", pos ? "lcp_units" : "source_units"},
                                  {"max_l2_error", max},
                                  {"mean_l2_error", sum / count},
                                  {"rms_l2_error", std::sqrt(squared / count)}};
    double bound = m["field_error_bounds"][group + "_xnyzip"]["abs"];
    auto b = m["field_error_bounds"][group + "_xnyzip"];
    report["xnyzip_l2_error_bound_consistency"][group] = {
        {"norm", "l2"},
        {"mode", b["mode"]},
        {"relative_error_bound", b["relative"]},
        {"range_for_relative", b["range"]},
        {"range_units", b["range_units"]},
        {"requested_l2_bound", bound},
        {"compressor_l2_bound", b["compressor_abs"]},
        {"preprocess_l2_allowance", b.value("preprocess_l2_max_abs", 0.)},
        {"observed_max_l2_error", max},
        {"satisfied", max <= bound + 1e-12 + 1e-6 * std::max(1., bound)}};
  }
  report["sizes"]["selected_particle_count"] = count;
  report["sizes"]["input_h5_file_bytes"] =
      m.value("input_h5_file_bytes", size_t(0));
  report["sizes"]["limit"] = m.value("limit", Json(nullptr));
  report["timing"]["metrics_wall_seconds"] =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return report;
}
} // namespace particle
