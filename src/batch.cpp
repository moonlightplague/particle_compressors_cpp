#include "particle/pipeline.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
namespace particle {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
std::vector<std::pair<int, std::string>> natural(const fs::path &p) {
  std::vector<std::pair<int, std::string>> v;
  std::istringstream s(p.filename().string());
  std::string part;
  while (std::getline(s, part, '.')) {
    bool digit = !part.empty() &&
                 std::all_of(part.begin(), part.end(),
                             [](unsigned char c) { return std::isdigit(c); });
    if (digit) {
      auto first = part.find_first_not_of('0');
      part = first == std::string::npos ? "0" : part.substr(first);
      part = std::string(30 - std::min<size_t>(30, part.size()), '0') + part;
    }
    v.emplace_back(digit ? 0 : 1, part);
  }
  return v;
}
std::vector<fs::path> discover(const fs::path &directory) {
  std::vector<fs::path> paths;
  for (const auto &e : fs::directory_iterator(directory)) {
    if (!e.is_regular_file())
      continue;
    auto name = e.path().filename().string();
    if (e.path().extension() == ".h5" ||
        (name.starts_with("dat_") && e.path().extension() != ".h5"))
      paths.push_back(e.path());
    if (name.starts_with("cfg_") &&
        !fs::exists(directory / ("dat_" + name.substr(4))) &&
        !fs::exists(directory / ("dat_" + name.substr(4) + ".h5"))) {
      std::ifstream f(e.path());
      std::string token;
      int64_t n = 0;
      for (int i = 0; i <= 9; ++i) {
        if (!(f >> token))
          throw std::runtime_error("Invalid native configuration");
        if (i == 9)
          n = std::stoll(token);
      }
      if (n)
        throw std::runtime_error(
            "Native configuration has no matching data file: " +
            e.path().string());
    }
  }
  std::sort(paths.begin(), paths.end(), [](const auto &a, const auto &b) {
    return natural(a) < natural(b);
  });
  if (paths.empty())
    throw std::runtime_error(
        "No HDF5 or native dat_* particle files found in directory: " +
        directory.string());
  return paths;
}
Json stable_attrs(Json a) {
  for (auto key : {"npart", "npart_total", "proc_size", "rank"})
    a.erase(key);
  return a;
}
int merge(const Options &o, const std::vector<fs::path> &files) {
  auto started = Clock::now();
  Snapshot joined;
  bool first = true;
  Json counts = Json::object();
  size_t bytes = 0;
  Json native = Json::array();
  for (const auto &path : files) {
    auto s = read_hdf5(path);
    if (!s.count())
      throw std::runtime_error("Cannot merge an empty snapshot");
    counts[path.filename().string()] = s.count();
    bytes += fs::file_size(path);
    if (s.metadata.contains("source"))
      native.push_back(s.metadata["source"]);
    if (first) {
      joined = std::move(s);
      first = false;
      continue;
    }
    if (json_text(stable_attrs(joined.metadata["root_attrs"])) !=
        json_text(stable_attrs(s.metadata["root_attrs"])))
      throw std::runtime_error("Root attributes differ between merge inputs");
    for (size_t i = 0; i < 7; ++i) {
      if (joined.data[i].type.name != s.data[i].type.name ||
          joined.metadata["fields"][fields[i]]["attrs"] !=
              s.metadata["fields"][fields[i]]["attrs"])
        throw std::runtime_error("Field schemas differ between merge inputs");
      joined.data[i].bytes.insert(joined.data[i].bytes.end(),
                                  s.data[i].bytes.begin(),
                                  s.data[i].bytes.end());
    }
  }
  auto order = stable_id_order(joined.data[0]);
  for (size_t i = 1; i < order.size(); ++i)
    if (joined.data[0].number(order[i - 1]) == joined.data[0].number(order[i]))
      throw std::runtime_error(
          "Input files contain overlapping or duplicate particle IDs");
  size_t n = joined.count();
  for (auto key : {"npart", "npart_total", "proc_size", "rank"}) {
    auto &attrs = joined.metadata["root_attrs"];
    if (!attrs.contains(key))
      continue;
    auto &a = attrs[key];
    if (!a["shape"].empty())
      throw std::runtime_error("Mutable merge attributes must be scalar");
    size_t value = std::string(key) == "rank"        ? 0
                   : std::string(key) == "proc_size" ? 1
                                                     : n;
    auto t = type_of(a["dtype"]);
    if (!t.floating && t.bytes < 8 &&
        value >= (uint64_t(1) << (t.bytes * 8 - (t.sign ? 1 : 0))))
      throw std::runtime_error("Merged count overflows root attribute dtype");
    a["value"] = value;
  }
  for (auto name : fields) {
    joined.metadata["fields"][name]["shape"] = {n};
    joined.metadata["fields"][name]["selected_shape"] = {n};
  }
  auto path = fs::path(o.work_dir) / "merged" / "merged.h5";
  require_output(path, o.force);
  fs::create_directories(path.parent_path());
  auto partial = path.parent_path() / ".merged.h5.partial";
  require_output(partial, o.force);
  try {
    write_hdf5(partial, joined, o.force);
    fs::rename(partial, path);
  } catch (...) {
    fs::remove(partial);
    throw;
  }
  Options next = o;
  next.input = path.string();
  next.merge = false;
  Json metadata = {
      {"enabled", true},
      {"input_directory", o.input},
      {"input_files", Json::array()},
      {"source_file_count", files.size()},
      {"source_particle_counts", counts},
      {"source_h5_file_bytes_total", bytes},
      {"merged_h5", path.string()},
      {"merged_h5_file_bytes", fs::file_size(path)},
      {"merged_particle_count", n},
      {"normalized_root_attributes",
       {"npart", "npart_total", "proc_size", "rank"}},
      {"id_overlap_check",
       {{"status", "passed"}, {"algorithm", "in_memory_sort"}, {"count", n}}},
      {"wall_seconds", seconds(started)}};
  for (auto &p : files)
    metadata["input_files"].push_back(p.string());
  if (!native.empty()) {
    metadata["native_sources"] = native;
    metadata["source_file_bytes_total"] = bytes;
    metadata.erase("source_h5_file_bytes_total");
  }
  next.merge_metadata = metadata;
  return run(next);
}
} // namespace
int run_directory(const Options &o) {
  auto paths = discover(o.input);
  if (o.merge)
    return merge(o, paths);
  fs::create_directories(o.work_dir);
  size_t workers = std::min<size_t>(
      paths.size(),
      o.file_workers
          ? o.file_workers
          : std::min(128u, std::max(1u, std::thread::hardware_concurrency())));
  std::map<pid_t, size_t> active;
  std::vector<Json> results(paths.size());
  size_t next = 0;
  auto started = Clock::now();
  auto launch = [&](size_t index) {
    Options child = o;
    child.input = paths[index].string();
    child.work_dir = (fs::path(o.work_dir) / paths[index].filename()).string();
    child.file_workers = 1;
    if (!child.field_workers)
      child.field_workers = std::min<size_t>(
          6,
          std::max<size_t>(1, std::thread::hardware_concurrency() / workers));
    std::cout.flush();
    std::cerr.flush();
    pid_t pid = fork();
    if (pid < 0)
      throw std::runtime_error("Cannot fork file worker");
    if (pid == 0) {
      auto t = Clock::now();
      int code = 0;
      Json r = {{"input_h5", child.input}, {"work_dir", child.work_dir}};
      try {
        code = run(child);
        r["report_path"] = (fs::path(child.work_dir) /
                            (child.metrics ? "metrics.json" : "manifest.json"))
                               .string();
        r["report"] =
            read_json(fs::path(child.work_dir) /
                      (child.metrics ? "metrics.json" : "manifest.json"));
      } catch (const H5::Exception &e) {
        r["error"] = e.getDetailMsg();
        code = 2;
      } catch (const std::exception &e) {
        r["error"] = e.what();
        code = 2;
      }
      r["wall_seconds"] = seconds(t);
      r["succeeded"] = code == 0;
      try {
        write_json(fs::path(child.work_dir) / ".batch-result.json", r);
      } catch (...) {
        code = 2;
      }
      std::cout.flush();
      std::cerr.flush();
      _exit(code);
    }
    active[pid] = index;
  };
  while (next < paths.size() || !active.empty()) {
    while (next < paths.size() && active.size() < workers)
      launch(next++);
    int status = 0;
    pid_t pid = waitpid(-1, &status, 0);
    if (pid < 0) {
      if (errno == EINTR)
        continue;
      throw std::runtime_error("Cannot wait for file worker");
    }
    if (!active.contains(pid))
      continue;
    size_t index = active.at(pid);
    active.erase(pid);
    auto path =
        fs::path(o.work_dir) / paths[index].filename() / ".batch-result.json";
    if (fs::exists(path)) {
      results[index] = read_json(path);
      fs::remove(path);
    } else
      results[index] = {
          {"input_h5", paths[index].string()},
          {"succeeded", false},
          {"error", "File worker terminated before writing its result"},
          {"wall_seconds", 0.}};
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
      results[index]["succeeded"] = false;
  }
  size_t success = 0, count = 0, original = 0, compressed = 0;
  bool complete = true;
  double total_seconds = 0;
  Json entries = Json::array(), stages = Json::object();
  std::vector<double> ratios, times;
  Json groups = Json::object();
  for (auto name : {"positions", "id", "velocities"})
    groups[name] = {{"original_bytes_total", 0}, {"compressed_bytes_total", 0}};
  bool groups_available = true;
  auto ratio = [](size_t original, size_t compressed) {
    return compressed ? double(original) / compressed
                      : std::numeric_limits<double>::infinity();
  };
  for (auto &r : results) {
    Json entry = r;
    entry.erase("report");
    entry.erase("succeeded");
    entry["status"] = r.value("succeeded", false) ? "success" : "failed";
    if (r.contains("report")) {
      auto &report = r["report"];
      auto sizes = report.value("sizes", Json::object());
      entry["particle_count"] = report.value(
          "count", sizes.value("selected_particle_count", size_t(0)));
      entry["selected_original_payload_bytes"] =
          sizes.value("selected_original_payload_bytes", size_t(0));
      entry["compressed_total_bytes"] =
          sizes.value("compressed_total_bytes", Json(nullptr));
      entry["payload_compression_ratio"] =
          sizes.value("payload_compression_ratio", Json(nullptr));
      entry["timing"] = report.value("timing", Json::object());
      Json quality = Json::object();
      auto report_fields = report.value("fields", Json::object());
      for (auto &[name, field] : report_fields.items()) {
        bool position = name == "x" || name == "y" || name == "z";
        auto display =
            position ? field.value("fixed_point_units", field) : field;
        if (!display.contains("max_absolute_error") ||
            !display.contains("mse") || !display.contains("psnr"))
          continue;
        quality[name] = {
            {"max_abs", display["max_absolute_error"]},
            {"mse", display["mse"]},
            {"psnr", display["psnr"]},
            {"units", position && !field.contains("fixed_point_units")
                          ? "lcp_units"
                          : "source_units"}};
      }
      if (!quality.empty())
        entry["quality_metrics"] = quality;
    }
    entries.push_back(entry);
    if (!r.value("succeeded", false)) {
      complete = false;
      continue;
    }
    ++success;
    total_seconds += r.value("wall_seconds", 0.);
    auto report = r.value("report", Json::object());
    count +=
        report.value("count", report.value("sizes", Json::object())
                                  .value("selected_particle_count", size_t(0)));
    auto sizes = report.value("sizes", Json::object());
    original += sizes.value("selected_original_payload_bytes", size_t(0));
    compressed += sizes.value("compressed_total_bytes", size_t(0));
    complete = complete && sizes.contains("compressed_total_bytes");
    times.push_back(r.value("wall_seconds", 0.));
    if (sizes.contains("payload_compression_ratio"))
      ratios.push_back(sizes["payload_compression_ratio"]);
    auto report_timing = report.value("timing", Json::object());
    for (auto &[name, value] : report_timing.items())
      if (value.is_number())
        stages[name] = stages.value(name, 0.) + value.get<double>();
    groups_available =
        groups_available && sizes.contains("compressed_components_bytes");
    for (size_t i = 0; i < fields.size(); ++i) {
      std::string group = i == 0 ? "id" : i < 4 ? "positions" : "velocities";
      auto field = report["fields"][fields[i]];
      auto dtype =
          field.value("dtype", field.value("original_dtype", std::string()));
      size_t n = report.value(
          "count", sizes.value("selected_particle_count", size_t(0)));
      groups[group]["original_bytes_total"] =
          groups[group]["original_bytes_total"].get<size_t>() +
          checked_bytes(n, type_of(dtype).bytes);
    }
    auto components =
        sizes.value("compressed_components_bytes", Json::object());
    for (auto &[path, value] : components.items()) {
      std::string group;
      for (auto prefix : {"compressed/positions.lcp",
                          "compressed/positions.xnyzip", "compressed/order.",
                          "compressed/x.", "compressed/y.", "compressed/z."})
        if (path.starts_with(prefix))
          group = "positions";
      for (auto prefix :
           {"compressed/velocities.lcp", "compressed/velocities.xnyzip",
            "compressed/velocity_order.", "compressed/velocity_block_ids.",
            "compressed/vx.", "compressed/vy.", "compressed/vz."})
        if (path.starts_with(prefix))
          group = "velocities";
      if (path.starts_with("compressed/id."))
        group = "id";
      if (!group.empty())
        groups[group]["compressed_bytes_total"] =
            groups[group]["compressed_bytes_total"].get<size_t>() +
            value.get<size_t>();
    }
  }
  complete = complete && success > 0;
  for (auto &group : groups) {
    if (complete && groups_available)
      group["compression_ratio"] =
          ratio(group["original_bytes_total"], group["compressed_bytes_total"]);
    else
      group = {{"original_bytes_total", nullptr},
               {"compressed_bytes_total", nullptr},
               {"compression_ratio", nullptr}};
  }
  auto statistics = [](std::vector<double> values) -> Json {
    if (values.empty())
      return Json::object();
    std::sort(values.begin(), values.end());
    size_t n = values.size();
    return {{"min", values.front()},
            {"max", values.back()},
            {"mean", std::accumulate(values.begin(), values.end(), 0.) / n},
            {"median",
             n % 2 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2}};
  };
  double wall = seconds(started);
  Json batch = {
      {"input_directory", o.input},
      {"command", o.command},
      {"workers", workers},
      {"summary",
       {{"discovered_files", paths.size()},
        {"successful_files", success},
        {"failed_files", paths.size() - success},
        {"total_particle_count", count}}},
      {"sizes",
       {{"selected_original_payload_bytes_total", original},
        {"compressed_total_bytes", complete ? Json(compressed) : Json(nullptr)},
        {"payload_compression_ratio",
         complete ? Json(ratio(original, compressed)) : Json(nullptr)},
        {"field_groups", groups}}},
      {"timing",
       {{"batch_wall_seconds", wall},
        {"file_wall_seconds_total", total_seconds},
        {"effective_parallel_speedup", total_seconds / wall},
        {"payload_throughput_mib_per_second",
         double(original) / (1024 * 1024) / wall},
        {"stage_seconds_total", stages}}},
      {"statistics",
       {{"per_file_payload_compression_ratio", statistics(ratios)},
        {"per_file_wall_seconds", statistics(times)}}},
      {"files", entries}};
  write_json(fs::path(o.work_dir) / "batch_metrics.json", batch);
  if (success != paths.size())
    throw std::runtime_error(std::to_string(paths.size() - success) +
                             " file pipelines failed; see batch_metrics.json");
  std::cout << "batch_metrics = "
            << (fs::path(o.work_dir) / "batch_metrics.json").string() << '\n';
  return 0;
}
} // namespace particle
