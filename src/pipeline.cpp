#include "particle/pipeline.hpp"
#include "particle/huffman.hpp"
#include "particle/layout.hpp"
#include "particle/runtime.hpp"
#include "particle/triplets.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace particle {
namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}
std::string extension(const std::string &c) {
  return c == "pcodec" ? "pco" : c == "sz3" ? "psz" : c;
}
Array read_array(const std::string &path, const Type &t, size_t count) {
  auto b = read_bytes(path);
  if (b.size() != checked_bytes(count, t.bytes))
    throw std::runtime_error("Raw field size mismatch: " + path);
  return {t, std::move(b)};
}
Json selection(const Options &o, bool position, double range) {
  std::string codec = position ? o.pos_compressor : o.vel_compressor;
  Json result = {{"mode", "lossless"},
                 {"abs", 0.0},
                 {"relative", nullptr},
                 {"range", range},
                 {"range_units", position ? "lcp_units" : "source_units"},
                 {"compressor_abs", 0.0}};
  if (codec == "pcodec")
    return result;
  auto absolute = position ? o.pos_abs_eb : o.vel_abs_eb;
  auto relative = position ? o.pos_rel_eb : o.vel_rel_eb;
  if (absolute && relative)
    throw std::runtime_error(
        "Field-specific absolute and relative error bounds cannot both be set");
  if (!absolute && !relative) {
    absolute = o.abs_eb;
    relative = o.rel_eb;
  }
  double bound;
  if (relative) {
    result["mode"] = "relative";
    result["relative"] = *relative;
    bound = *relative * range;
  } else if (absolute) {
    result["mode"] = "absolute";
    bound = *absolute;
  } else
    throw std::runtime_error("No error bound selected");
  if (bound < 0 || !std::isfinite(bound))
    throw std::runtime_error("Error bounds must be finite and non-negative");
  result["abs"] = bound;
  result["compressor_abs"] = bound;
  return result;
}

} // namespace
Json preprocess(const Options &o) {
  validate(o);
  auto started = Clock::now();
  fs::path work = o.work_dir;
  auto source = read_hdf5(o.input, o.limit);
  size_t count = source.count();
  if (!count)
    throw std::runtime_error("Not a valid input array.");
  fs::create_directories(work);
  require_output(work / "manifest.json", o.force);
  double scale = 1;
  std::string scale_mode = o.position_scale;
  Json attr = nullptr;
  auto attrs = source.metadata["root_attrs"];
  if (o.position_scale == "value") {
    if (!o.position_scale_value)
      throw std::runtime_error(
          "--position-scale value requires --position-scale-value");
    scale = *o.position_scale_value;
  } else if (o.position_scale == "attr" ||
             (o.position_scale == "auto" && !source.data[1].type.floating &&
              attrs.contains(o.position_scale_attr))) {
    if (!attrs.contains(o.position_scale_attr))
      throw std::runtime_error("Position scale attribute missing");
    auto v = attrs[o.position_scale_attr]["value"];
    while (v.is_array() && v.size() == 1)
      v = v[0];
    scale = v.get<double>();
    attr = o.position_scale_attr;
    if (o.position_scale == "auto") {
      if (scale > 0)
        scale_mode = "auto_attr";
      else {
        scale = 1;
        scale_mode = "auto_raw";
        attr = nullptr;
      }
    }
  } else if (o.position_scale == "auto")
    scale_mode = "auto_raw";
  else if (o.position_scale != "raw")
    throw std::runtime_error("Invalid position scale mode");
  if (!(scale > 0) || !std::isfinite(scale))
    throw std::runtime_error("Position scale must be finite and positive");
  Json m = {{"format_version", 2},
            {"input_h5", o.input},
            {"input_h5_file_bytes", fs::file_size(o.input)},
            {"count", count},
            {"limit", o.limit ? Json(o.limit) : Json(nullptr)},
            {"fields", source.metadata["fields"]},
            {"root_attrs", attrs},
            {"position_scale",
             {{"mode", scale_mode}, {"value", scale}, {"attr", attr}}},
            {"compressors",
             {{"positions", o.pos_compressor},
              {"velocities", o.vel_compressor},
              {"lossless", "pcodec"}}},
            {"tools",
             {{"lcp", o.lcp},
              {"pcodec", nullptr},
              {"pysz", nullptr},
              {"pyszo", nullptr}}},
            {"velocity_chunking",
             {{"chunk_size", o.vel_chunk_size},
              {"enabled", o.vel_chunk_size != 0},
              {"configured_workers", o.vel_chunk_workers}}},
            {"blockwise_order", {{"enabled", false}, {"field", nullptr}}},
            {"order_dtype", "int32"},
            {"compressed_fields", Json::object()}};
  if (!o.merge_metadata.empty()) {
    m["merge"] = o.merge_metadata;
    m["timing"]["merge_wall_seconds"] = o.merge_metadata.at("wall_seconds");
  }
  if (source.metadata.contains("source")) {
    m["source"] = source.metadata["source"];
    m["input_format"] = "cfg_dat_field_major_v1";
    m["input_file"] = o.input;
    m["input_file_bytes"] = fs::file_size(o.input);
    m.erase("input_h5_file_bytes");
    fs::create_directories(work / "input_adapters");
    auto adapter = work / "input_adapters" /
                   (fs::path(o.input).filename().string() + ".h5");
    write_hdf5(adapter, source, true);
    m["input_h5"] = adapter.string();
  }
  Json raw = Json::object(), compressed = Json::object(),
       bounds = Json::object(), stats = Json::object();
  size_t payload = 0;
  double pos_min_bound = std::numeric_limits<double>::infinity(),
         vel_max_bound = 0, pos_diagonal = 0, pos_cast_squared = 0;
  double vector_requested = 0;
  std::array<Array, 7> exported = source.data;
  for (size_t i = 0; i < 7; ++i) {
    std::string name = fields[i];
    auto &a = source.data[i];
    bool pos = i > 0 && i < 4;
    std::string codec = i == 0 ? "pcodec"
                        : pos  ? o.pos_compressor
                               : o.vel_compressor;
    payload += a.bytes.size();
    double lo = std::numeric_limits<double>::infinity(), hi = -lo, cast = 0;
    long double rawlo = a.number(0), rawhi = rawlo;
    if (pos && codec != "pcodec")
      exported[i] = {type_of("float32"),
                     std::vector<uint8_t>(checked_bytes(count, 4))};
    for (size_t j = 0; j < count; ++j) {
      long double v = a.number(j);
      rawlo = std::min(rawlo, v);
      rawhi = std::max(rawhi, v);
      double x = static_cast<double>(v) / (pos ? scale : 1);
      if (i >= 4 && (codec == "lcp" || codec == "xnyzip"))
        cast = std::max(cast, std::abs(double(float(x)) - x));
      if (!std::isfinite(x))
        throw std::runtime_error(
            "Nonfinite input is not supported by preprocessing");
      lo = std::min(lo, x);
      hi = std::max(hi, x);
      if (pos && codec != "pcodec") {
        float f = static_cast<float>(x);
        exported[i].set(j, f);
        cast = std::max(cast, std::abs(static_cast<double>(f) - x));
      }
    }
    auto filename =
        name + "." + (pos && codec != "pcodec" ? "f32" : a.type.name) + ".raw";
    raw[name] = (work / "preprocessed" / filename).string();
    compressed[name] =
        (work / "compressed" / (name + "." + extension(codec))).string();
    if (i == 0) {
      Json min = a.type.sign ? Json(static_cast<int64_t>(rawlo))
                             : Json(static_cast<uint64_t>(rawlo));
      Json max = a.type.sign ? Json(static_cast<int64_t>(rawhi))
                             : Json(static_cast<uint64_t>(rawhi));
      stats["id"] = {{"source_dtype", a.type.name},
                     {"min", min},
                     {"max", max},
                     {"pcodec_dtype", a.type.name}};
      bounds[name] = {{"mode", "lossless"},
                      {"abs", o.id_abs_eb},
                      {"relative", nullptr},
                      {"range", static_cast<double>(rawhi - rawlo)},
                      {"range_units", "source_units"},
                      {"compressor_abs", 0.0}};
      continue;
    }
    bounds[name] = selection(o, pos, hi - lo);
    double requested = bounds[name]["abs"];
    double prep = pos ? cast + (!a.type.floating ? .5 / scale : 0) : 0;
    double cb = requested;
    if (pos && bounds[name]["mode"] == "relative")
      cb = std::max(0., requested - prep);
    bounds[name]["compressor_abs"] = cb;
    if (pos) {
      stats["positions"][name] = {
          {"scale", scale},
          {"preprocess_cast_max_abs_in_lcp_units", cast},
          {"preprocess_cast_max_abs_in_original_fixed_point_units",
           cast * scale},
          {"min_in_lcp_units", lo},
          {"max_in_lcp_units", hi},
          {"range_in_lcp_units", hi - lo},
          {"min_in_int_units", static_cast<double>(std::trunc(rawlo))},
          {"max_in_int_units", static_cast<double>(std::trunc(rawhi))},
          {"range_in_int_units",
           static_cast<double>(std::trunc(rawhi) - std::trunc(rawlo))}};
      pos_min_bound = std::min(pos_min_bound, cb);
      pos_diagonal += (hi - lo) * (hi - lo);
      pos_cast_squared += prep * prep;
      vector_requested = requested;
    } else {
      stats["velocities"][name] = {{"float_min", lo},
                                   {"float_max", hi},
                                   {"float_range", hi - lo},
                                   {"preprocess_cast_max_abs", cast}};
      vel_max_bound = std::max(vel_max_bound, requested);
    }
  }
  double diagonal = std::sqrt(pos_diagonal);
  bool relative = bounds["x"]["mode"] == "relative";
  if (relative)
    vector_requested = bounds["x"]["relative"].get<double>() * diagonal;
  double vector_preprocess = (relative || o.pos_compressor == "xnyzip")
                                 ? std::sqrt(pos_cast_squared)
                                 : 0;
  double vector_bound = std::max(0., vector_requested - vector_preprocess);
  bounds["positions_xnyzip"] = {{"mode", bounds["x"]["mode"]},
                                {"abs", vector_requested},
                                {"relative", bounds["x"]["relative"]},
                                {"range", diagonal},
                                {"range_units", "lcp_units_bbox_diagonal"},
                                {"compressor_abs", vector_bound},
                                {"preprocess_l2_max_abs", vector_preprocess}};
  m["error_bounds"] = {{"positions_lcp_abs", pos_min_bound},
                       {"positions_xnyzip_abs", vector_bound},
                       {"velocities_sz3_abs", vel_max_bound},
                       {"id_sz3_abs", o.id_abs_eb}};
  if (o.pos_compressor == "xnyzip" && vector_bound <= 0)
    throw std::runtime_error(
        "No positive XnYZip position budget after preprocessing");
  if (o.vel_compressor == "lcp" || o.vel_compressor == "xnyzip") {
    double minimum = std::numeric_limits<double>::infinity(), diagonal2 = 0,
           cast2 = 0;
    for (size_t i = 4; i < 7; ++i) {
      auto &b = bounds[fields[i]];
      double requested = b["abs"], range = b["range"],
             cast = stats["velocities"][fields[i]]["preprocess_cast_max_abs"];
      minimum = std::min(minimum, std::max(0., requested - cast));
      diagonal2 += range * range;
      cast2 += cast * cast;
    }
    if (o.vel_compressor == "lcp") {
      for (size_t i = 4; i < 7; ++i)
        bounds[fields[i]]["compressor_abs"] = minimum;
      m["error_bounds"]["velocities_lcp_abs"] = minimum;
    } else {
      double requested =
          bounds["vx"]["mode"] == "relative"
              ? bounds["vx"]["relative"].get<double>() * std::sqrt(diagonal2)
              : bounds["vx"]["abs"].get<double>();
      double remaining = std::max(0., requested - std::sqrt(cast2));
      if (remaining <= 0)
        throw std::runtime_error(
            "No positive XnYZip velocity budget after preprocessing");
      bounds["velocities_xnyzip"] = {
          {"mode", bounds["vx"]["mode"]},
          {"abs", requested},
          {"relative", bounds["vx"]["relative"]},
          {"range", std::sqrt(diagonal2)},
          {"range_units", "source_units_bbox_diagonal"},
          {"compressor_abs", remaining},
          {"preprocess_l2_max_abs", std::sqrt(cast2)}};
      m["error_bounds"]["velocities_xnyzip_abs"] = remaining;
    }
  }
  for (bool position : {true, false}) {
    auto codec = position ? o.pos_compressor : o.vel_compressor;
    if (codec == "lcp" || codec == "xnyzip") {
      size_t first = position ? 1 : 4;
      for (size_t i = first; i < first + 3; ++i)
        compressed.erase(fields[i]);
      compressed[position ? "positions" : "velocities"] =
          (work / "compressed" /
           ((position ? std::string("positions.")
                      : std::string("velocities.")) +
            codec))
              .string();
      if (!position)
        compressed["velocity_order"] =
            (work / "compressed" / "velocity_order.pco").string();
    }
  }
  m["field_error_bounds"] = bounds;
  m["preprocess"] = stats;
  m["sizes"] = {{"selected_original_payload_bytes", payload}};
  m["artifacts"] = {{"preprocessed", raw}, {"compressed", compressed}};
  // Validate every output before replacing an existing package.
  for (auto &p : raw)
    require_output(p.get<std::string>(), o.force);
  if (o.force && fs::exists(work / "compressed"))
    fs::remove_all(work / "compressed");
  fs::create_directories(work / "compressed");
  for (size_t i = 0; i < 7; ++i)
    write_bytes(raw[fields[i]].get<std::string>(), exported[i].bytes, o.force);
  m["timing"]["preprocess_wall_seconds"] = elapsed(started);
  write_json(work / "manifest.json", m);
  return m;
}
void update_sizes(Json &m, const fs::path &work) {
  Json components = Json::object();
  if (fs::exists(work / "compressed"))
    for (auto &p : fs::recursive_directory_iterator(work / "compressed"))
      if (p.is_regular_file())
        components[fs::relative(p.path(), work).string()] = p.file_size();
  components["manifest.json"] = 0;
  for (int i = 0; i < 10; ++i) {
    size_t total = 0;
    for (auto &v : components)
      total += v.get<size_t>();
    auto &s = m["sizes"];
    s["compressed_components_bytes"] = components;
    s["compressed_total_bytes"] = total;
    s["payload_compression_ratio"] =
        total
            ? double(s["selected_original_payload_bytes"].get<size_t>()) / total
            : 0;
    if (m.contains("input_h5_file_bytes"))
      s["h5_file_to_compressed_ratio"] =
          total ? double(m["input_h5_file_bytes"].get<size_t>()) / total : 0;
    size_t size = json_text(m).size();
    if (components["manifest.json"] == size)
      break;
    components["manifest.json"] = size;
  }
}
void compress(const Options &o, Json &m) {
  auto started = Clock::now();
  size_t count = m.at("count");
  std::vector<size_t> order;
  auto &raw = m["artifacts"]["preprocessed"];
  std::array<Array, 7> data;
  for (size_t i = 0; i < 7; ++i) {
    auto t = type_of(i > 0 && i < 4 && o.pos_compressor != "pcodec"
                         ? "float32"
                         : m["fields"][fields[i]]["dtype"].get<std::string>());
    data[i] = read_array(raw[fields[i]], t, count);
  }
  auto is_triplet = [](const std::string &c) {
    return c == "lcp" || c == "xnyzip";
  };
  auto save_triplet = [&](bool position, TripletEncoded &c) {
    std::string group = position ? "positions" : "velocities";
    auto path = m["artifacts"]["compressed"][group].get<std::string>();
    write_bytes(path, c.bytes, o.force);
    auto f = c.metadata;
    f["field"] = group;
    f["path"] = path;
    f["bytes"] = c.bytes.size();
    size_t first = position ? 1 : 4;
    for (size_t i = first; i < first + 3; ++i)
      f["source_dtypes"][fields[i]] = m["fields"][fields[i]]["dtype"];
    if (f["codec"] == "xnyzip") {
      f["interleaved_fields"] = {fields[first], fields[first + 1],
                                 fields[first + 2]};
      m["error_bounds"][group + "_xnyzip_abs"] = f["l2_error_bound"];
      m["field_error_bounds"][group + "_xnyzip"]["compressor_abs"] =
          f["l2_error_bound"];
    }
    m["compressed_fields"][group] = f;
  };
  auto canonical_started = Clock::now();
  std::string mapping = "original_row";
  std::optional<Lattice> lattice;
  Json structured;
  Triplet canonical_positions;
  std::vector<size_t> hybrid;
  bool lattice_possible =
      o.lattice_layout &&
      (o.vel_compressor == "sz3" || o.vel_compressor == "szo") &&
      (is_triplet(o.pos_compressor) || o.pos_compressor == "sz3" ||
       o.pos_compressor == "szo");
  if (o.xnyzip_structure_aware) {
    m["structured_layout"] = {{"requested", true}, {"enabled", false}};
    try {
      if (o.pos_compressor != "xnyzip" ||
          (o.vel_compressor != "szo" && o.vel_compressor != "xnyzip"))
        throw std::runtime_error(
            "Requires XnYZip positions and SZO or XnYZip velocities");
      if (!m["root_attrs"].contains("nsidemesh"))
        throw std::runtime_error("Root attribute nsidemesh is unavailable");
      auto geometry = infer_lattice(
          data[0], {data[1], data[2], data[3]},
          m["root_attrs"]["nsidemesh"]["value"].get<size_t>(), 0, true);
      structured = structured_metadata(geometry, o.xnyzip_velocity_cell_bits);
      structured["requested"] = true;
      m["structured_layout"] = structured;
    } catch (const std::runtime_error &e) {
      m["structured_layout"]["reason"] = e.what();
    }
  }

  if (is_triplet(o.pos_compressor)) {
    double b =
        m["error_bounds"][o.pos_compressor == "lcp" ? "positions_lcp_abs"
                                                    : "positions_xnyzip_abs"];
    auto c = encode_triplet({data[1], data[2], data[3]}, o.pos_compressor, b, 0,
                            !structured.is_null());
    order = c.order;
    if (!structured.is_null())
      canonical_positions = decode_triplet(c.bytes, c.metadata);
    save_triplet(true, c);
    mapping = o.pos_compressor + "_position_sorted";
  } else if (o.sort || lattice_possible) {
    order = stable_id_order(data[0]);
    mapping = "id_sorted";
  }
  if (!order.empty())
    for (auto &a : data)
      a = a.gather(order);
  m["timing"]["canonical_order_wall_seconds"] = elapsed(canonical_started);
  auto lattice_started = Clock::now();
  if (o.lattice_layout) {
    m["lattice_layout"] = {{"requested", true},
                           {"enabled", false},
                           {"minimum_occupancy", o.lattice_min_occupancy},
                           {"axis_search", o.lattice_axis_search},
                           {"field_scope", is_triplet(o.pos_compressor)
                                               ? "velocities"
                                               : "positions_and_velocities"}};
    try {
      if (!lattice_possible)
        throw std::runtime_error("Lattice layout requires fieldwise SZ3/SZO "
                                 "velocities and compatible positions");
      if (!m["root_attrs"].contains("nsidemesh"))
        throw std::runtime_error("Root attribute nsidemesh is unavailable");
      lattice =
          infer_lattice(data[0], {data[1], data[2], data[3]},
                        m["root_attrs"]["nsidemesh"]["value"].get<size_t>(),
                        o.lattice_min_occupancy);
      m["lattice_layout"].update(lattice->metadata());
    } catch (const std::runtime_error &e) {
      m["lattice_layout"]["reason"] = e.what();
    }
  }
  if (!structured.is_null()) {
    hybrid = hybrid_order(data[0], canonical_positions, structured);
    for (size_t i = 4; i < 7; ++i)
      data[i] = data[i].gather(hybrid);
  }

  m["timing"]["lattice_prepare_wall_seconds"] = elapsed(lattice_started);
  bool reordered = !order.empty();
  m["ordering"]["reconstructed_rows"] = {
      {"mapping", mapping},
      {"original_row_order_restored", !reordered},
      {"canonical_field", is_triplet(o.pos_compressor) ? Json("positions")
                          : mapping == "id_sorted"     ? Json("id")
                                                       : Json(nullptr)},
      {"canonical_lcp_field",
       o.pos_compressor == "lcp" ? Json("positions") : Json(nullptr)},
      {"lcp_permutation_stored", false},
      {"temporary_permutation_artifact", nullptr},
      {"temporary_permutation_dtype", nullptr}};
  if (reordered) {
    auto dtype = type_of(mapping == "id_sorted"      ? "int64"
                         : o.pos_compressor == "lcp" ? "int32"
                                                     : "uint64");
    std::string key =
        mapping == "id_sorted" ? "id_sort_order" : "position_order";
    auto path = fs::path(o.work_dir) / "preprocessed" /
                (key + "." + dtype.name + ".raw");
    Array indices{dtype,
                  std::vector<uint8_t>(checked_bytes(count, dtype.bytes))};
    for (size_t i = 0; i < count; ++i)
      indices.set(i, order[i]);
    write_bytes(path, indices.bytes, o.force);
    m["artifacts"]["preprocessed"][key] = path.string();
    m["ordering"]["reconstructed_rows"]["temporary_permutation_artifact"] = key;
    m["ordering"]["reconstructed_rows"]["temporary_permutation_dtype"] =
        dtype.name;
  }
  m["particle_sort"] = {
      {"requested", o.sort || o.lattice_layout},
      {"stable", mapping == "id_sorted"},
      {"enabled", mapping == "id_sorted"},
      {"key", mapping == "id_sorted" ? Json("id") : Json(nullptr)},
      {"direction",
       mapping == "id_sorted" ? Json("ascending") : Json(nullptr)}};
  auto field_payload = [&](const std::string &name, const Array &a,
                           const std::string &codec, double bound,
                           unsigned level = 8) {
    size_t encoded = 0;
    auto bytes = encode(a, codec, bound, encoded, level);
    auto path = m["artifacts"]["compressed"][name].get<std::string>();
    write_bytes(path, bytes, o.force);
    Json f = {{"field", name},
              {"codec", codec == "sz3"           ? "pysz"
                        : codec == "szo_lorenzo" ? "szo"
                                                 : codec},
              {"dtype", a.type.name},
              {"count", a.size()},
              {"path", path},
              {"bytes", bytes.size()}};
    if (codec != "pcodec") {
      f["abs_error_bound"] = bound;
      f["encoded_count"] = encoded;
    }
    return f;
  };
  auto save_field = [&](const std::string &name, const Array &a,
                        const std::string &codec, double bound,
                        unsigned level = 8) {
    m["compressed_fields"][name] = field_payload(name, a, codec, bound, level);
  };
  auto id_started = Clock::now();
  if (!structured.is_null()) {
    auto codes = hilbert_ids(data[0], structured, false, data[0].type);
    save_field("id", codes, "pcodec", 0, 12);
    auto &f = m["compressed_fields"]["id"];
    f.update({{"codec", hilbert_codec},
              {"dtype", data[0].type.name},
              {"encoded_dtype", "uint32"},
              {"transform", "physical_axis_hilbert_3d"},
              {"lossless_backend", "pcodec"},
              {"pcodec_compression_level", 12},
              {"structured_layout", structured}});
  } else
    save_field("id", data[0], "pcodec", 0);
  m["timing"]["id_compress_wall_seconds"] = elapsed(id_started);
  auto fields_started = Clock::now();
  std::vector<size_t> jobs;
  for (size_t i = 1; i < 7; ++i)
    if (!(i < 4 ? is_triplet(o.pos_compressor) : is_triplet(o.vel_compressor)))
      jobs.push_back(i);
  size_t workers = worker_count(o.field_workers, jobs.size());
  m["runtime"]["compression_field_workers"] = workers;
  auto results = parallel_map(jobs.size(), workers, [&](size_t index) {
    size_t i = jobs[index];
    bool pos = i < 4;
    std::string name = fields[i],
                codec = pos ? o.pos_compressor : o.vel_compressor;
    if (lattice) {
      Json local = m;
      return compress_lattice_field(o, local, name, data[i], *lattice,
                                    pos ? i - 1 : i - 4, pos);
    }
    bool hybrid_field = i >= 4 && !structured.is_null();
    auto f = field_payload(
        name, data[i], hybrid_field && codec == "szo" ? "szo_lorenzo" : codec,
        m["field_error_bounds"][name]["compressor_abs"]);
    if (hybrid_field) {
      f["spatial_layout"] = hybrid_name;
      f["compression_profile"] = "lorenzo_1d";
    }
    return f;
  });
  for (size_t j = 0; j < jobs.size(); ++j) {
    auto &f = results[j];
    std::string name = fields[jobs[j]];
    m["compressed_fields"][name] = f;
    if (f.contains("lattice_wrap_field"))
      m["artifacts"]["compressed"][name + "_lattice_wrap"] =
          f["lattice_wrap_field"]["path"];
  }
  if (is_triplet(o.vel_compressor)) {
    double b =
        m["error_bounds"][o.vel_compressor == "lcp" ? "velocities_lcp_abs"
                                                    : "velocities_xnyzip_abs"];
    auto c = o.blockwise_ord
                 ? encode_blockwise({data[4], data[5], data[6]}, b)
                 : encode_triplet({data[4], data[5], data[6]}, o.vel_compressor,
                                  b, o.vel_chunk_size, false,
                                  o.vel_chunk_workers, false);
    save_triplet(false, c);
    if (!structured.is_null())
      m["compressed_fields"]["velocities"]["spatial_layout"] = hybrid_name;
    if (o.blockwise_ord) {
      save_field("velocity_order", c.packed_order, "pcodec", 0);
      auto &f = m["compressed_fields"]["velocity_order"];
      f.update({{"order_encoding", "lcp_blockwise_packed"},
                {"word_bits", 32},
                {"packed_word_count", c.packed_order.size()},
                {"particle_count", count},
                {"uncompressed_bytes", c.packed_order.bytes.size()},
                {"index_scope", "block_local_packed"},
                {"chunk_size", 0},
                {"chunk_count", 1},
                {"applied_during_lcp_decompression", true},
                {"block_id_field", "velocity_block_ids"}});
      Json metadata;
      auto bytes = huffman_encode(c.block_ids, metadata);
      Array huffman{type_of("uint8"), std::move(bytes)};
      m["artifacts"]["compressed"]["velocity_block_ids"] =
          (fs::path(o.work_dir) / "compressed" / "velocity_block_ids.pco")
              .string();
      save_field("velocity_block_ids", huffman, "pcodec", 0);
      m["compressed_fields"]["velocity_block_ids"].update(
          {{"preprocessor", metadata},
           {"decoded_dtype", c.block_ids.type.name},
           {"decoded_count", count},
           {"decoded_bytes", c.block_ids.bytes.size()},
           {"applied_during_lcp_decompression", true}});
      m["blockwise_order"] = {{"enabled", true}, {"field", "velocities"}};
    } else {
      auto type = type_of(o.vel_compressor == "lcp" ? "int32" : "uint64");
      Array indices{type,
                    std::vector<uint8_t>(checked_bytes(count, type.bytes))};
      for (size_t i = 0; i < count; ++i)
        indices.set(i, c.order.at(i));
      save_field("velocity_order", indices, "pcodec", 0);
      m["compressed_fields"]["velocity_order"]["chunk_size"] = o.vel_chunk_size;
    }
  }
  m["timing"]["lossy_fields_wall_seconds"] = elapsed(fields_started);
  for (auto group : {"id", "positions", "velocities"})
    m["ordering"][group] = {{"mapping", mapping}};
  if (is_triplet(o.pos_compressor)) {
    m["ordering"]["reconstructed_rows"]["position_permutation_stored"] = false;
    m["ordering"]["id"]["replaces_" + o.pos_compressor + "_position_order"] =
        true;
    if (o.pos_compressor == "xnyzip") {
      m["ordering"]["reconstructed_rows"]["canonical_xnyzip_field"] =
          "positions";
      m["ordering"]["reconstructed_rows"]["xnyzip_permutation_stored"] = false;
    }
  }
  if (lattice) {
    if (!is_triplet(o.pos_compressor))
      m["ordering"]["positions"]["spatial_layout"] = lattice_name;
    m["ordering"]["velocities"]["spatial_layout"] = lattice_name;
  }
  if (!structured.is_null())
    m["ordering"]["velocities"] = {{"mapping", hybrid_name},
                                   {"reconstructed_mapping", mapping}};
  if (is_triplet(o.vel_compressor)) {
    std::string velocity_mapping =
        o.vel_compressor == "lcp"
            ? "lcp_velocity_sorted_index_to_lcp_position_sorted_row"
            : "xnyzip_velocity_sorted_index_to_" +
                  (structured.is_null() ? mapping : std::string(hybrid_name)) +
                  "_row";
    auto scope = o.blockwise_ord    ? "block_local_packed"
                 : o.vel_chunk_size ? "chunk_local"
                                    : "global";
    auto &order_metadata = m["ordering"]["velocities"];
    order_metadata = {{"mapping", velocity_mapping},
                      {"field", "velocity_order"},
                      {"index_scope", scope},
                      {"chunk_size", o.vel_chunk_size}};
    if (o.vel_compressor == "lcp")
      order_metadata["applied_during_lcp_decompression"] = o.blockwise_ord;
    if (o.blockwise_ord)
      order_metadata["block_id_field"] = "velocity_block_ids";
    if (!structured.is_null())
      order_metadata["reconstructed_mapping"] = mapping;
    auto &sidecar = m["compressed_fields"]["velocity_order"];
    sidecar["order_mapping"] = velocity_mapping;
    sidecar["index_scope"] = scope;
  }
  m["order_dtype"] = o.pos_compressor == "xnyzip" ? "uint64" : "int32";
  m["format_version"] =
      o.vel_compressor == "xnyzip" && o.vel_chunk_size                 ? 6
      : (o.pos_compressor == "xnyzip" || o.vel_compressor == "xnyzip") ? 5
      : o.vel_chunk_size                                               ? 4
                                                                       : 3;
  if (o.blockwise_ord)
    m["format_version"] = 7;
  if (lattice)
    m["format_version"] = 8;
  if (!structured.is_null())
    m["format_version"] = 9;
  size_t chunks = o.vel_chunk_size ? count / o.vel_chunk_size +
                                         (count % o.vel_chunk_size != 0)
                                   : 1;
  m["velocity_chunking"]["chunk_count"] = chunks;
  m["velocity_chunking"]["effective_workers"] =
      o.vel_chunk_size ? worker_count(o.vel_chunk_workers,
                                      std::numeric_limits<size_t>::max(), 16)
                       : 1;
  m["timing"]["compress_wall_seconds"] = elapsed(started);
  update_sizes(m, o.work_dir);
  write_json(fs::path(o.work_dir) / "manifest.json", m);
}
void decompress(const Options &o, Json &m) {
  auto started = Clock::now();
  size_t count = m.at("count");
  auto work = fs::path(o.work_dir);
  require_output(work / "reconstructed.h5", o.force);
  if (!m.contains("compressed_fields") || m["compressed_fields"].empty())
    throw std::runtime_error("Manifest does not contain compressed_fields");
  // Older LCP manifests identify the position stream only in artifacts.
  if (!m["compressed_fields"].contains("positions") &&
      !m["compressed_fields"].contains("x") &&
      m["artifacts"]["compressed"].contains("positions")) {
    auto path = m["artifacts"]["compressed"]["positions"].get<std::string>();
    if (!path.ends_with(".xnyzip"))
      m["compressed_fields"]["positions"] = {{"field", "positions"},
                                             {"codec", "lcp"},
                                             {"dtype", "float32"},
                                             {"count", count},
                                             {"path", path}};
  }
  if (m["compressed_fields"]
              .value("velocities", Json::object())
              .value("codec", "") == "xnyzip" &&
      !m["compressed_fields"].contains("velocity_order"))
    throw std::runtime_error(
        "XnYZip velocity package is missing velocity_order");
  Snapshot s;
  s.metadata = {{"fields", m["fields"]},
                {"root_attrs", m.value("root_attrs", Json::object())}};
  Json paths = Json::object();
  std::array<Array, 7> decoded;
  std::function<Array(const Json &)> load_field = [&](const Json &f) -> Array {
    size_t n = f.at("count");
    std::string codec = f.value("codec", "pcodec");
    auto bytes = read_bytes(f.at("path").get<std::string>());
    if (codec == hilbert_codec) {
      if (!m.value("structured_layout", Json::object()).value("enabled", false))
        throw std::runtime_error(
            "Structured ID field lacks enabled layout metadata");
      auto codes = decode(bytes, type_of(f.value("encoded_dtype", "uint32")), n,
                          n, "pcodec");
      return hilbert_ids(codes, f.at("structured_layout"), true,
                         type_of(f.at("dtype")));
    }
    if (f.value("spatial_layout", "") == lattice_name) {
      auto dense = decode_shaped(bytes, f);
      auto l = load_lattice(decoded[0], m.at("lattice_layout"));
      std::string name = f.at("field");
      size_t axis = name == "x" || name == "vx"   ? 0
                    : name == "y" || name == "vy" ? 1
                                                  : 2;
      std::optional<Array> wrap;
      if (f.contains("lattice_wrap_field"))
        wrap = load_field(f["lattice_wrap_field"]);
      return lattice_decode(dense, l, axis,
                            f.value("lattice_transform", "identity"),
                            wrap ? &*wrap : nullptr);
    }
    return decode(bytes, type_of(f.at("dtype")), n, f.value("encoded_count", n),
                  codec);
  };
  decoded[0] = load_field(m["compressed_fields"].at("id"));
  std::vector<size_t> field_jobs;
  for (size_t i = 1; i < 7; ++i)
    if (m["compressed_fields"].contains(fields[i]))
      field_jobs.push_back(i);
  auto field_started = Clock::now();
  m["runtime"]["decompression_field_workers"] =
      worker_count(o.field_workers, field_jobs.size());
  auto field_results = parallel_map(
      field_jobs.size(), worker_count(o.field_workers, field_jobs.size()),
      [&](size_t j) {
        return load_field(m["compressed_fields"].at(fields[field_jobs[j]]));
      });
  m["timing"]["lossy_fields_decompress_wall_seconds"] = elapsed(field_started);
  for (size_t j = 0; j < field_jobs.size(); ++j)
    decoded[field_jobs[j]] = std::move(field_results[j]);
  for (bool pos : {true, false}) {
    std::string group = pos ? "positions" : "velocities";
    if (m["compressed_fields"].contains(group)) {
      auto f = m["compressed_fields"][group];
      Triplet a;
      if (!pos &&
          m["compressed_fields"]
                  .value("velocity_order", Json::object())
                  .value("order_encoding", "") == "lcp_blockwise_packed") {
        auto packed = load_field(m["compressed_fields"].at("velocity_order"));
        auto huffman =
            load_field(m["compressed_fields"].at("velocity_block_ids"));
        auto ids = huffman_decode(huffman.bytes, count);
        a = decode_blockwise(read_bytes(f.at("path").get<std::string>()),
                             packed, ids);
      } else
        a = decode_triplet(read_bytes(f.at("path").get<std::string>()), f,
                           o.vel_chunk_workers);
      for (size_t k = 0; k < 3; ++k)
        decoded[(pos ? 1 : 4) + k] = std::move(a[k]);
    }
  }
  auto restore = [&](size_t first, const std::string &key) {
    if (!m["compressed_fields"].contains(key))
      return;
    auto f = m["compressed_fields"][key];
    if (f.value("applied_during_lcp_decompression", false))
      return;
    auto indices = load_field(f);
    size_t chunk = f.value("chunk_size", size_t(0));
    std::vector<size_t> inverse(count);
    std::vector<bool> seen(count, false);
    for (size_t i = 0; i < count; ++i) {
      auto value = indices.number(i);
      size_t base = chunk ? (i / chunk) * chunk : 0;
      size_t local_count = chunk ? std::min(chunk, count - base) : count;
      if (value < 0 || value >= local_count)
        throw std::runtime_error("Invalid order index");
      size_t target = base + static_cast<size_t>(value);
      if (seen[target])
        throw std::runtime_error("Duplicate order index");
      seen[target] = true;
      inverse[target] = i;
    }
    for (size_t k = 0; k < 3; ++k)
      decoded[first + k] = decoded[first + k].gather(inverse);
  };
  restore(1, "order");
  restore(4, "velocity_order");
  if (m.value("structured_layout", Json::object()).value("enabled", false)) {
    auto layout = m["structured_layout"];
    auto f = m["compressed_fields"].value("id", Json::object());
    if (f.value("codec", "") != hilbert_codec)
      throw std::runtime_error("Missing structured IDs");
    for (auto key : {"side", "id_base", "position_digit_axes", "lattice_bits",
                     "velocity_cell_bits"})
      if (f.at("structured_layout").at(key) != layout.at(key))
        throw std::runtime_error("Structured ID layout disagrees with package");
    auto order =
        hybrid_order(decoded[0], {decoded[1], decoded[2], decoded[3]}, layout);
    std::vector<size_t> inverse(count);
    for (size_t i = 0; i < count; ++i)
      inverse[order[i]] = i;
    for (size_t i = 4; i < 7; ++i)
      decoded[i] = decoded[i].gather(inverse);
  }

  for (size_t i = 0; i < 7; ++i) {
    std::string name = fields[i];
    auto a = std::move(decoded[i]);
    if (a.size() != count)
      throw std::runtime_error("Missing or invalid decoded field: " + name);
    auto path = work / "decompressed" / (name + "." + a.type.name + ".raw");
    write_bytes(path, a.bytes, o.force);
    paths[name] = path.string();
    auto target = type_of(m["fields"][name]["dtype"]);
    bool lossless = m["compressed_fields"].contains(name) &&
                    m["compressed_fields"][name]["codec"] == "pcodec";
    if (i > 0 && i < 4 && !lossless) {
      double scale = m["position_scale"]["value"];
      Array b{target, std::vector<uint8_t>(checked_bytes(count, target.bytes))};
      for (size_t j = 0; j < count; ++j) {
        long double v = static_cast<double>(a.number(j)) * scale;
        if (!target.floating) {
          v = std::nearbyint(v);
          long double low =
              target.sign ? -std::ldexp(1.L, target.bytes * 8 - 1) : 0;
          long double high =
              std::ldexp(1.L, target.bytes * 8 - (target.sign ? 1 : 0)) - 1;
          v = std::clamp(v, low, high);
        }
        b.set(j, v);
      }
      a = std::move(b);
    } else if (i >= 4 && m["compressed_fields"].contains("velocities"))
      a = a.convert(target);
    else if (target.name != a.type.name)
      throw std::runtime_error("Lossless dtype mismatch");
    s.data[i] = std::move(a);
  }
  auto recombine_started = Clock::now();
  write_hdf5(work / "reconstructed.h5", s, o.force);
  m["timing"]["recombine_h5_wall_seconds"] = elapsed(recombine_started);
  m["sizes"]["reconstructed_h5_file_bytes"] =
      fs::file_size(work / "reconstructed.h5");
  m["artifacts"]["decompressed"] = paths;
  m["artifacts"]["reconstructed_h5"] = (work / "reconstructed.h5").string();
  m["timing"]["decompress_and_recombine_wall_seconds"] = elapsed(started);
  update_sizes(m, work);
  write_json(work / "manifest.json", m);
}
int run(const Options &o) {
  auto start = Clock::now();
  if (!o.input.empty() && fs::is_directory(o.input))
    return run_directory(o);
  auto work = fs::path(o.work_dir);
  Json m;
  if (o.command == "decompress")
    m = read_json(work / "manifest.json");
  else {
    if (o.command == "roundtrip" && fs::exists(work / "metrics.json")) {
      require_output(work / "metrics.json", o.force);
      fs::remove(work / "metrics.json");
    }
    m = preprocess(o);
  }
  if (o.command == "compress" || o.command == "roundtrip")
    compress(o, m);
  if (o.command == "decompress" || o.command == "roundtrip") {
    decompress(o, m);
    if (o.metrics) {
      auto report = compute_metrics(o, m);
      report["timing"]["roundtrip_wall_seconds"] = elapsed(start);
      write_json(work / "metrics.json", report);
      std::cout << "metrics = " << (work / "metrics.json").string() << '\n';
    }
    if (o.clean_raw) {
      fs::remove_all(work / "preprocessed");
      fs::remove_all(work / "decompressed");
    }
    if (o.command == "roundtrip")
      m["timing"]["roundtrip_wall_seconds"] = elapsed(start);
    update_sizes(m, work);
    write_json(work / "manifest.json", m);
    std::cout << "reconstructed_h5 = "
              << m["artifacts"]["reconstructed_h5"].get<std::string>() << '\n';
  }
  std::cout << "package_dir = " << work.string()
            << "\nmanifest = " << (work / "manifest.json").string() << '\n';
  if (m["sizes"].contains("payload_compression_ratio"))
    std::cout << "payload_CR = " << m["sizes"]["payload_compression_ratio"]
              << '\n';
  return 0;
}
} // namespace particle
