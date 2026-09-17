#include "particle/pipeline.hpp"
#include <bit>
#include <fstream>
#include <stdexcept>
namespace particle {
Snapshot read_native(const fs::path &path, int64_t limit) {
  if (std::endian::native != std::endian::little)
    throw std::runtime_error("Native snapshots require a little-endian host");
  auto filename = path.filename().string();
  if (!filename.starts_with("dat_"))
    throw std::runtime_error("Expected HDF5 or native dat_* input");
  auto config = path.parent_path() / ("cfg_" + filename.substr(4));
  std::ifstream f(config);
  if (!f)
    throw std::runtime_error("Missing native configuration: " +
                             config.string());
  std::vector<std::string> tokens;
  std::string token;
  while (f >> token)
    tokens.push_back(token);
  if (tokens.size() != 43)
    throw std::runtime_error("Native configuration must contain 43 values");
  std::array<int64_t, 43> values{};
  for (size_t i = 0; i < 43; ++i) {
    size_t used = 0;
    if (i >= 3 && i < 9)
      std::stod(tokens[i], &used);
    else
      values[i] = std::stoll(tokens[i], &used);
    if (used != tokens[i].size())
      throw std::runtime_error("Invalid native configuration value");
  }
  auto rank = values[0], procs = values[1], total = values[2], n = values[9],
       mesh = values[41], bitwidth = values[42];
  if (procs <= 0 || rank < 0 || rank >= procs || n <= 0 || total < n ||
      mesh <= 0 || bitwidth <= 0)
    throw std::runtime_error("Invalid native snapshot header");
  if (fs::file_size(path) != checked_bytes(n, 32))
    throw std::runtime_error("Native snapshot byte count mismatch");
  size_t count = limit > 0 ? std::min<size_t>(n, limit) : size_t(n);
  Snapshot s;
  s.metadata = {{"fields", Json::object()}, {"root_attrs", Json::object()}};
  auto addattr = [&](const std::string &name, int64_t v,
                     const std::string &dtype = "int32") {
    s.metadata["root_attrs"][name] = {
        {"dtype", dtype}, {"shape", Json::array()}, {"value", v}};
  };
  addattr("bitwidth", bitwidth);
  addattr("npart", count);
  addattr("npart_total", total, "uint64");
  addattr("nsidemesh", mesh);
  addattr("proc_size", procs);
  addattr("rank", rank);
  Json source = {{"format", "cfg_dat_field_major_v1"},
                 {"config_file", config.string()},
                 {"data_file", path.string()},
                 {"data_file_bytes", fs::file_size(path)},
                 {"storage_order", "field_major"},
                 {"byte_order", "little_endian"},
                 {"bytes_per_particle", 32},
                 {"fields", Json::object()}};
  std::ifstream input(path, std::ios::binary);
  size_t offset = 0;
  for (size_t index : {1, 2, 3, 4, 5, 6, 0}) {
    std::string name = index == 0 ? "id"
                       : index < 4
                           ? "pos" + std::string(fields[index])
                           : "vel" + std::string(fields[index]).substr(1);
    auto t = type_of(index == 0 ? "uint64" : index < 4 ? "int32" : "float32");
    s.data[index] = {t, std::vector<uint8_t>(checked_bytes(count, t.bytes))};
    input.seekg(offset);
    if (!input.read(reinterpret_cast<char *>(s.data[index].bytes.data()),
                    s.data[index].bytes.size()))
      throw std::runtime_error("Native snapshot read failed");
    s.metadata["fields"][fields[index]] = {{"h5_path", name},
                                           {"dtype", t.name},
                                           {"shape", {n}},
                                           {"selected_shape", {count}},
                                           {"attrs", Json::object()}};
    source["fields"][name] = {{"byte_offset", offset},
                              {"byte_count", checked_bytes(n, t.bytes)},
                              {"dtype", t.name},
                              {"shape", {n}}};
    offset += checked_bytes(n, t.bytes);
  }
  s.metadata["source"] = source;
  return s;
}
} // namespace particle
