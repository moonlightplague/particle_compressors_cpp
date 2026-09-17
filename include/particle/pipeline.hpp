#pragma once
#include <H5Cpp.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "json.hpp"
namespace particle {
using Json = nlohmann::ordered_json;
namespace fs = std::filesystem;
inline constexpr std::array<const char *, 7> fields{"id", "x",  "y", "z",
                                                    "vx", "vy", "vz"};
struct Type {
  std::string name;
  size_t bytes = 0;
  uint8_t pco = 0;
  bool floating = false, sign = false;
  H5::DataType h5() const;
};
Type type_of(const std::string &name);
Type type_of(const H5::DataType &type);
size_t checked_bytes(size_t count, size_t width);
struct Array {
  Type type;
  std::vector<uint8_t> bytes;
  size_t size() const { return type.bytes ? bytes.size() / type.bytes : 0; }
  long double number(size_t i) const;
  void set(size_t i, long double value);
  Array convert(const Type &target) const;
  Array gather(const std::vector<size_t> &order) const;
};
struct Snapshot {
  std::array<Array, 7> data;
  Json metadata;
  size_t count() const { return data[0].size(); }
};
struct Options {
  Json merge_metadata = Json::object();
  std::string command, input, work_dir = "particle_pipeline_runs", config;
  std::string pos_compressor = "lcp", vel_compressor = "sz3",
              lossless = "pcodec", lcp, xnyzip;
  std::string position_scale = "auto", position_scale_attr = "bitwidth";
  std::optional<double> abs_eb, rel_eb = 1e-3, pos_abs_eb, pos_rel_eb,
                                vel_abs_eb, vel_rel_eb, position_scale_value;
  double id_abs_eb = 0, lattice_min_occupancy = .8;
  int64_t limit = 0, vel_chunk_size = 0;
  int field_workers = 0, file_workers = 0, vel_chunk_workers = 0,
      xnyzip_velocity_cell_bits = 7;
  bool force = false, clean_raw = false, sort = false, metrics = false,
       merge = false;
  bool lattice_layout = false, lattice_axis_search = true,
       xnyzip_structure_aware = false, blockwise_ord = false;
};
Options parse_cli(int argc, char **argv);
void validate(const Options &o);
Snapshot read_native(const fs::path &path, int64_t limit = 0);
Json compute_metrics(const Options &o, const Json &manifest);
int run_directory(const Options &o);
Snapshot read_hdf5(const fs::path &path, int64_t limit = 0);
void write_hdf5(const fs::path &path, const Snapshot &snapshot, bool force);
Json attributes(const H5::H5Object &object);
void apply_attributes(H5::H5Object &object, const Json &attrs);
std::vector<uint8_t> read_bytes(const fs::path &path);
void write_bytes(const fs::path &path, const std::vector<uint8_t> &data,
                 bool force);
void require_output(const fs::path &path, bool force);
std::string json_text(const Json &value);
Json read_json(const fs::path &path);
void write_json(const fs::path &path, const Json &data);
std::vector<uint8_t> encode(const Array &data, const std::string &codec,
                            double bound, size_t &encoded_count,
                            unsigned level = 8);
Array decode(const std::vector<uint8_t> &payload, const Type &type,
             size_t count, size_t encoded_count, const std::string &codec);
std::vector<size_t> stable_id_order(const Array &ids);
Json preprocess(const Options &o);
void compress(const Options &o, Json &manifest);
void decompress(const Options &o, Json &manifest);
void update_sizes(Json &manifest, const fs::path &work);
int run(const Options &o);
}  // namespace particle
