#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>

#include "CLI11.hpp"
#include "particle/pipeline.hpp"
namespace particle {
namespace {
std::string absolute_path(std::string path,
                          const fs::path &base = fs::current_path()) {
  if (path.starts_with("~/")) {
    const char *home = std::getenv("HOME");
    if (home) path = std::string(home) + path.substr(1);
  }
  fs::path p(path);
  return fs::absolute(p.is_absolute() ? p : base / p)
      .lexically_normal()
      .string();
}
void load_config(Options &o) {
  auto path = fs::path(o.config);
  if (!fs::is_regular_file(path))
    throw std::runtime_error("Config file does not exist: " + o.config);
  auto root = YAML::LoadFile(o.config);
  if (root.IsNull()) return;
  if (!root.IsMap()) throw std::runtime_error("Config must contain a mapping");
  for (auto it : root)
    if (it.first.as<std::string>() != "advanced")
      throw std::runtime_error("Unknown config section: " +
                               it.first.as<std::string>());
  auto a = root["advanced"];
  if (!a) return;
  if (!a.IsMap())
    throw std::runtime_error("Config section advanced must be a mapping");
  const std::set<std::string> keys{"lcp",
                                   "xnyzip",
                                   "abs_eb",
                                   "rel_eb",
                                   "pos_abs_eb",
                                   "pos_rel_eb",
                                   "vel_abs_eb",
                                   "vel_rel_eb",
                                   "vel_chunk_size",
                                   "vel_chunk_workers",
                                   "id_abs_eb",
                                   "position_scale",
                                   "position_scale_attr",
                                   "position_scale_value",
                                   "pos_compressor",
                                   "vel_compressor",
                                   "lossless",
                                   "lattice_layout",
                                   "lattice_min_occupancy",
                                   "lattice_axis_search",
                                   "xnyzip_structure_aware",
                                   "xnyzip_velocity_cell_bits"};
  for (auto it : a)
    if (!keys.contains(it.first.as<std::string>()))
      throw std::runtime_error("Unknown advanced config key: " +
                               it.first.as<std::string>());
#define OPTIONAL(K)                   \
  if (a[#K])                          \
  o.K = a[#K].IsNull() ? std::nullopt \
                       : std::optional<double>(a[#K].as<double>())
  OPTIONAL(abs_eb);
  OPTIONAL(rel_eb);
  OPTIONAL(pos_abs_eb);
  OPTIONAL(pos_rel_eb);
  OPTIONAL(vel_abs_eb);
  OPTIONAL(vel_rel_eb);
  OPTIONAL(position_scale_value);
#undef OPTIONAL
#define VALUE(K, T) \
  if (a[#K]) o.K = a[#K].as<T>()
  VALUE(lcp, std::string);
  VALUE(xnyzip, std::string);
  VALUE(position_scale, std::string);
  VALUE(position_scale_attr, std::string);
  VALUE(pos_compressor, std::string);
  VALUE(vel_compressor, std::string);
  VALUE(lossless, std::string);
  VALUE(vel_chunk_size, int64_t);
  VALUE(vel_chunk_workers, int);
  VALUE(id_abs_eb, double);
  VALUE(lattice_min_occupancy, double);
  VALUE(xnyzip_velocity_cell_bits, int);
  VALUE(lattice_layout, bool);
  VALUE(lattice_axis_search, bool);
  VALUE(xnyzip_structure_aware, bool);
#undef VALUE
  o.lcp = absolute_path(o.lcp, path.parent_path());
  o.xnyzip = absolute_path(o.xnyzip, path.parent_path());
}
}  // namespace
Options parse_cli(int argc, char **argv) {
  Options o;
  o.config = PARTICLE_DEFAULT_CONFIG;
  o.lcp = "tools/LCP/build/bin/lcp";
  o.xnyzip = "tools/XnYZip/build/XnYZip";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--config" && i + 1 < argc)
      o.config = argv[++i];
    else if (a.starts_with("--config="))
      o.config = a.substr(9);
  }
  o.config = absolute_path(o.config);
  load_config(o);
  CLI::App app{"Particle compression pipeline (C++ native implementation)"};
  app.require_subcommand(1);
  app.add_option("--config", o.config);
  for (auto name : {"preprocess", "compress", "decompress", "roundtrip"}) {
    auto *sub = app.add_subcommand(name);
    sub->callback([&o, name] { o.command = name; });
    sub->add_option("--config", o.config);
    sub->add_option("--work-dir", o.work_dir);
    sub->add_option("--lcp", o.lcp);
    sub->add_option("--xnyzip", o.xnyzip);
    sub->add_option("--vel-chunk-workers", o.vel_chunk_workers)
        ->check(CLI::NonNegativeNumber);
    sub->add_option("--field-workers", o.field_workers)
        ->check(CLI::NonNegativeNumber);
    sub->add_flag("--force", o.force);
    sub->add_flag("--clean-raw", o.clean_raw);
    if (std::string(name) == "decompress") continue;
    sub->add_option("input_h5", o.input)->required();
    sub->add_option("--file-workers", o.file_workers)
        ->check(CLI::NonNegativeNumber);
    sub->add_flag("--merge", o.merge);
    if (std::string(name) == "roundtrip") sub->add_flag("--metrics", o.metrics);
#define OPT(K, N) sub->add_option(N, o.K)
    OPT(abs_eb, "--abs-eb");
    OPT(rel_eb, "--rel-eb");
    OPT(pos_abs_eb, "--pos-abs-eb");
    OPT(pos_rel_eb, "--pos-rel-eb");
    OPT(vel_abs_eb, "--vel-abs-eb");
    OPT(vel_rel_eb, "--vel-rel-eb");
    OPT(id_abs_eb, "--id-abs-eb");
    OPT(position_scale_value, "--position-scale-value");
    OPT(position_scale_attr, "--position-scale-attr");
    OPT(vel_chunk_size, "--vel-chunk-size")->check(CLI::NonNegativeNumber);
    OPT(limit, "--limit")->check(CLI::PositiveNumber);
    OPT(lattice_min_occupancy, "--lattice-min-occupancy");
    OPT(xnyzip_velocity_cell_bits, "--xnyzip-velocity-cell-bits");
    OPT(position_scale, "--position-scale")
        ->check(CLI::IsMember({"auto", "raw", "attr", "value"}));
    OPT(pos_compressor, "--pos-compressor")
        ->check(CLI::IsMember({"lcp", "xnyzip", "pcodec", "sz3", "szo"}));
    OPT(vel_compressor, "--vel-compressor")
        ->check(CLI::IsMember({"lcp", "xnyzip", "pcodec", "sz3", "szo"}));
    OPT(lossless, "--lossless")->check(CLI::IsMember({"pcodec"}));
#undef OPT
    sub->add_flag("--sort", o.sort);
    sub->add_flag("--blockwise-ord", o.blockwise_ord);
    sub->add_flag("--lattice-layout,!--no-lattice-layout", o.lattice_layout);
    sub->add_flag("--lattice-axis-search,!--no-lattice-axis-search",
                  o.lattice_axis_search);
    sub->add_flag("--xnyzip-structure-aware,!--no-xnyzip-structure-aware",
                  o.xnyzip_structure_aware);
  }
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    int code = app.exit(e);
    if (code == 0) throw code;
    throw std::runtime_error("Invalid command line");
  }
  o.work_dir = absolute_path(o.work_dir);
  if (!o.input.empty()) o.input = absolute_path(o.input);
  validate(o);
  return o;
}
void validate(const Options &o) {
  if (o.command == "decompress") return;
  if (o.lattice_layout && o.xnyzip_structure_aware)
    throw std::runtime_error(
        "--xnyzip-structure-aware cannot be combined with --lattice-layout");
  auto known = [](const std::string &c) {
    return c == "lcp" || c == "xnyzip" || c == "szo" || c == "sz3" ||
           c == "pcodec";
  };
  if (!known(o.pos_compressor) || !known(o.vel_compressor) ||
      o.lossless != "pcodec")
    throw std::runtime_error("Invalid configured compressor");
  if (o.vel_compressor == "lcp" && o.pos_compressor != "lcp")
    throw std::runtime_error(
        "--vel-compressor lcp requires --pos-compressor lcp.");
  if (o.vel_compressor == "xnyzip" && o.pos_compressor != "lcp" &&
      o.pos_compressor != "xnyzip")
    throw std::runtime_error(
        "--vel-compressor xnyzip requires --pos-compressor lcp or xnyzip.");
  if (o.vel_chunk_size < 0 || o.vel_chunk_workers < 0 || o.field_workers < 0 ||
      o.file_workers < 0)
    throw std::runtime_error(
        "Worker counts and chunk size must be non-negative");
  if (o.vel_chunk_size && o.vel_compressor != "lcp" &&
      o.vel_compressor != "xnyzip")
    throw std::runtime_error(
        "--vel-chunk-size requires a triplet velocity compressor");
  if (o.blockwise_ord && (o.pos_compressor != "lcp" ||
                          o.vel_compressor != "lcp" || o.vel_chunk_size))
    throw std::runtime_error(
        "--blockwise-ord requires unchunked LCP positions and velocities");
  if (!(o.lattice_min_occupancy > 0 && o.lattice_min_occupancy <= 1) ||
      o.xnyzip_velocity_cell_bits < 1 || o.xnyzip_velocity_cell_bits > 10)
    throw std::runtime_error("Invalid layout configuration");
  if (o.merge && !fs::is_directory(o.input))
    throw std::runtime_error("--merge requires input_h5 to be a directory.");
}
}  // namespace particle
