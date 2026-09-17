#include "particle/layout.hpp"
#include <algorithm>
#include <cmath>
namespace particle {
Json compress_lattice_field(const Options &o, Json &m, const std::string &name,
                            const Array &data, const Lattice &lattice,
                            size_t axis, bool position) {
  std::string codec = position ? o.pos_compressor : o.vel_compressor;
  double requested =
             m["field_error_bounds"][name][position ? "compressor_abs" : "abs"],
         guard = 0, error = 0;
  Array wraps;
  bool residual = position;
  auto dense = lattice_encode(data, lattice, axis, residual, wraps, error);
  if (position) {
    double magnitude = 1;
    for (size_t i = 0; i < data.size(); ++i)
      magnitude = std::max(magnitude, std::abs(double(data.number(i))));
    guard =
        std::max(error, 2. * std::numeric_limits<float>::epsilon() * magnitude);
    if (requested - guard <= 0) {
      residual = false;
      guard = 0;
      dense = lattice_encode(data, lattice, axis, false, wraps, error);
    }
  }
  double bound = requested - guard;
  auto path = m["artifacts"]["compressed"][name].get<std::string>();
  Json f = {{"field", name},
            {"codec", codec == "sz3" ? "pysz" : codec},
            {"dtype", data.type.name},
            {"count", data.size()},
            {"path", path},
            {"abs_error_bound", bound},
            {"spatial_layout", lattice_name},
            {"lattice_transform", residual ? residual_name : "identity"}};
  auto payload = encode_shaped(dense, codec, bound, lattice.shape,
                               o.lattice_axis_search, f);
  write_bytes(path, payload, o.force);
  f["bytes"] = payload.size();
  if (position) {
    f["requested_compressor_abs"] = requested;
    f["transform_roundoff_guard"] = guard;
    f["transform_roundtrip_max_abs"] = error;
  }
  if (residual) {
    size_t encoded;
    auto bytes = encode(wraps, "pcodec", 0, encoded);
    auto wrap_path =
        fs::path(o.work_dir) / "compressed" / (name + ".lattice-wrap.pco");
    write_bytes(wrap_path, bytes, o.force);
    Json sidecar = {{"field", name + "_lattice_wrap"},
                    {"codec", "pcodec"},
                    {"dtype", "int8"},
                    {"count", data.size()},
                    {"path", wrap_path.string()},
                    {"bytes", bytes.size()}};
    f["lattice_wrap_field"] = sidecar;
    m["artifacts"]["compressed"][name + "_lattice_wrap"] = wrap_path.string();
  }
  return f;
}
} // namespace particle
