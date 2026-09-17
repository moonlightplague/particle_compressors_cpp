#pragma once
#include "particle/triplets.hpp"
namespace particle {
inline constexpr const char *lattice_name = "periodic_dense_id_lattice_v1";
inline constexpr const char *residual_name = "periodic_lattice_residual_v1";
inline constexpr const char *hybrid_name = "eulerian_lagrangian_hybrid_v1";
inline constexpr const char *hilbert_codec = "lattice_hilbert_pcodec_v1";
struct Lattice {
  size_t side = 0;
  int base = 0;
  std::array<size_t, 3> starts{}, shape{}, axes{0, 1, 2};
  std::vector<std::array<size_t, 3>> coordinates;
  std::vector<size_t> indices;
  size_t dense_count() const;
  Json metadata() const;
};
Lattice infer_lattice(const Array &ids, const Triplet &pos, size_t side,
                      double occupancy, bool structured = false);
Lattice load_lattice(const Array &ids, const Json &meta);
Array lattice_encode(const Array &values, const Lattice &layout, size_t axis,
                     bool residual, Array &wraps, double &error);
Json compress_lattice_field(const Options &o, Json &manifest,
                            const std::string &name, const Array &data,
                            const Lattice &lattice, size_t axis, bool position);
Array lattice_decode(const Array &dense, const Lattice &layout, size_t axis,
                     const std::string &transform, const Array *wraps);
Json structured_metadata(const Lattice &layout, int cell_bits);
Array hilbert_ids(const Array &ids, const Json &layout, bool inverse,
                  const Type &target);
std::vector<size_t> hybrid_order(const Array &ids, const Triplet &decoded,
                                 const Json &layout);
std::vector<uint8_t> encode_shaped(const Array &array, const std::string &codec,
                                   double bound,
                                   const std::array<size_t, 3> &shape,
                                   bool search, Json &metadata);
Array decode_shaped(const std::vector<uint8_t> &bytes, const Json &field);
std::vector<uint8_t> encode_dimensions(const Array &a, const std::string &codec,
                                       double bound,
                                       const std::array<size_t, 3> &shape);
}  // namespace particle
