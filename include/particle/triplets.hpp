#pragma once
#include "particle/pipeline.hpp"
namespace particle {
using Triplet = std::array<Array, 3>;
struct TripletEncoded {
  std::vector<uint8_t> bytes;
  std::vector<size_t> order;
  Json metadata;
  Array packed_order, block_ids;
};
TripletEncoded encode_triplet(const Triplet &data, const std::string &codec,
                              double bound, size_t chunk_size = 0,
                              bool hilbert_curve = false, int workers = 0,
                              bool position = true);
TripletEncoded encode_blockwise(const Triplet &data, double bound);
Triplet decode_blockwise(const std::vector<uint8_t> &bytes, const Array &packed,
                         const Array &block_ids);
Triplet decode_triplet(const std::vector<uint8_t> &bytes, const Json &metadata,
                       int workers = 0);
}  // namespace particle
