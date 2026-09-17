#pragma once
#include "particle/pipeline.hpp"
namespace particle {
std::vector<uint8_t> huffman_encode(const Array &values, Json &metadata);
Array huffman_decode(const std::vector<uint8_t> &bytes, size_t count);
}  // namespace particle
