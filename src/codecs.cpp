#include "cpcodec.h"
#include "particle/pipeline.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
extern "C" {
unsigned char *particle_szo_lorenzo(int, void *, size_t, double, size_t *);
unsigned char *particle_sz3_compress(int, void *, size_t *, int, double, double,
                                     double, size_t, size_t, size_t, size_t,
                                     size_t);
void *particle_sz3_decompress(int, unsigned char *, size_t, size_t, size_t,
                              size_t, size_t, size_t);
unsigned char *particle_szo_compress(int, void *, size_t *, int, double, double,
                                     double, size_t, size_t, size_t, size_t,
                                     size_t);
void *particle_szo_decompress(int, unsigned char *, size_t, size_t, size_t,
                              size_t, size_t, size_t);
}
namespace particle {
std::vector<uint8_t> encode(const Array &a, const std::string &codec,
                            double bound, size_t &encoded_count,
                            unsigned level) {
  encoded_count = a.size();

  if (codec == "pcodec") {
    size_t capacity = pco_standalone_guarantee_file_size(a.size(), a.type.pco);
    if (!capacity)
      throw std::runtime_error("pcodec rejected dtype or size");
    std::vector<uint8_t> output(capacity);
    size_t written = 0;
    PcoChunkConfig config{level, 0};
    if (pco_standalone_simple_compress_into(
            a.bytes.empty() ? static_cast<const void *>(&config)
                            : a.bytes.data(),
            a.size(), a.type.pco, &config, output.data(), output.size(),
            &written) != PcoSuccess)
      throw std::runtime_error("pcodec compression failed");
    output.resize(written);
    return output;
  }
  if (codec != "sz3" && codec != "pysz" && codec != "szo" &&
      codec != "szo_lorenzo")
    throw std::runtime_error("Unsupported field codec: " + codec);
  if (!a.size())
    throw std::runtime_error("Empty lossy field");
  if (!a.type.floating)
    throw std::runtime_error("Lossy codecs require floating-point fields");
  Array padded = a;
  encoded_count = std::max<size_t>(10000, a.size());
  padded.bytes.resize(checked_bytes(encoded_count, a.type.bytes));
  for (size_t i = a.size(); i < encoded_count; ++i)
    std::memcpy(padded.bytes.data() + i * a.type.bytes,
                a.bytes.data() + (a.size() - 1) * a.type.bytes, a.type.bytes);
  size_t n = 0;
  if (codec == "szo_lorenzo") {
    std::unique_ptr<unsigned char, decltype(&std::free)> p(
        particle_szo_lorenzo(a.type.bytes == 4 ? 0 : 1, padded.bytes.data(),
                             encoded_count, bound, &n),
        std::free);
    if (!p)
      throw std::runtime_error("SZO compression failed");
    return {p.get(), p.get() + n};
  }
  auto fn = codec == "szo" ? particle_szo_compress : particle_sz3_compress;
  std::unique_ptr<unsigned char, decltype(&std::free)> p(
      fn(a.type.bytes == 4 ? 0 : 1, padded.bytes.data(), &n, 0, bound, 0, 0, 0,
         0, 0, 0, encoded_count),
      std::free);
  if (!p)
    throw std::runtime_error(codec + " compression failed");
  return {p.get(), p.get() + n};
}
Array decode(const std::vector<uint8_t> &payload, const Type &type,
             size_t count, size_t encoded_count, const std::string &codec) {
  if (encoded_count < count)
    throw std::runtime_error("Encoded count is smaller than selected count");
  if (payload.empty())
    throw std::runtime_error("Empty compressed stream");
  Array out{type, std::vector<uint8_t>(checked_bytes(count, type.bytes))};
  if (codec == "pcodec") {
    size_t written = 0;
    if (pco_standalone_simple_decompress_into(
            payload.data(), payload.size(), type.pco,
            out.bytes.empty() ? static_cast<void *>(&written)
                              : out.bytes.data(),
            count, &written) != PcoSuccess ||
        written != count)
      throw std::runtime_error(
          "pcodec decompression dtype/count mismatch or invalid stream");
    return out;
  }
  if (codec != "pysz" && codec != "sz3" && codec != "szo")
    throw std::runtime_error("Unsupported field codec: " + codec);
  if (!type.floating)
    throw std::runtime_error("Lossy codecs require floating-point fields");
  auto fn = codec == "szo" ? particle_szo_decompress : particle_sz3_decompress;
  std::unique_ptr<void, decltype(&std::free)> p(
      fn(type.bytes == 4 ? 0 : 1, const_cast<unsigned char *>(payload.data()),
         payload.size(), 0, 0, 0, 0, encoded_count),
      std::free);
  if (!p)
    throw std::runtime_error(codec + " decompression failed");
  std::memcpy(out.bytes.data(), p.get(), out.bytes.size());
  return out;
}
} // namespace particle
namespace particle {
std::vector<uint8_t> encode_dimensions(const Array &a, const std::string &codec,
                                       double bound,
                                       const std::array<size_t, 3> &shape) {
  if (!a.type.floating || (codec != "sz3" && codec != "szo"))
    throw std::runtime_error("Invalid shaped codec or dtype");
  if (checked_bytes(checked_bytes(shape[0], shape[1]), shape[2]) != a.size())
    throw std::runtime_error("Shaped array count mismatch");
  auto fn = codec == "szo" ? particle_szo_compress : particle_sz3_compress;
  size_t n = 0;
  std::unique_ptr<unsigned char, decltype(&std::free)> p(
      fn(a.type.bytes == 4 ? 0 : 1, const_cast<uint8_t *>(a.bytes.data()), &n,
         0, bound, 0, 0, 0, 0, shape[0], shape[1], shape[2]),
      std::free);
  if (!p)
    throw std::runtime_error("Shaped compression failed");
  return {p.get(), p.get() + n};
}
} // namespace particle
