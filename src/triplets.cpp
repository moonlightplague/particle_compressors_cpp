#include "particle/triplets.hpp"

#include <SZ3/api/lcp.hpp>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <xnyzip/xnyzip.hpp>

#include "particle/runtime.hpp"
namespace particle {
namespace {
std::array<std::vector<float>, 3> floats(const Triplet &t) {
  std::array<std::vector<float>, 3> v;
  for (size_t k = 0; k < 3; ++k) {
    v[k].resize(t[k].size());
    for (size_t i = 0; i < v[k].size(); ++i)
      v[k][i] = static_cast<float>(t[k].number(i));
  }
  return v;
}
Triplet arrays(const std::array<std::vector<float>, 3> &v) {
  Triplet t;
  for (size_t k = 0; k < 3; ++k) {
    t[k] = {type_of("float32"),
            std::vector<uint8_t>(checked_bytes(v[k].size(), 4))};
    std::memcpy(t[k].bytes.data(), v[k].data(), t[k].bytes.size());
  }
  return t;
}
void append64(std::vector<uint8_t> &b, uint64_t v) {
  for (int i = 0; i < 8; ++i) b.push_back(static_cast<uint8_t>(v >> (i * 8)));
}
uint64_t take64(const std::vector<uint8_t> &b, size_t &p) {
  if (p > b.size() || b.size() - p < 8)
    throw std::runtime_error("Truncated chunk container");
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v |= uint64_t(b[p++]) << (i * 8);
  return v;
}
XnYZip::Options xoptions(const Json &f) {
  XnYZip::Options o;
  o.l2_bound = f.at("l2_error_bound").get<float>();
  o.direct_threshold = 0;
  o.quantizer = f.value("quantizer", "to") == "cube"
                    ? XnYZip::Quantizer::cube
                    : XnYZip::Quantizer::truncated_octahedron;
  return o;
}
Triplet decode_native(const std::vector<uint8_t> &b, const Json &f,
                      bool batch = false) {
  if (f.at("codec") == "lcp") {
    auto d = lcp::decompress<float>(
        b.data(), b.size(),
        batch ? lcp::Mode::independent_frames : lcp::Mode::frame);
    return arrays({d.x, d.y, d.z});
  }
  auto v = XnYZip::decompress(b, xoptions(f));
  if (v.size() % 3) throw std::runtime_error("Invalid decoded triplet");
  std::array<std::vector<float>, 3> a;
  for (size_t k = 0; k < 3; ++k) {
    a[k].resize(v.size() / 3);
    for (size_t i = 0; i < a[k].size(); ++i) a[k][i] = v[3 * i + k];
  }
  return arrays(a);
}
TripletEncoded encode_native(const Triplet &t, const std::string &codec,
                             double bound, bool hilbert_curve = false,
                             bool position = true) {
  TripletEncoded out;
  size_t n = t[0].size();
  auto v = floats(t);
  out.metadata = {{"codec", codec}, {"dtype", "float32"}, {"count", n}};
  if (codec == "lcp") {
    lcp::Options o;
    o.absolute_error = bound;
    auto c = lcp::compress(v[0], v[1], v[2], o);
    out.bytes = std::move(c.bytes);
    out.order = std::move(c.order);
    out.metadata["abs_error_bound"] = bound;
    return out;
  }
  if (codec != "xnyzip") throw std::runtime_error("Unsupported triplet codec");
  std::vector<float> xyz(checked_bytes(n, 3));
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < 3; ++k) xyz[3 * i + k] = v[k][i];
  XnYZip::Options options;
  options.l2_bound = static_cast<float>(bound);
  options.direct_threshold = 0;
  options.curve =
      hilbert_curve ? XnYZip::Curve::hilbert : XnYZip::Curve::z_order;
  double requested = bound;
  bool accepted = false;
  double maximum = 0;
  int attempt = 0;
  for (; attempt < 6; ++attempt) {
    auto c = XnYZip::compress(xyz, options);
    auto restored = XnYZip::decompress(c.bytes, options);
    if (restored.size() != xyz.size())
      throw std::runtime_error("XnYZip count mismatch");
    maximum = 0;
    for (size_t i = 0; i < n; ++i) {
      if (c.order[i] >= n)
        throw std::runtime_error("Invalid XnYZip permutation");
      double squared = 0;
      for (size_t k = 0; k < 3; ++k) {
        double e = double(restored[3 * i + k]) - xyz[3 * c.order[i] + k];
        squared += e * e;
      }
      maximum = std::isfinite(squared)
                    ? std::max(maximum, std::sqrt(squared))
                    : std::numeric_limits<double>::infinity();
    }
    if (std::isfinite(maximum) && maximum <= requested) {
      out.bytes = std::move(c.bytes);
      out.order.assign(c.order.begin(), c.order.end());
      accepted = true;
      break;
    }
    if (position &&
        options.quantizer == XnYZip::Quantizer::truncated_octahedron) {
      options.quantizer = XnYZip::Quantizer::cube;
      continue;
    }
    double margin = std::max(.01 * requested, 2 * (maximum - requested));
    double next = options.l2_bound - margin;
    if (!(next > 0) || !std::isfinite(next)) break;
    options.l2_bound = static_cast<float>(next);
  }
  if (!accepted)
    throw std::runtime_error("XnYZip failed L2 validation after codec retries");
  out.metadata.update(
      {{"input_layout", "triplet_interleaved"},
       {"native_order_dtype", "uint64"},
       {"error_bound_norm", "l2"},
       {"l2_error_bound", options.l2_bound},
       {"quantizer",
        options.quantizer == XnYZip::Quantizer::cube ? "cube" : "to"},
       {"curve", hilbert_curve ? "-h" : "-z"},
       {"storage_mode", "-rle"},
       {"direct_threshold", 0},
       {"validated_max_l2_error", maximum},
       {"validation_l2_bound", requested},
       {"compression_attempts", attempt + 1}});
  return out;
}
}  // namespace
TripletEncoded encode_triplet(const Triplet &t, const std::string &codec,
                              double bound, size_t chunk, bool hilbert_curve,
                              int workers, bool position) {
  if (!chunk) return encode_native(t, codec, bound, hilbert_curve, position);
  size_t n = t[0].size();
  TripletEncoded result;
  size_t chunks = n / chunk + (n % chunk != 0);
  std::string magic = codec == "lcp" ? std::string("LCPCHK2\0", 8)
                                     : std::string("XNYCHK1\0", 8);
  result.bytes.assign(magic.begin(), magic.end());
  append64(result.bytes, n);
  append64(result.bytes, chunk);
  append64(result.bytes, chunks);
  // A single set of XnYZip quantizer options must cover the whole container.
  // Retry decisions cannot differ per chunk because legacy metadata is global.
  auto chunks_encoded = parallel_map(
      chunks, worker_count(workers, chunks, 16), [&](size_t index) {
        size_t start = index * chunk, size = std::min(chunk, n - start);
        Triplet part;
        for (size_t k = 0; k < 3; ++k) {
          part[k].type = t[k].type;
          auto begin = t[k].bytes.begin() + start * t[k].type.bytes;
          part[k].bytes.assign(begin, begin + size * t[k].type.bytes);
        }
        return encode_native(part, codec, bound, false, false);
      });
  for (size_t index = 0; index < chunks; ++index) {
    size_t start = index * chunk, size = std::min(chunk, n - start);
    auto &c = chunks_encoded[index];
    if (start == 0)
      result.metadata = c.metadata;
    else if (codec == "xnyzip" &&
             (c.metadata["quantizer"] != result.metadata["quantizer"] ||
              c.metadata["l2_error_bound"] !=
                  result.metadata["l2_error_bound"]))
      throw std::runtime_error(
          "XnYZip chunks require inconsistent retry settings");
    append64(result.bytes, size);
    append64(result.bytes, c.bytes.size());
    result.bytes.insert(result.bytes.end(), c.bytes.begin(), c.bytes.end());
    result.order.insert(result.order.end(), c.order.begin(), c.order.end());
  }
  result.metadata.update(
      {{"count", n},
       {"chunk_size", chunk},
       {"chunk_count", chunks},
       {"container", codec == "lcp" ? "chunked_lcp_v2" : "chunked_xnyzip_v1"}});
  return result;
}
Triplet decode_triplet(const std::vector<uint8_t> &bytes, const Json &f,
                       int workers) {
  size_t count = f.at("count"), chunk = f.value("chunk_size", size_t(0));
  if (!chunk) {
    auto t = decode_native(bytes, f);
    if (t[0].size() != count)
      throw std::runtime_error("Decoded triplet count mismatch");
    return t;
  }
  std::string codec = f.at("codec"), magic = codec == "lcp"
                                                 ? std::string("LCPCHK2\0", 8)
                                                 : std::string("XNYCHK1\0", 8);
  if (bytes.size() < 32 ||
      !std::equal(magic.begin(), magic.end(), bytes.begin()))
    throw std::runtime_error("Invalid chunk magic");
  size_t p = 8;
  auto stored_count = take64(bytes, p), stored_chunk = take64(bytes, p),
       segments = take64(bytes, p);
  if (stored_count != count || stored_chunk != chunk || !segments ||
      segments > count)
    throw std::runtime_error("Invalid chunk header");
  Triplet out;
  for (auto &a : out) a.type = type_of("float32");
  struct Segment {
    size_t offset, length, count;
    bool batch;
  };
  std::vector<Segment> entries;
  size_t total = 0;
  for (size_t j = 0; j < segments; ++j) {
    size_t n = take64(bytes, p), len = take64(bytes, p);
    if (!n || n > count - total || len > bytes.size() - p || !len)
      throw std::runtime_error("Invalid chunk entry");
    bool batch = codec == "lcp" && n > chunk;
    if (batch && n % chunk)
      throw std::runtime_error("Invalid LCP batch segment");
    entries.push_back({p, len, n, batch});
    p += len;
    total += n;
  }
  if (total != count || p != bytes.size())
    throw std::runtime_error("Truncated or trailing chunk data");
  auto decoded = parallel_map(
      entries.size(), worker_count(workers, entries.size(), 16),
      [&](size_t index) {
        auto e = entries[index];
        std::vector<uint8_t> payload(bytes.begin() + e.offset,
                                     bytes.begin() + e.offset + e.length);
        auto t = decode_native(payload, f, e.batch);
        if (t[0].size() != e.count)
          throw std::runtime_error("Chunk count mismatch");
        return t;
      });
  for (auto &t : decoded)
    for (size_t k = 0; k < 3; ++k)
      out[k].bytes.insert(out[k].bytes.end(), t[k].bytes.begin(),
                          t[k].bytes.end());
  return out;
}
}  // namespace particle
namespace particle {
TripletEncoded encode_blockwise(const Triplet &t, double bound) {
  size_t n = t[0].size();
  auto v = floats(t);
  lcp::detail::Codec<float> codec{SZ3::HuffmanEncoder<size_t>(),
                                  SZ3::Lossless_zstd()};
  SZ3::Config config(n);
  config.absErrorBound = bound;
  std::vector<size_t> local(n), blocks(n), counts;
  size_t size = 0;
  std::unique_ptr<unsigned char[]> payload(codec.compress(
      config, v[0].data(), v[1].data(), v[2].data(), size, nullptr, 3, 0, 0, 0,
      local.data(), &counts, blocks.data()));
  TripletEncoded c;
  c.bytes.assign(payload.get(), payload.get() + size);
  c.metadata = {{"codec", "lcp"},
                {"dtype", "float32"},
                {"count", n},
                {"abs_error_bound", bound}};
  size_t bits = 0;
  for (auto count : counts)
    bits += checked_bytes(count, std::bit_width(count - 1));
  c.packed_order = {
      type_of("uint32"),
      std::vector<uint8_t>(checked_bytes(bits / 32 + (bits % 32 != 0), 4), 0)};
  size_t bit = 0, index = 0;
  for (auto count : counts) {
    size_t width = std::bit_width(count - 1);
    for (size_t i = 0; i < count; ++i, ++index) {
      if (local[index] >= count)
        throw std::runtime_error("Invalid block-local order");
      for (size_t k = 0; k < width; ++k, ++bit)
        c.packed_order.bytes[bit / 8] |= ((local[index] >> k) & 1) << (bit % 8);
    }
  }
  size_t maximum = *std::max_element(blocks.begin(), blocks.end());
  auto type = type_of(maximum <= UINT32_MAX ? "uint32" : "uint64");
  c.block_ids = {type, std::vector<uint8_t>(checked_bytes(n, type.bytes))};
  for (size_t i = 0; i < n; ++i) c.block_ids.set(i, blocks[i]);
  return c;
}
Triplet decode_blockwise(const std::vector<uint8_t> &bytes, const Array &packed,
                         const Array &ids) {
  lcp::detail::check_frame(bytes.data(), bytes.size());
  lcp::detail::Codec<float> codec{SZ3::HuffmanEncoder<size_t>(),
                                  SZ3::Lossless_zstd()};
  float *x = nullptr, *y = nullptr, *z = nullptr;
  size_t n = 0;
  std::vector<size_t> counts, blocks;
  const auto *p = bytes.data();
  codec.decompressWithoutAllocateMemory(p, x, y, z, n, bytes.size(), &counts,
                                        &blocks);
  std::unique_ptr<float[]> owner(x);
  if (n != ids.size() || counts.size() != blocks.size())
    throw std::runtime_error("Blockwise count mismatch");
  std::map<size_t, std::vector<size_t>> original;
  for (size_t i = 0; i < n; ++i)
    original[static_cast<size_t>(ids.number(i))].push_back(i);
  if (original.size() != blocks.size())
    throw std::runtime_error("Block-ID set mismatch");
  size_t bits = 0;
  for (auto count : counts) {
    if (!count) throw std::runtime_error("Empty LCP block");
    bits += checked_bytes(count, std::bit_width(count - 1));
  }
  if (packed.type.name != "uint32" ||
      packed.bytes.size() != checked_bytes(bits / 32 + (bits % 32 != 0), 4))
    throw std::runtime_error("Packed order size mismatch");
  std::array<std::vector<float>, 3> out;
  for (auto &a : out) a.resize(n);
  size_t bit = 0, index = 0;
  for (size_t block = 0; block < counts.size(); ++block) {
    size_t count = counts[block], width = std::bit_width(count - 1);
    auto it = original.find(blocks[block]);
    if (it == original.end() || it->second.size() != count)
      throw std::runtime_error("Block population mismatch");
    std::vector<bool> seen(count, false);
    for (size_t i = 0; i < count; ++i, ++index) {
      size_t local = 0;
      for (size_t k = 0; k < width; ++k, ++bit)
        local |= size_t((packed.bytes[bit / 8] >> (bit % 8)) & 1) << k;
      if (local >= count || seen[local])
        throw std::runtime_error("Invalid block-local permutation");
      seen[local] = true;
      size_t row = it->second[local];
      out[0][row] = x[index];
      out[1][row] = y[index];
      out[2][row] = z[index];
    }
  }
  if (index != n)
    throw std::runtime_error("LCP block counts do not sum to particle count");
  return arrays(out);
}
}  // namespace particle
