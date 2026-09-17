#include "particle/layout.hpp"
#include <algorithm>
#include <cstring>
#include <stdexcept>
namespace particle {
namespace {
size_t product(const std::array<size_t, 3> &s) {
  if (!s[0] || !s[1] || !s[2])
    throw std::runtime_error("Invalid dense shape");
  return checked_bytes(checked_bytes(s[0], s[1]), s[2]);
}
std::vector<size_t> transpose_order(const std::array<size_t, 3> &shape,
                                    const std::array<size_t, 3> &permutation,
                                    const std::array<bool, 3> &flips) {
  auto p = permutation;
  std::sort(p.begin(), p.end());
  if (p != std::array<size_t, 3>{0, 1, 2})
    throw std::runtime_error("Invalid axis permutation");
  std::array<size_t, 3> dims{shape[permutation[0]], shape[permutation[1]],
                             shape[permutation[2]]};
  std::vector<size_t> order(product(dims));
  for (size_t i = 0; i < order.size(); ++i) {
    std::array<size_t, 3> c{i / (dims[1] * dims[2]), (i / dims[2]) % dims[1],
                            i % dims[2]},
        base;
    for (size_t k = 0; k < 3; ++k)
      base[permutation[k]] = flips[k] ? dims[k] - c[k] - 1 : c[k];
    order[i] = (base[0] * shape[1] + base[1]) * shape[2] + base[2];
  }
  return order;
}
} // namespace
std::vector<uint8_t> encode_shaped(const Array &a, const std::string &codec,
                                   double bound,
                                   const std::array<size_t, 3> &shape,
                                   bool search, Json &m) {
  if (a.size() != product(shape))
    throw std::runtime_error("Invalid shaped array size");
  std::vector<uint8_t> best;
  std::array<size_t, 3> selected{0, 1, 2}, permutation{0, 1, 2};
  std::array<bool, 3> selected_flips{};
  auto trial = [&](const std::array<size_t, 3> &p,
                   const std::array<bool, 3> &flips) {
    auto data = a.gather(transpose_order(shape, p, flips));
    std::array<size_t, 3> dims{shape[p[0]], shape[p[1]], shape[p[2]]};
    auto bytes = encode_dimensions(data, codec, bound, dims);
    if (best.empty() || bytes.size() < best.size()) {
      best = std::move(bytes);
      selected = p;
      selected_flips = flips;
    }
  };
  do {
    trial(permutation, {});
    if (!search)
      break;
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  if (search) {
    auto p = selected;
    for (int bits = 1; bits < 8; ++bits)
      trial(p, {bool(bits & 4), bool(bits & 2), bool(bits & 1)});
  }
  m.update({{"encoded_count", a.size()},
            {"base_shape", shape},
            {"encoded_shape",
             {shape[selected[0]], shape[selected[1]], shape[selected[2]]}},
            {"axis_permutation", selected},
            {"axis_flips", selected_flips},
            {"axis_search", search}});
  return best;
}
Array decode_shaped(const std::vector<uint8_t> &bytes, const Json &f) {
  auto shape = f.at("base_shape").get<std::array<size_t, 3>>(),
       encoded = f.at("encoded_shape").get<std::array<size_t, 3>>(),
       permutation = f.at("axis_permutation").get<std::array<size_t, 3>>();
  auto flips = f.value("axis_flips", std::array<bool, 3>{});
  size_t n = f.at("encoded_count");
  if (product(shape) != n || product(encoded) != n)
    throw std::runtime_error("Shaped metadata count mismatch");
  auto order = transpose_order(shape, permutation, flips);
  for (size_t k = 0; k < 3; ++k)
    if (encoded[k] != shape[permutation[k]])
      throw std::runtime_error("Encoded dimensions disagree with permutation");
  auto a = decode(bytes, type_of(f.at("dtype")), n, n, f.at("codec"));
  std::vector<size_t> inverse(n);
  for (size_t i = 0; i < n; ++i)
    inverse[order[i]] = i;
  return a.gather(inverse);
}
} // namespace particle
