#include "particle/layout.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
namespace particle {
namespace {
double modulo(double v, double side) { return v - std::floor(v / side) * side; }
std::array<size_t, 3> coordinates(long double id, size_t side, int base) {
  size_t total = checked_bytes(checked_bytes(side, side), side);
  if (id < base || id >= static_cast<long double>(total) + base)
    throw std::runtime_error("IDs do not fit zero/one-based lattice");
  uint64_t v = static_cast<uint64_t>(id) - base;
  return {v / (side * side), (v / side) % side, v % side};
}
size_t index_of(const std::array<size_t, 3> &c, const Lattice &l) {
  std::array<size_t, 3> local;
  for (size_t k = 0; k < 3; ++k) {
    local[k] = (c[k] + l.side - l.starts[k]) % l.side;
    if (local[k] >= l.shape[k])
      throw std::runtime_error("Lattice coordinate exceeds dense shape");
  }
  return (local[0] * l.shape[1] + local[1]) * l.shape[2] + local[2];
}
uint32_t morton(const std::array<size_t, 3> &c) {
  uint32_t m = 0;
  for (int k = 0; k < 3; ++k) {
    if (c[k] >= 1024)
      throw std::runtime_error("Morton coordinate exceeds 10 bits");
    for (int b = 0; b < 10; ++b)
      m |= ((c[k] >> b) & 1) << (3 * b + k);
  }
  return m;
}
std::array<size_t, 3> unmorton(uint32_t m) {
  std::array<size_t, 3> c{};
  for (int k = 0; k < 3; ++k)
    for (int b = 0; b < 10; ++b)
      c[k] |= ((m >> (3 * b + k)) & 1) << b;
  return c;
}
uint32_t hilbert(uint32_t v, int bits, bool inverse) {
  if (bits < 1 || bits > 10)
    throw std::runtime_error("Invalid Hilbert bit count");
  if (inverse) {
    v ^= (v & 0x92492492U) >> 1;
    v ^= (v >> 1) & 0x92492492U;
  }
  if (bits > 1) {
    uint32_t block = bits * 3 - 3, h = (v >> block) & 7, shift = 0, sign = 0;
    while (block) {
      block -= 3;
      h <<= 2;
      uint32_t m = (0x20212021U >> h) & 3;
      shift = (0x48U >> (inverse ? 4 - shift + m : 7 - shift - m)) & 3;
      sign = (sign | (sign << 3)) >> m;
      sign = (sign ^ (0x53560300U >> h)) & 7;
      if (inverse) {
        h = (v >> block) & 7;
        m = h ^ sign;
        m = ((m | (m << 3)) >> shift) & 7;
      } else {
        m = (v >> block) & 7;
        h = ((m | (m << 3)) >> shift) & 7;
        h ^= sign;
      }
      v ^= (m ^ h) << block;
    }
  }
  if (!inverse) {
    v ^= (v >> 1) & 0x92492492U;
    v ^= (v & 0x92492492U) >> 1;
  }
  return v;
}
void validate_structured(const Json &j) {
  if (!j.value("enabled", false) || j.value("name", "") != hybrid_name)
    throw std::runtime_error("Invalid structured layout");
  size_t side = j.at("side");
  int bits = j.at("lattice_bits"), cell = j.at("velocity_cell_bits"),
      base = j.at("id_base");
  auto axes = j.at("position_digit_axes").get<std::array<size_t, 3>>();
  std::sort(axes.begin(), axes.end());
  if (side <= 1 || side > 1024 || bits != int(std::bit_width(side - 1)) ||
      cell < 1 || cell > 10 || base < 0 || base > 1 ||
      axes != std::array<size_t, 3>{0, 1, 2})
    throw std::runtime_error("Inconsistent structured geometry");
}
} // namespace
size_t Lattice::dense_count() const {
  return checked_bytes(checked_bytes(shape[0], shape[1]), shape[2]);
}
Json Lattice::metadata() const {
  return {{"name", lattice_name},
          {"enabled", true},
          {"side", side},
          {"id_base", base},
          {"periodic_starts", starts},
          {"dense_shape", shape},
          {"dense_count", dense_count()},
          {"particle_count", indices.size()},
          {"occupancy", double(indices.size()) / dense_count()},
          {"position_digit_axes", axes},
          {"hole_fill", "linear_flat_index"},
          {"implicit_full_lattice", false}};
}
Lattice infer_lattice(const Array &ids, const Triplet &pos, size_t side,
                      double occupancy, bool structured) {
  if (!side || ids.size() == 0)
    throw std::runtime_error("Lattice side/count must be positive");
  for (const auto &a : pos)
    for (size_t i = 0; i < a.size(); ++i)
      if (a.number(i) < 0 || a.number(i) >= 1 || !std::isfinite(a.number(i)))
        throw std::runtime_error(
            "Normalized positions are outside periodic [0, 1) domain");
  auto order = stable_id_order(ids);
  for (size_t i = 1; i < order.size(); ++i)
    if (ids.number(order[i - 1]) == ids.number(order[i]))
      throw std::runtime_error("Particles contain duplicate lattice IDs");
  size_t total = checked_bytes(checked_bytes(side, side), side);
  Lattice l;
  l.side = side;
  double best = std::numeric_limits<double>::infinity();
  bool found = false;
  size_t sample = std::min<size_t>(200000, ids.size());
  for (int base : {0, 1}) {
    if (ids.number(order.front()) < base ||
        ids.number(order.back()) >= static_cast<long double>(total) + base)
      continue;
    std::array<size_t, 3> axes{0, 1, 2};
    do {
      double score = 0;
      for (size_t k = 0; k < 3; ++k) {
        double squared = 0;
        for (size_t j = 0; j < sample; ++j) {
          size_t row =
              order[sample == 1 ? 0 : j * (ids.size() - 1) / (sample - 1)];
          auto c = coordinates(ids.number(row), side, base);
          double d =
              modulo(double(pos[k].number(row)) * side - c[axes[k]] + side / 2.,
                     side) -
              side / 2.;
          squared += d * d;
        }
        score += squared / sample;
      }
      if (score < best) {
        best = score;
        l.base = base;
        l.axes = axes;
        found = true;
      }
    } while (std::next_permutation(axes.begin(), axes.end()));
  }
  if (!found)
    throw std::runtime_error(
        "Particle IDs fit neither zero-based nor one-based lattice");
  for (size_t i = 0; i < ids.size(); ++i)
    l.coordinates.push_back(coordinates(ids.number(i), side, l.base));
  for (size_t k = 0; k < 3; ++k) {
    std::set<size_t> present;
    for (auto &c : l.coordinates)
      present.insert(c[k]);
    if (present.size() == side) {
      l.starts[k] = 0;
      l.shape[k] = side;
      continue;
    }
    std::vector<size_t> v(present.begin(), present.end());
    size_t gap = 0, start = v.front();
    for (size_t i = 0; i < v.size(); ++i) {
      size_t following = v[(i + 1) % v.size()],
             g = (following + side - v[i]) % side;
      if (g > gap) {
        gap = g;
        start = following;
      }
    }
    size_t span = 1;
    for (auto x : v)
      span = std::max(span, (x + side - start) % side + 1);
    l.starts[k] = start;
    l.shape[k] = span;
  }
  for (auto &c : l.coordinates)
    l.indices.push_back(index_of(c, l));
  if (!structured) {
    if (l.dense_count() < 10000)
      throw std::runtime_error("Dense box has fewer than 10000 values");
    if (double(ids.size()) / l.dense_count() < occupancy)
      throw std::runtime_error(
          "Dense lattice occupancy is below configured minimum");
  }
  return l;
}
Lattice load_lattice(const Array &ids, const Json &j) {
  if (j.value("name", "") != lattice_name || !j.value("enabled", false))
    throw std::runtime_error("Invalid lattice metadata");
  Lattice l;
  l.side = j.at("side");
  l.base = j.at("id_base");
  l.starts = j.at("periodic_starts").get<std::array<size_t, 3>>();
  l.shape = j.at("dense_shape").get<std::array<size_t, 3>>();
  l.axes = j.at("position_digit_axes").get<std::array<size_t, 3>>();
  auto axes = l.axes;
  std::sort(axes.begin(), axes.end());
  if (!l.side || l.base < 0 || l.base > 1 ||
      axes != std::array<size_t, 3>{0, 1, 2})
    throw std::runtime_error("Invalid lattice geometry");
  for (size_t k = 0; k < 3; ++k)
    if (l.starts[k] >= l.side || !l.shape[k] || l.shape[k] > l.side)
      throw std::runtime_error("Invalid dense shape");
  if (j.at("dense_count").get<size_t>() != l.dense_count() ||
      j.at("particle_count").get<size_t>() != ids.size())
    throw std::runtime_error("Lattice count mismatch");
  std::set<size_t> seen;
  for (size_t i = 0; i < ids.size(); ++i) {
    auto c = coordinates(ids.number(i), l.side, l.base);
    size_t index = index_of(c, l);
    if (!seen.insert(index).second)
      throw std::runtime_error("Duplicate lattice cell");
    l.coordinates.push_back(c);
    l.indices.push_back(index);
  }
  return l;
}
Array lattice_encode(const Array &a, const Lattice &l, size_t axis,
                     bool residual, Array &wraps, double &error) {
  Array transformed = a;
  error = 0;
  if (residual) {
    transformed = a.convert(type_of("float32"));
    wraps = {type_of("int8"), std::vector<uint8_t>(a.size())};
    for (size_t i = 0; i < a.size(); ++i) {
      double predictor = double(l.coordinates[i][l.axes[axis]]) / l.side;
      double r = double(a.number(i)) - predictor, w = modulo(r + .5, 1) - .5,
             offset = std::nearbyint(r - w);
      if (std::abs(offset) > 1)
        throw std::runtime_error("Unexpected lattice wrap offset");
      transformed.set(i, w);
      wraps.set(i, offset);
      double restored =
          float(double(transformed.number(i)) + predictor + offset);
      error = std::max(error, std::abs(restored - double(a.number(i))));
    }
  }
  Array dense{transformed.type, std::vector<uint8_t>(checked_bytes(
                                    l.dense_count(), transformed.type.bytes))};
  std::vector<size_t> order(a.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](size_t x, size_t y) { return l.indices[x] < l.indices[y]; });
  for (size_t i = 0; i < a.size(); ++i)
    std::memcpy(dense.bytes.data() + l.indices[i] * dense.type.bytes,
                transformed.bytes.data() + i * dense.type.bytes,
                dense.type.bytes);
  size_t next = 0;
  for (size_t index = 0; index < dense.size(); ++index) {
    while (next < order.size() && l.indices[order[next]] < index)
      ++next;
    if (next < order.size() && l.indices[order[next]] == index)
      continue;
    double value;
    if (!next)
      value = transformed.number(order[0]);
    else if (next == order.size())
      value = transformed.number(order.back());
    else {
      size_t left = order[next - 1], right = order[next];
      double x = transformed.number(left), y = transformed.number(right);
      value = x + (y - x) / double(l.indices[right] - l.indices[left]) *
                      double(index - l.indices[left]);
    }
    dense.set(index, value);
  }
  return dense;
}
Array lattice_decode(const Array &dense, const Lattice &l, size_t axis,
                     const std::string &transform, const Array *wraps) {
  if (dense.size() != l.dense_count())
    throw std::runtime_error("Dense decoded count mismatch");
  Array a = dense.gather(l.indices);
  if (transform == "identity")
    return a;
  if (transform != residual_name)
    throw std::runtime_error("Unknown lattice transform");
  if (wraps && wraps->size() != a.size())
    throw std::runtime_error("Wrap sidecar count mismatch");
  for (size_t i = 0; i < a.size(); ++i) {
    double value =
        double(a.number(i)) + double(l.coordinates[i][l.axes[axis]]) / l.side;
    value = wraps ? value + double(wraps->number(i)) : modulo(value, 1);
    a.set(i, value);
  }
  return a;
}
Json structured_metadata(const Lattice &l, int cells) {
  int bits = std::bit_width(l.side - 1);
  Json j = {{"name", hybrid_name},
            {"enabled", true},
            {"side", l.side},
            {"id_base", l.base},
            {"position_digit_axes", l.axes},
            {"lattice_bits", bits},
            {"velocity_cell_bits", cells},
            {"id_transform", hilbert_codec},
            {"velocity_order", "coarse_eulerian_morton_then_lagrangian_morton"},
            {"permutation_sidecar", false},
            {"position_source", "decoded_xnyzip_float32"}};
  validate_structured(j);
  return j;
}
Array hilbert_ids(const Array &a, const Json &j, bool inverse,
                  const Type &target) {
  validate_structured(j);
  size_t side = j["side"];
  int base = j["id_base"], bits = j["lattice_bits"];
  auto axes = j["position_digit_axes"].get<std::array<size_t, 3>>();
  auto type = inverse ? target : type_of("uint32");
  Array out{type, std::vector<uint8_t>(checked_bytes(a.size(), type.bytes))};
  for (size_t i = 0; i < a.size(); ++i) {
    if (!inverse) {
      auto digits = coordinates(a.number(i), side, base);
      std::array<size_t, 3> physical{digits[axes[0]], digits[axes[1]],
                                     digits[axes[2]]};
      out.set(i, hilbert(morton(physical), bits, false));
    } else {
      auto code = a.number(i);
      if (code < 0 || code >= uint64_t(1) << (3 * bits))
        throw std::runtime_error("Hilbert code exceeds capacity");
      auto physical =
          unmorton(hilbert(static_cast<uint32_t>(code), bits, true));
      std::array<size_t, 3> digits{};
      for (size_t k = 0; k < 3; ++k) {
        if (physical[k] >= side)
          throw std::runtime_error("Hilbert code outside lattice");
        digits[axes[k]] = physical[k];
      }
      out.set(i, (digits[0] * side + digits[1]) * side + digits[2] + base);
    }
  }
  return out;
}
std::vector<size_t> hybrid_order(const Array &ids, const Triplet &pos,
                                 const Json &j) {
  validate_structured(j);
  size_t side = j["side"], scale = size_t(1)
                                   << j["velocity_cell_bits"].get<int>();
  int base = j["id_base"], bits = j["lattice_bits"];
  auto axes = j["position_digit_axes"].get<std::array<size_t, 3>>();
  std::vector<uint64_t> keys(ids.size());
  std::vector<size_t> order(ids.size());
  std::iota(order.begin(), order.end(), 0);
  for (size_t i = 0; i < ids.size(); ++i) {
    auto c = coordinates(ids.number(i), side, base);
    std::array<size_t, 3> physical{c[axes[0]], c[axes[1]], c[axes[2]]}, cells;
    for (size_t k = 0; k < 3; ++k) {
      double v = pos[k].number(i);
      if (!std::isfinite(v))
        throw std::runtime_error("Nonfinite decoded position");
      cells[k] = std::min(scale - 1, size_t(std::floor(modulo(v, 1) * scale)));
    }
    keys[i] = (uint64_t(morton(cells)) << (3 * bits)) | morton(physical);
  }
  std::stable_sort(order.begin(), order.end(),
                   [&](size_t a, size_t b) { return keys[a] < keys[b]; });
  return order;
}
} // namespace particle
