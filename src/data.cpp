#include "particle/pipeline.hpp"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace particle {
Type type_of(const std::string &n) {
  if (n.size() == 3 &&
      (n[0] == '>' || n[0] == '<' || n[0] == '=' || n[0] == '|') &&
      (n[1] == 'i' || n[1] == 'u' || n[1] == 'f') && n[2] >= '1' && n[2] <= '8')
    return type_of(std::string(n[1] == 'i'   ? "int"
                               : n[1] == 'u' ? "uint"
                                             : "float") +
                   std::to_string((n[2] - '0') * 8));
#define T(N, B, P, F, S, H)                                                    \
  if (n == N)                                                                  \
  return {N, B, P, F, S}
  T("uint8", 1, 10, false, false, NATIVE_UINT8);
  T("int8", 1, 11, false, true, NATIVE_INT8);
  T("uint16", 2, 7, false, false, NATIVE_UINT16);
  T("int16", 2, 8, false, true, NATIVE_INT16);
  T("uint32", 4, 1, false, false, NATIVE_UINT32);
  T("int32", 4, 3, false, true, NATIVE_INT32);
  T("uint64", 8, 2, false, false, NATIVE_UINT64);
  T("int64", 8, 4, false, true, NATIVE_INT64);
  T("float32", 4, 5, true, true, NATIVE_FLOAT);
  T("float64", 8, 6, true, true, NATIVE_DOUBLE);
#undef T
  throw std::runtime_error("Unsupported numeric dtype: " + n);
}
H5::DataType Type::h5() const {
  const auto &n = name;
#define T(N, B, P, F, S, H)                                                    \
  if (n == N)                                                                  \
  return H5::PredType::H
  T("uint8", 1, 10, false, false, NATIVE_UINT8);
  T("int8", 1, 11, false, true, NATIVE_INT8);
  T("uint16", 2, 7, false, false, NATIVE_UINT16);
  T("int16", 2, 8, false, true, NATIVE_INT16);
  T("uint32", 4, 1, false, false, NATIVE_UINT32);
  T("int32", 4, 3, false, true, NATIVE_INT32);
  T("uint64", 8, 2, false, false, NATIVE_UINT64);
  T("int64", 8, 4, false, true, NATIVE_INT64);
  T("float32", 4, 5, true, true, NATIVE_FLOAT);
  T("float64", 8, 6, true, true, NATIVE_DOUBLE);
#undef T
  throw std::runtime_error("Unsupported numeric dtype: " + n);
}
Type type_of(const H5::DataType &t) {
  auto cls = t.getClass();
  if (cls == H5T_ENUM)
    return type_of(t.getSuper());
  if (cls != H5T_INTEGER && cls != H5T_FLOAT)
    throw std::runtime_error("Expected integer or floating-point dataset");
  bool sign = cls == H5T_INTEGER && H5Tget_sign(t.getId()) == H5T_SGN_2;
  return type_of(std::string(cls == H5T_FLOAT ? "float"
                             : sign           ? "int"
                                              : "uint") +
                 std::to_string(t.getSize() * 8));
}
size_t checked_bytes(size_t n, size_t w) {
  if (w && n > std::numeric_limits<size_t>::max() / w)
    throw std::runtime_error("Array size overflow");
  return n * w;
}
template <class T> long double get(const uint8_t *p) {
  T v;
  std::memcpy(&v, p, sizeof v);
  return v;
}
template <class T> void put(uint8_t *p, long double v) {
  T value = static_cast<T>(v);
  std::memcpy(p, &value, sizeof value);
}
long double Array::number(size_t i) const {
  if (i >= size())
    throw std::out_of_range("Array index");
  auto p = bytes.data() + i * type.bytes;
#define G(N, T)                                                                \
  if (type.name == N)                                                          \
  return get<T>(p)
  G("uint8", uint8_t);
  G("int8", int8_t);
  G("uint16", uint16_t);
  G("int16", int16_t);
  G("uint32", uint32_t);
  G("int32", int32_t);
  G("uint64", uint64_t);
  G("int64", int64_t);
  G("float32", float);
  G("float64", double);
#undef G
  throw std::runtime_error("Unsupported dtype");
}
void Array::set(size_t i, long double v) {
  if (i >= size())
    throw std::out_of_range("Array index");
  auto p = bytes.data() + i * type.bytes;
#define P(N, T)                                                                \
  if (type.name == N) {                                                        \
    put<T>(p, v);                                                              \
    return;                                                                    \
  }
  P("uint8", uint8_t);
  P("int8", int8_t);
  P("uint16", uint16_t);
  P("int16", int16_t);
  P("uint32", uint32_t);
  P("int32", int32_t);
  P("uint64", uint64_t);
  P("int64", int64_t);
  P("float32", float);
  P("float64", double);
#undef P
  throw std::runtime_error("Unsupported dtype");
}
Array Array::convert(const Type &t) const {
  if (t.name == type.name)
    return *this;
  Array out{t, std::vector<uint8_t>(checked_bytes(size(), t.bytes))};
  for (size_t i = 0; i < size(); ++i)
    out.set(i, number(i));
  return out;
}
Array Array::gather(const std::vector<size_t> &order) const {
  Array out{type,
            std::vector<uint8_t>(checked_bytes(order.size(), type.bytes))};
  for (size_t i = 0; i < order.size(); ++i) {
    if (order[i] >= size())
      throw std::runtime_error("Invalid permutation");
    std::memcpy(out.bytes.data() + i * type.bytes,
                bytes.data() + order[i] * type.bytes, type.bytes);
  }
  return out;
}
std::vector<size_t> stable_id_order(const Array &ids) {
  std::vector<size_t> order(ids.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return ids.number(a) < ids.number(b);
  });
  return order;
}
std::vector<uint8_t> read_bytes(const fs::path &p) {
  std::ifstream f(p, std::ios::binary | std::ios::ate);
  if (!f)
    throw std::runtime_error("Cannot read " + p.string());
  auto n = f.tellg();
  if (n < 0)
    throw std::runtime_error("Cannot size " + p.string());
  std::vector<uint8_t> b(static_cast<size_t>(n));
  f.seekg(0);
  if (!f.read(reinterpret_cast<char *>(b.data()), b.size()))
    throw std::runtime_error("Truncated read: " + p.string());
  return b;
}
void require_output(const fs::path &p, bool force) {
  if (fs::exists(p) && !force)
    throw std::runtime_error(p.string() +
                             " already exists. Use --force to overwrite.");
}
void write_bytes(const fs::path &p, const std::vector<uint8_t> &b, bool force) {
  require_output(p, force);
  fs::create_directories(p.parent_path());
  std::ofstream f(p, std::ios::binary);
  if (!f.write(reinterpret_cast<const char *>(b.data()), b.size()))
    throw std::runtime_error("Cannot write " + p.string());
}
} // namespace particle
