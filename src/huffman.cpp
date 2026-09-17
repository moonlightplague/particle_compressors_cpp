#include "particle/huffman.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
namespace particle {
namespace {
void put(std::vector<uint8_t> &out, uint64_t value, size_t n) {
  for (size_t i = 0; i < n; ++i)
    out.push_back(value >> (8 * i));
}
uint64_t get(const std::vector<uint8_t> &in, size_t &p, size_t n) {
  if (p > in.size() || n > in.size() - p)
    throw std::runtime_error("Truncated Huffman stream");
  uint64_t v = 0;
  for (size_t i = 0; i < n; ++i)
    v |= uint64_t(in[p++]) << (i * 8);
  return v;
}
std::map<uint64_t, std::string>
canonical(const std::map<uint64_t, size_t> &lengths) {
  std::vector<std::pair<size_t, uint64_t>> order;
  for (auto [v, n] : lengths) {
    if (!n || n > 255)
      throw std::runtime_error("Invalid Huffman length");
    order.emplace_back(n, v);
  }
  std::sort(order.begin(), order.end());
  std::map<uint64_t, std::string> codes;
  std::string bits;
  bool overflow = false;
  for (auto [n, v] : order) {
    if (overflow || n < bits.size())
      throw std::runtime_error("Oversubscribed Huffman table");
    bits.resize(n, '0');
    codes[v] = bits;
    size_t i = bits.size();
    while (i && bits[i - 1] == '1')
      bits[--i] = '0';
    if (i)
      bits[i - 1] = '1';
    else
      overflow = true;
  }
  return codes;
}
} // namespace
std::vector<uint8_t> huffman_encode(const Array &a, Json &meta) {
  if (a.type.name != "uint32" && a.type.name != "uint64")
    throw std::runtime_error("Huffman requires uint32/uint64");
  uint64_t mask = a.type.bytes == 4 ? UINT32_MAX : UINT64_MAX, previous = 0;
  std::vector<uint64_t> delta;
  std::map<uint64_t, size_t> frequencies;
  for (size_t i = 0; i < a.size(); ++i) {
    uint64_t v = a.number(i), d = (v - previous) & mask;
    previous = v;
    delta.push_back(d);
    ++frequencies[d];
  }
  struct Node {
    size_t frequency;
    uint64_t minimum;
    size_t serial;
    int left = -1, right = -1;
  };
  std::vector<Node> nodes;
  using Item = std::tuple<size_t, uint64_t, size_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
  for (auto [v, n] : frequencies) {
    size_t index = nodes.size();
    nodes.push_back({n, v, index});
    heap.emplace(n, v, index);
  }
  while (heap.size() > 1) {
    auto [a_freq, a_min, a_index] = heap.top();
    heap.pop();
    auto [b_freq, b_min, b_index] = heap.top();
    heap.pop();
    size_t index = nodes.size();
    nodes.push_back({a_freq + b_freq, std::min(a_min, b_min), index,
                     int(a_index), int(b_index)});
    heap.emplace(a_freq + b_freq, std::min(a_min, b_min), index);
  }
  std::map<uint64_t, size_t> lengths;
  std::vector<std::pair<size_t, size_t>> stack;
  if (!heap.empty())
    stack.emplace_back(std::get<2>(heap.top()), 0);
  while (!stack.empty()) {
    auto [index, depth] = stack.back();
    stack.pop_back();
    auto &n = nodes[index];
    if (n.left < 0)
      lengths[n.minimum] = std::max<size_t>(depth, 1);
    else {
      stack.emplace_back(n.right, depth + 1);
      stack.emplace_back(n.left, depth + 1);
    }
  }
  auto codes = canonical(lengths);
  size_t bits = 0;
  for (auto [v, n] : frequencies)
    bits += checked_bytes(n, lengths[v]);
  std::vector<uint8_t> out{'L', 'C', 'P', 'H', 'U', 'F', '3', 0};
  put(out, a.size(), 8);
  put(out, bits, 8);
  put(out, lengths.size(), 4);
  put(out, a.type.bytes * 8, 1);
  put(out, 0, 3);
  for (auto [v, n] : lengths) {
    put(out, v, 8);
    put(out, n, 1);
  }
  size_t start = out.size();
  out.resize(start + bits / 8 + (bits % 8 != 0), 0);
  size_t bit = 0;
  for (auto v : delta)
    for (char c : codes[v]) {
      out[start + bit / 8] |= (c == '1') << (7 - bit % 8);
      ++bit;
    }
  meta = {{"codec", "canonical_huffman"},
          {"symbol_dtype", a.type.name},
          {"value_count", a.size()},
          {"unique_symbol_count", lengths.size()},
          {"encoded_bit_count", bits},
          {"encoded_bytes", out.size()},
          {"transform", a.type.name + "_delta_modulo"},
          {"format", "LCPHUF3"}};
  return out;
}
Array huffman_decode(const std::vector<uint8_t> &b, size_t count) {
  if (b.size() < 8)
    throw std::runtime_error("Truncated Huffman header");
  std::string magic(b.begin(), b.begin() + 8);
  bool v3 = magic == std::string("LCPHUF3\0", 8);
  if (!v3 && magic != std::string("LCPHUF2\0", 8))
    throw std::runtime_error("Invalid Huffman magic");
  size_t p = 8, n = get(b, p, 8), bits = get(b, p, 8), symbols = get(b, p, 4),
         width = 32;
  if (v3) {
    width = get(b, p, 1);
    get(b, p, 3);
  }
  if (n != count || symbols > n || (width != 32 && width != 64))
    throw std::runtime_error("Invalid Huffman header");
  std::map<uint64_t, size_t> lengths;
  for (size_t i = 0; i < symbols; ++i) {
    uint64_t v = get(b, p, v3 ? 8 : 4);
    size_t length = get(b, p, 1);
    if (lengths.contains(v) || (width == 32 && v > UINT32_MAX))
      throw std::runtime_error("Invalid Huffman symbol");
    lengths[v] = length;
  }
  if (bool(n) != bool(symbols) || b.size() - p != bits / 8 + (bits % 8 != 0))
    throw std::runtime_error("Huffman payload size mismatch");
  auto codes = canonical(lengths);
  struct Node {
    int next[2]{-1, -1};
    bool leaf = false;
    uint64_t value = 0;
  };
  std::vector<Node> trie(1);
  for (auto &[value, code] : codes) {
    size_t node = 0;
    for (char ch : code) {
      int bit = ch - '0';
      if (trie[node].leaf)
        throw std::runtime_error("Invalid Huffman prefix");
      if (trie[node].next[bit] < 0) {
        trie[node].next[bit] = trie.size();
        trie.emplace_back();
      }
      node = trie[node].next[bit];
    }
    trie[node].leaf = true;
    trie[node].value = value;
  }
  auto type = type_of(width == 32 ? "uint32" : "uint64");
  Array out{type, std::vector<uint8_t>(checked_bytes(count, type.bytes))};
  size_t node = 0, written = 0;
  uint64_t previous = 0, mask = width == 32 ? UINT32_MAX : UINT64_MAX;
  for (size_t i = 0; i < bits; ++i) {
    int bit = (b[p + i / 8] >> (7 - i % 8)) & 1;
    int next = trie[node].next[bit];
    if (next < 0)
      throw std::runtime_error("Invalid Huffman code");
    node = next;
    if (trie[node].leaf) {
      if (written >= count)
        throw std::runtime_error("Excess Huffman values");
      previous = (previous + trie[node].value) & mask;
      out.set(written++, previous);
      node = 0;
    }
  }
  if (written != count || node)
    throw std::runtime_error("Incomplete Huffman data");
  return out;
}
} // namespace particle
