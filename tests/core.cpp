#include "particle/huffman.hpp"
#include "particle/layout.hpp"
#include "particle/pipeline.hpp"
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace particle;
void check(bool b, const char *message) {
  if (!b)
    throw std::runtime_error(message);
}
int main() {
  try {
    for (auto name : {"uint8", "int8", "uint16", "int16", "uint32", "int32",
                      "uint64", "int64", "float32", "float64"}) {
      auto t = type_of(name);
      Array a{t, std::vector<uint8_t>(31 * t.bytes)};
      for (size_t i = 0; i < a.size(); ++i)
        a.set(i, t.floating ? double(i) / 7 : i);
      size_t n;
      auto b = encode(a, "pcodec", 0, n);
      check(decode(b, t, a.size(), n, "pcodec").bytes == a.bytes,
            "pcodec roundtrip mismatch");
    }
    Array ids{type_of("uint64"), std::vector<uint8_t>(24)};
    ids.set(0, 9007199254740995.L);
    ids.set(1, 9007199254740993.L);
    ids.set(2, 9007199254740995.L);
    check(stable_id_order(ids) == std::vector<size_t>({1, 0, 2}),
          "stable wide-ID sorting failed");
    for (auto codec : {"sz3", "szo"}) {
      Array a{type_of("float32"), std::vector<uint8_t>(100 * 4)};
      for (size_t i = 0; i < a.size(); ++i)
        a.set(i, i * .123);
      size_t n;
      auto b = encode(a, codec, .001, n);
      auto d = decode(b, a.type, a.size(), n, codec);
      for (size_t i = 0; i < a.size(); ++i)
        check(std::abs(a.number(i) - d.number(i)) <= .00101,
              "lossy bound exceeded");
    }
    for (auto name : {"uint32", "uint64"}) {
      auto type = type_of(name);
      Array a{type, std::vector<uint8_t>(6 * type.bytes)};
      long double max = type.bytes == 4 ? static_cast<long double>(UINT32_MAX)
                                        : static_cast<long double>(UINT64_MAX);
      std::array<long double, 6> values{0, max, 5, 5, 0, 13};
      for (size_t i = 0; i < 6; ++i)
        a.set(i, values[i]);
      Json meta;
      auto bytes = huffman_encode(a, meta);
      auto golden =
          type.bytes == 4
              ? "4c4350485546330006000000000000000e0000000000000005000000200000"
                "000000000000000000020600000000000000030d0000000000000003fbffff"
                "ff0000000002ffffffff00000000022c3c"
              : "4c4350485546330006000000000000000e0000000000000005000000400000"
                "000000000000000000020600000000000000030d0000000000000003fbffff"
                "ffffffffff02ffffffffffffffff022c3c";
      std::vector<uint8_t> expected;
      std::string text = golden;
      for (size_t i = 0; i < text.size(); i += 2)
        expected.push_back(std::stoul(text.substr(i, 2), nullptr, 16));
      check(bytes == expected, "Huffman differs from Python golden");
      check(huffman_decode(bytes, 6).bytes == a.bytes,
            "Huffman modulo delta mismatch");
      bytes.pop_back();
      bool rejected = false;
      try {
        huffman_decode(bytes, 6);
      } catch (const std::exception &) {
        rejected = true;
      }
      check(rejected, "Truncated Huffman accepted");
    }
    Array empty{type_of("uint32"), {}};
    size_t encoded;
    auto stream = encode(empty, "pcodec", 0, encoded);
    check(decode(stream, empty.type, 0, 0, "pcodec").size() == 0,
          "Empty packed order failed");
    Array special{type_of("float64"), std::vector<uint8_t>(32)};
    std::array<uint64_t, 4> bits{0, 0x8000000000000000ULL,
                                 0x7ff8000000001234ULL, 0x7ff0000000000000ULL};
    std::memcpy(special.bytes.data(), bits.data(), 32);
    stream = encode(special, "pcodec", 0, encoded);
    check(decode(stream, special.type, 4, 4, "pcodec").bytes == special.bytes,
          "Lossless float bits changed");
    std::cout << "core checks passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
