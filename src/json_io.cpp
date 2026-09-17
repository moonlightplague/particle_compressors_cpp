#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "particle/pipeline.hpp"
namespace particle {
namespace {
// Python's reports use JSON's commonly implemented NaN/Infinity extension.
void render(std::ostream &out, const Json &j, size_t depth) {
  auto indent = [&](size_t n) { out << std::string(n * 2, ' '); };
  if (j.is_number_float() && !std::isfinite(j.get<double>())) {
    double v = j.get<double>();
    out << (std::isnan(v) ? "NaN" : v < 0 ? "-Infinity" : "Infinity");
  } else if (j.is_object()) {
    out << '{';
    std::vector<std::string> keys;
    for (auto it = j.begin(); it != j.end(); ++it) keys.push_back(it.key());
    std::sort(keys.begin(), keys.end());
    for (size_t i = 0; i < keys.size(); ++i) {
      out << (i ? ",\n" : "\n");
      indent(depth + 1);
      out << Json(keys[i]).dump() << ": ";
      render(out, j.at(keys[i]), depth + 1);
    }
    if (!keys.empty()) {
      out << '\n';
      indent(depth);
    }
    out << '}';
  } else if (j.is_array()) {
    out << '[';
    for (size_t i = 0; i < j.size(); ++i) {
      out << (i ? ",\n" : "\n");
      indent(depth + 1);
      render(out, j[i], depth + 1);
    }
    if (!j.empty()) {
      out << '\n';
      indent(depth);
    }
    out << ']';
  } else
    out << j.dump();
}
class Reader {
  const std::string &text;
  size_t pos = 0;
  void space() {
    while (pos < text.size() &&
           std::isspace(static_cast<unsigned char>(text[pos])))
      ++pos;
  }
  char peek() {
    space();
    if (pos == text.size()) throw std::runtime_error("Truncated JSON");
    return text[pos];
  }
  void expect(char c) {
    if (peek() != c) throw std::runtime_error("Malformed JSON");
    ++pos;
  }
  Json string() {
    size_t begin = pos;
    expect('"');
    bool escape = false;
    while (pos < text.size()) {
      char c = text[pos++];
      if (!escape && c == '"')
        return Json::parse(text.substr(begin, pos - begin));
      if (!escape && c == '\\')
        escape = true;
      else
        escape = false;
    }
    throw std::runtime_error("Unterminated JSON string");
  }
  Json value(size_t depth) {
    if (depth > 512) throw std::runtime_error("JSON nesting exceeds 512");
    char c = peek();
    if (c == '"') return string();
    if (c == '{' || c == '[') {
      ++pos;
      Json j = c == '{' ? Json::object() : Json::array();
      char close = c == '{' ? '}' : ']';
      if (peek() == close) {
        ++pos;
        return j;
      }
      while (true) {
        if (c == '{') {
          if (peek() != '"')
            throw std::runtime_error("JSON key is not a string");
          std::string key = string();
          expect(':');
          j[key] = value(depth + 1);
        } else
          j.push_back(value(depth + 1));
        if (peek() == close) {
          ++pos;
          break;
        }
        expect(',');
      }
      return j;
    }
    size_t begin = pos;
    while (pos < text.size() && text[pos] != ',' && text[pos] != ']' &&
           text[pos] != '}' &&
           !std::isspace(static_cast<unsigned char>(text[pos])))
      ++pos;
    auto token = text.substr(begin, pos - begin);
    if (token == "Infinity") return std::numeric_limits<double>::infinity();
    if (token == "-Infinity") return -std::numeric_limits<double>::infinity();
    if (token == "NaN") return std::numeric_limits<double>::quiet_NaN();
    return Json::parse(token);
  }

 public:
  explicit Reader(const std::string &s) : text(s) {}
  Json read() {
    Json result = value(0);
    space();
    if (pos != text.size()) throw std::runtime_error("Trailing JSON data");
    return result;
  }
};
}  // namespace
std::string json_text(const Json &j) {
  std::ostringstream out;
  render(out, j, 0);
  out << '\n';
  return out.str();
}
Json read_json(const fs::path &p) {
  std::ifstream f(p);
  if (!f) throw std::runtime_error("Cannot read JSON: " + p.string());
  std::string text((std::istreambuf_iterator<char>(f)), {});
  return Reader(text).read();
}
void write_json(const fs::path &p, const Json &j) {
  fs::create_directories(p.parent_path());
  std::ofstream f(p);
  f << json_text(j);
  if (!f) throw std::runtime_error("Cannot write " + p.string());
}
}  // namespace particle
