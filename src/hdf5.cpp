#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <stdexcept>

#include "particle/pipeline.hpp"
namespace particle {
namespace {
bool boolean_type(const H5::DataType &type) {
  if (type.getClass() != H5T_ENUM || type.getSize() != 1 ||
      H5Tget_nmembers(type.getId()) != 2)
    return false;
  for (unsigned i = 0; i < 2; ++i) {
    char *raw = H5Tget_member_name(type.getId(), i);
    if (!raw) return false;
    std::string name(raw);
    H5free_memory(raw);
    uint8_t value = 0;
    if (H5Tget_member_value(type.getId(), i, &value) < 0 ||
        !((name == "FALSE" && value == 0) || (name == "TRUE" && value == 1)))
      return false;
  }
  return true;
}
std::string storage_name(const H5::DataType &t) {
  auto type = type_of(t);
  if (type.bytes > 1 && H5Tget_order(t.getId()) == H5T_ORDER_BE)
    return std::string(">") +
           (type.floating ? "f"
            : type.sign   ? "i"
                          : "u") +
           std::to_string(type.bytes);
  return type.name;
}
H5::DataType storage_type(const std::string &name) {
  auto t = type_of(name).h5();
  if (name.starts_with(">")) {
    H5::DataType copy;
    copy.copy(t);
    if (H5Tset_order(copy.getId(), H5T_ORDER_BE) < 0)
      throw std::runtime_error("Cannot set HDF5 byte order");
    return copy;
  }
  return t;
}
Json shaped(const Json &flat, const std::vector<hsize_t> &dims, size_t depth,
            size_t &offset) {
  if (depth == dims.size()) return flat.at(offset++);
  Json a = Json::array();
  for (hsize_t i = 0; i < dims[depth]; ++i)
    a.push_back(shaped(flat, dims, depth + 1, offset));
  return a;
}
void flatten(const Json &j, Json &a) {
  if (j.is_array())
    for (const auto &v : j) flatten(v, a);
  else
    a.push_back(j);
}
std::vector<hsize_t> dimensions(const H5::DataSpace &s) {
  int rank = s.getSimpleExtentNdims();
  std::vector<hsize_t> d(rank);
  if (rank) s.getSimpleExtentDims(d.data());
  return d;
}
void visit(const H5::Group &g, const std::string &prefix,
           std::map<std::string, std::string> &available) {
  std::vector<std::string> names;
  for (hsize_t i = 0; i < g.getNumObjs(); ++i)
    names.push_back(g.getObjnameByIdx(i));
  std::sort(names.begin(), names.end());
  for (auto &name : names) {
    std::string path = prefix + name;
    H5O_info2_t info{};
    if (H5Oget_info_by_name3(g.getId(), name.c_str(), &info, H5O_INFO_BASIC,
                             H5P_DEFAULT) < 0)
      throw std::runtime_error("Cannot inspect HDF5 object");
    if (info.type == H5O_TYPE_GROUP)
      visit(g.openGroup(name), path + "/", available);
    else if (info.type == H5O_TYPE_DATASET) {
      auto lower = name;
      std::transform(lower.begin(), lower.end(), lower.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      available[lower] = path;
    }
  }
}
}  // namespace
Json attributes(const H5::H5Object &obj) {
  Json result = Json::object();
  for (int i = 0; i < obj.getNumAttrs(); ++i) {
    auto a = obj.openAttribute(i);
    auto dt = a.getDataType();
    auto space = a.getSpace();
    auto dims = dimensions(space);
    size_t n = space.getSimpleExtentNpoints();
    Json flat = Json::array();
    std::string dtype;
    if (dt.getClass() == H5T_STRING) {
      auto st = a.getStrType();
      dtype =
          st.isVariableStr() ? "object" : "|S" + std::to_string(st.getSize());
      if (st.isVariableStr()) {
        std::vector<char *> values(n, nullptr);
        a.read(st, values.data());
        for (auto p : values) flat.push_back(p ? p : "");
        H5Treclaim(st.getId(), space.getId(), H5P_DEFAULT, values.data());
      } else {
        std::vector<char> b(checked_bytes(n, st.getSize()));
        a.read(st, b.data());
        for (size_t k = 0; k < n; ++k) {
          auto p = b.data() + k * st.getSize();
          flat.push_back(std::string(p, strnlen(p, st.getSize())));
        }
      }
    } else if (boolean_type(dt)) {
      std::vector<uint8_t> values(n);
      a.read(dt, values.data());
      dtype = "bool";
      for (auto value : values) flat.push_back(value != 0);
    } else {
      auto t = type_of(dt);
      dtype = storage_name(dt);
      Array b{t, std::vector<uint8_t>(checked_bytes(n, t.bytes))};
      a.read(t.h5(), b.bytes.data());
      for (size_t k = 0; k < n; ++k) {
        if (t.floating)
          flat.push_back(static_cast<double>(b.number(k)));
        else if (t.sign)
          flat.push_back(static_cast<int64_t>(b.number(k)));
        else
          flat.push_back(static_cast<uint64_t>(b.number(k)));
      }
    }
    size_t offset = 0;
    result[a.getName()] = {{"dtype", dtype},
                           {"shape", dims},
                           {"value", shaped(flat, dims, 0, offset)}};
  }
  return result;
}
void apply_attributes(H5::H5Object &obj, const Json &attrs) {
  for (auto it = attrs.begin(); it != attrs.end(); ++it) {
    const auto &v = it.value();
    auto dims = v.value("shape", std::vector<hsize_t>{});
    H5::DataSpace space = dims.empty()
                              ? H5::DataSpace(H5S_SCALAR)
                              : H5::DataSpace(dims.size(), dims.data());
    Json flat = Json::array();
    flatten(v.at("value"), flat);
    auto name = v.at("dtype").get<std::string>();
    if (name == "object" || name.find('S') != std::string::npos ||
        name.find('U') != std::string::npos) {
      bool variable = name == "object" || name.find('U') != std::string::npos;
      size_t width =
          variable ? H5T_VARIABLE : std::stoul(name.substr(name.find('S') + 1));
      H5::StrType t(H5::PredType::C_S1, width);
      t.setCset(H5T_CSET_UTF8);
      auto a = obj.createAttribute(it.key(), t, space);
      if (variable) {
        std::vector<std::string> strings;
        for (auto &x : flat) strings.push_back(x.get<std::string>());
        std::vector<const char *> ptrs;
        for (auto &s : strings) ptrs.push_back(s.c_str());
        a.write(t, ptrs.data());
      } else {
        std::vector<char> b(checked_bytes(flat.size(), width), 0);
        for (size_t i = 0; i < flat.size(); ++i) {
          auto s = flat[i].get<std::string>();
          std::memcpy(b.data() + i * width, s.data(),
                      std::min(width, s.size()));
        }
        a.write(t, b.data());
      }
    } else if (name == "bool") {
      H5::IntType base(H5::PredType::NATIVE_INT8);
      H5::EnumType type(base);
      int8_t no = 0, yes = 1;
      type.insert("FALSE", &no);
      type.insert("TRUE", &yes);
      std::vector<int8_t> values;
      for (const auto &value : flat) values.push_back(value.get<bool>());
      auto attribute = obj.createAttribute(it.key(), type, space);
      attribute.write(type, values.data());
    } else {
      auto t = type_of(name);
      Array a{t, std::vector<uint8_t>(checked_bytes(flat.size(), t.bytes))};
      for (size_t i = 0; i < flat.size(); ++i) {
        long double x = flat[i].is_number_unsigned()
                            ? static_cast<long double>(flat[i].get<uint64_t>())
                        : flat[i].is_number_integer()
                            ? static_cast<long double>(flat[i].get<int64_t>())
                            : flat[i].get<double>();
        a.set(i, x);
      }
      auto attr = obj.createAttribute(it.key(), storage_type(name), space);
      attr.write(t.h5(), a.bytes.data());
    }
  }
}
Snapshot read_hdf5(const fs::path &path, int64_t limit) {
  if (H5Fis_hdf5(path.c_str()) <= 0) return read_native(path, limit);
  H5::H5File file(path.string(), H5F_ACC_RDONLY);
  std::map<std::string, std::string> available;
  visit(file.openGroup("/"), "", available);
  static const std::array<std::array<const char *, 3>, 7> aliases{
      {{{"id", "particle_id", "pid"}},
       {{"x", "posx", "position_x"}},
       {{"y", "posy", "position_y"}},
       {{"z", "posz", "position_z"}},
       {{"vx", "velx", "velocity_x"}},
       {{"vy", "vely", "velocity_y"}},
       {{"vz", "velz", "velocity_z"}}}};
  Snapshot s;
  size_t total = 0, count = 0;
  s.metadata = {{"fields", Json::object()}, {"root_attrs", attributes(file)}};
  for (size_t i = 0; i < 7; ++i) {
    std::string name;
    for (auto alias : aliases[i])
      if (available.contains(alias)) {
        name = available.at(alias);
        break;
      }
    if (name.empty())
      throw std::runtime_error("Could not find dataset for logical field " +
                               std::string(fields[i]));
    auto ds = file.openDataSet(name);
    auto dims = dimensions(ds.getSpace());
    if (dims.size() != 1)
      throw std::runtime_error("Particle fields must be one-dimensional");
    if (i == 0) {
      total = dims[0];
      count = limit > 0 ? std::min(total, static_cast<size_t>(limit)) : total;
    } else if (dims[0] != total)
      throw std::runtime_error("Particle fields do not have the same length");
    auto t = type_of(ds.getDataType());
    if (i == 0 && t.floating)
      throw std::runtime_error("Particle IDs must have an integer dtype");
    if (i >= 4 && !t.floating)
      throw std::runtime_error("Velocity fields must be float32/float64");
    s.data[i] = {t, std::vector<uint8_t>(checked_bytes(count, t.bytes))};
    if (count) {
      auto space = ds.getSpace();
      hsize_t start = 0, n = count;
      space.selectHyperslab(H5S_SELECT_SET, &n, &start);
      H5::DataSpace mem(1, &n);
      ds.read(s.data[i].bytes.data(), t.h5(), mem, space);
    }
    s.metadata["fields"][fields[i]] = {
        {"h5_path", name},
        {"dtype", storage_name(ds.getDataType())},
        {"shape", dims},
        {"selected_shape", {count}},
        {"attrs", attributes(ds)}};
  }
  if (limit > 0 && s.metadata["root_attrs"].contains("npart")) {
    auto &a = s.metadata["root_attrs"]["npart"];
    a["shape"] = Json::array();
    a["value"] = count;
  }
  return s;
}
void write_hdf5(const fs::path &path, const Snapshot &s, bool force) {
  require_output(path, force);
  H5::H5File file(path.string(), H5F_ACC_TRUNC);
  apply_attributes(file, s.metadata.value("root_attrs", Json::object()));
  hsize_t n = s.count();
  H5::DataSpace space(1, &n);
  for (size_t i = 0; i < 7; ++i) {
    auto f = s.metadata.at("fields").at(fields[i]);
    std::string name = f.at("h5_path");
    size_t pos = 0;
    while ((pos = name.find('/', pos)) != std::string::npos) {
      auto group = name.substr(0, pos++);
      if (!group.empty() &&
          H5Lexists(file.getId(), group.c_str(), H5P_DEFAULT) <= 0)
        file.createGroup(group);
    }
    auto ds = file.createDataSet(name, storage_type(f.at("dtype")), space);
    if (n) ds.write(s.data[i].bytes.data(), s.data[i].type.h5());
    apply_attributes(ds, f.value("attrs", Json::object()));
  }
}
}  // namespace particle
