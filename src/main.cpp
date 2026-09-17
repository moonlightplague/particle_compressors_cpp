#include <iostream>

#include "particle/pipeline.hpp"
int main(int argc, char **argv) {
  H5::Exception::dontPrint();
  try {
    return particle::run(particle::parse_cli(argc, argv));
  } catch (int code) {
    return code;
  } catch (const H5::Exception &e) {
    std::cerr << "error: HDF5: " << e.getDetailMsg() << '\n';
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
  }
  return 2;
}
