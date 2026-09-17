#include <SZo/api/sz.hpp>
#include <cstdlib>
#include <cstring>
#include <memory>
extern "C" unsigned char *particle_szo_lorenzo(int dtype, void *data,
                                               size_t count, double bound,
                                               size_t *size) {
  SZo::Config config(count);
  config.errorBoundMode = SZo::EB_ABS;
  config.absErrorBound = bound;
  config.cmprAlgo = SZo::ALGO_LORENZO_REG;
  config.lorenzo = true;
  config.lorenzo2 = false;
  config.regression = false;
  config.regression2 = false;
  char *encoded =
      dtype == 0
          ? SZ_compress<float>(config, static_cast<float *>(data), *size)
          : SZ_compress<double>(config, static_cast<double *>(data), *size);
  std::unique_ptr<char[]> owner(encoded);
  auto *output = static_cast<unsigned char *>(std::malloc(*size));
  if (!output)
    return nullptr;
  std::memcpy(output, encoded, *size);
  return output;
}
