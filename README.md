# Particle compressors — C++ port

C++20 implementation of the production particle pipeline in `origin/`, with a
reusable `particle_pipeline` library and the `particle_compressors` executable.
The reference is a read-only symlink; builds and synthetic fixtures stay here.
See [translation status](docs/translation-status.md) for verified coverage and
compatibility differences and validation limits.

## Build

Dependencies: CMake 3.20+, a C++20 compiler, Rust/Cargo, HDF5 with its C++ API,
Zstandard, OpenMP, fmt, Eigen3, and TBB. Initialize the five repository
submodules before building. yaml-cpp is built from its submodule.

```sh
git submodule update --init --recursive
bash build_hdf5.sh  # if $HOME/.local/hdf5 is not already installed
./build.sh
```

`HDF5_ROOT` selects another HDF5 installation; `BUILD_JOBS` controls compilation
parallelism (default 2). fmt and Eigen can be installed system-wide or supplied
through `PARTICLE_FMT_INCLUDE`, `PARTICLE_FMT_LIBRARY`, and
`PARTICLE_EIGEN_INCLUDE` CMake variables. Existing Conan cache installations are
also discovered. `build.sh` recognizes a workspace-local Rust installation in
`.deps/cargo` and `.deps/rustup`, or uses Cargo on PATH. It does not install Rust.

LCP and SZO use the repository's native code. pcodec is compiled as a Rust static
library through its C ABI. The installed Python `pysz` reference uses SZ3 3.3.2;
LCP bundles 3.1.8. CMake downloads a pinned, hash-verified SZ3 3.3.2 source archive
for the fieldwise C adapter to preserve Python stream interoperability.
LCP's compression algorithm remains unchanged. A build-local header overlay fixes
an upstream radix-sort allocator mismatch found by ASan. Native symbol names are
isolated during compilation.

```sh
./build/particle_compressors roundtrip input.h5 --work-dir runs/example \
  --pos-compressor szo --vel-compressor szo --metrics
./build/particle_compressors preprocess input.h5 --work-dir runs/raw
./build/particle_compressors compress input.h5 --work-dir runs/package \
  --pos-compressor lcp --vel-compressor lcp --vel-chunk-size 100000
./build/particle_compressors decompress --work-dir runs/package --force
```

Use `--config path.yaml` for advanced defaults and `--help` on any command for
options. Relative tool paths in configuration are resolved from that file.
Lossless positions and velocities use `--pos-compressor pcodec
--vel-compressor pcodec`. Lossy position codecs establish their documented
canonical order; IDs remain attached to the same particle. `--sort` applies a
stable ID sort to fieldwise pipelines. Native `dat_*`/`cfg_*` inputs, directory
batches, and `--merge` are supported. Merge rejects duplicate IDs and incompatible
schemas. `--clean-raw` removes intermediate raw directories after reconstruction.

LCP and XnYZip are called through linked native APIs. Their historic executable
path options remain parseable and recorded, but do not select alternate codec
binaries. Packages use upstream formats and require matching native platforms;
legacy native codec payloads should come from trusted sources.

## Test

The default build runs native core checks. Enable the differential suite with
a Python environment containing the reference dependencies:

```sh
./build.sh -DPARTICLE_REFERENCE_PYTHON=/path/to/python
ctest --test-dir build --output-on-failure
```

These tests generate small deterministic inputs in the build directory and
exercise Python/C++ encoding and decoding, ordering, layouts, native inputs,
merging, batch failures, and cleanup. They do not run large datasets.

The application holds selected fields and merged inputs in memory. Size large
jobs accordingly. The status document records the small benchmark and validation
limits; large snapshot streaming was not benchmarked.

For application sanitizer checks:

```sh
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug \
  -DHDF5_ROOT="$HOME/.local/hdf5" -DPARTICLE_SANITIZERS=ON \
  -DPARTICLE_REFERENCE_PYTHON=/path/to/python
cmake --build build-sanitize --parallel 2
ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build-sanitize --output-on-failure
```

If using workspace-local Rust, export CARGO_HOME, RUSTUP_HOME and PATH as in
`build.sh` before configuring a separate build directory. The six suites cover
native core checks plus differential field, triplet, layout, runtime and contract
checks. All 96 available reference unit tests also passed for the recorded revision.
