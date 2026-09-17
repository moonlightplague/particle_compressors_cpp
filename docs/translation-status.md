# Translation status

The production C++ CLI backbone is implemented. The matrix records verified
coverage; compatibility differences and validation limits are listed below.

## Baseline (2026-09-17)

Compatibility target: preserve all four production commands, their options and
semantic metadata; decode Python-accepted packages and emit Python-decodable
packages. No new package version is introduced. `origin/` remains read-only.

Python reference: `3dfb753bf03ee70d377f4d901f25a6bf04d0b09c`.
Its working tree has pre-existing submodule pointer differences for LCP,
XnYZip, and pcodec. The current C++ repository started clean.

| Dependency | Checked-out revision |
| --- | --- |
| LCP | 576dcecb77b4074e69bbead9100f651d81c4d387 |
| SZo | 4694f8b93044233c6f9ecc5532d80a461f16e120 |
| XnYZip | dad0daea4c920f1fde7fa06b17683cfb38488eac |
| pcodec | 96443fb091b8e669e2ff63273b3bc58676d89bd2 |
| yaml-cpp | 8eb618e889c003457d73f7c54777603835677ff4 |

HDF5 2.1.1 is installed in `$HOME/.local/hdf5` using the repository's
`build_hdf5.sh`. C++20/CMake, CLI11, ordered nlohmann JSON, yaml-cpp, and
caller-owned pcodec C buffers form the foundation. SZ3 and SZO C entry points
are renamed at compilation to avoid their identical exported symbols. Rust
was absent; a workspace-local toolchain is installed under ignored `.deps/`.
No codec algorithm or submodule revision is modified.

## Package baseline

Preprocessing emits version 2. Compression emits 3 (ordinary), 4 (LCP chunks),
5 (XnYZip), 6 (XnYZip chunks), 7 (blockwise order), 8 (lattice), or 9
(structured). Python selects decoding primarily by field metadata rather than
version alone. Legacy fallbacks include int64 order dtype, inferred codec names,
and absent encoded_count. Missing compressed_fields is rejected. Order sidecars
are exact-width integer streams; chunk headers and entries use little-endian
uint64 lengths. Native triplet streams retain upstream platform restrictions.

## Feature and acceptance matrix

| Feature | C++ component | Status | Acceptance |
| --- | --- | --- | --- |
| CMake, library, CLI executable | CMakeLists.txt, cmake/, src/main.cpp | verified | Release build and CTest |
| Four commands, options and configuration | src/cli.cpp, src/pipeline.cpp | verified | default/config/CLI precedence, scaling, limits, invalid combinations; preprocessing JSON parity |
| HDF5 aliases, nested paths, numeric/boolean/fixed string attributes | src/hdf5.cpp | verified | cross-decoded paths, dtypes, attributes including big-endian storage |
| Numeric widths and float bits | src/data.cpp, src/codecs.cpp | verified | ten pcodec dtypes; signed zero, NaN payload, infinity |
| pcodec, SZ3/SZO fields, padding | src/codecs.cpp | verified | 24 encode/decode directions with and without sorting |
| Scaling, field error selection, stable sorting | src/pipeline.cpp | verified | int32 scaling, float32/64 velocities, IDs >2^53, duplicate-ID core case |
| LCP/XnYZip canonical triplets | src/triplets.cpp | verified | canonical IDs and scalar/L2 budgets in four directions |
| Velocity chunk containers | src/triplets.cpp | verified | 47-value chunks, short tail, workers 1/3, corrupt framing |
| Blockwise delta-Huffman order | src/triplets.cpp, src/huffman.cpp | verified | four directions and Python uint32/64 binary goldens |
| Dense/sparse/velocity-only lattice | src/layout.cpp, src/layout_pipeline.cpp, src/shaped.cpp | verified | periodic axis permutation, shaped and LCP positions, four directions |
| Structured Hilbert IDs and hybrid velocity order | src/layout.cpp, src/pipeline.cpp | verified | dense/sparse, fieldwise and XnYZip velocities, four directions |
| Native cfg/dat input | src/native.cpp | verified | field-major data and attributes, cross-decoding |
| Directory discovery, merge, batch failures | src/batch.cpp | verified | native merge, duplicate rejection, mixed valid/invalid batch |
| Worker execution | include/particle/runtime.hpp, src/batch.cpp | verified | serial/parallel fields/chunks/files; HDF5 confined to serial I/O |
| Metrics and batch aggregation | src/metrics.cpp, src/batch.cpp | verified | numerical metrics and batch JSON compared to Python on identical reports |
| Force, cleanup, saved package decode | src/pipeline.cpp | verified | no raw intermediates during cross-decode; overwrite rejection |
| Legacy metadata, invalid-input edges | src/pipeline.cpp | verified | absent compressor/integer codec metadata, LCP artifact inference, empty/nonfinite rejection |
| Performance and sanitizers | build-sanitize/, tests/ | verified | ASan/UBSan application checks and small repeated benchmark |

“Verified” describes the listed checks, not exhaustive compatibility for all
possible inputs. Experiments and historical analysis are outside scope.

## Dependency and acceptance findings

- Python pysz 1.0.3 writes SZ3 3.3.2 streams. LCP's vendored 3.1.8 decoder is
  incompatible. Fieldwise SZ3 uses pinned upstream commit
  `0ebe6fa26b0aa2e274d13b78df81a81e2db1b268`, downloaded and hash-verified by
  CMake. LCP continues to use its own headers.
- AddressSanitizer found an upstream LCP `new[]`/`free` mismatch in
  `SZ3/utils/RadixSort.hpp`. CMake creates a build-local header overlay with
  `delete[]`; sorting and stream algorithms remain unchanged. No submodule
  file is edited. Sanitizer revalidation passes.
- HDF5 objects are never constructed or copied by codec worker threads;
  numeric type descriptors are plain values and I/O uses serial RAII handles.
- Python JSON's Infinity, -Infinity and NaN extensions are accepted and emitted;
  ordinary nlohmann serialization would incorrectly turn these into null.
- Python's variable-length scalar Unicode attribute reconstruction fails in
  this reference environment. C++ handles it as a variable-length string;
  four-direction attribute fixtures therefore use fixed byte strings.
- LCP/XnYZip executable path options are parsed and recorded; linked native APIs
  replace subprocess codecs, so these options do not select another executable.
  Missing historical executables are not an error for linked native codecs.
- Selected fields and merged inputs are held in memory. Large-snapshot peak
  memory is not yet characterized. Native payloads retain upstream architecture
  and malformed-payload limitations; package framing is checked separately.
- Preprocessing manifests are compared semantically with only tools, timings
  and artifact paths excluded. Ordering, particle-sort and chunk metadata are
  compared across triplet/layout writers. Field metrics, error-bound summaries
  and batch aggregation are compared with Python on identical packages.
  Full compressed byte equality is not required: native codec settings and
  equivalent ordering may produce different streams; package sizes include JSON.
- Some incidental details differ: native APIs replace subprocess diagnostics;
  C++ rejects negative worker options during parsing even for preprocessing,
  where Python postpones worker validation. Work directories may contain fewer
  temporary adapter/interleaved files and internal timing subdivisions differ.
- No Python import facade is supplied. All four CLI commands and package families
  are implemented; experiments, historical scripts and large-data benchmarks
  are outside this translation. Exhaustive malformed native-stream fuzzing,
  thread sanitization and multi-platform validation were not performed.

## helpers.py compatibility inventory

The Python facade exports constants and aliases (`PYSZ_MIN_VALUES`,
`SZO_MIN_VALUES`), ToolPaths, JSON/raw I/O and worker/tool helpers, manifest size
and order queries, field/component ratios, metric accumulator operations,
comparison ordering, reporting, `parse_tool_stdout`, and
`manifest_path_from_manifest`. These belong in the C++ constants/models,
data/runtime, manifest and metrics components. This port does not provide a
Python import facade. CLI/package parity is the primary interface target.
The C++ API exposes pipeline, data, layout, codec and metric operations;
Python-only console/subprocess utility signatures are not exported.

## Validation

`./build.sh -DPARTICLE_REFERENCE_PYTHON=/home/moonlightplague/miniconda3/envs/h5/bin/python`
configures dependencies, builds and runs CTest. Reference runs use `python -B`,
PYTHONDONTWRITEBYTECODE=1, PYTHONPATH pointing to origin, and TMPDIR under build.

- Release: all six CTest suites pass, including 24 fieldwise, 32 triplet/blockwise
  and 20 layout encode/decode directions, plus native/runtime cases.
- Python reference: **96 tests passed** (2026-09-17).
- Batch reports are compared recursively against Python's builder using the same
  manifests, times and paths, including a failed child.
- Sanitizers: all six suites pass with ASan and UBSan on application code.
  Dependency libraries are not fully instrumented. Leak checking is disabled
  because this pass targets allocation safety, bounds and undefined behavior;
  allocation/deallocation mismatch checking remains enabled.
- Logs and generated fixtures stay in ignored `.deps/` and `build*/`.

## Small performance comparison

Reproduce with the reference Python environment:

```sh
python -B tests/benchmark.py --binary build/particle_compressors \
  --reference origin --work build/benchmark
```

20,000 synthetic particles, seed 343, float32 coordinates/velocities, uint64 IDs
above 2^53, three runs, relative bounds 0.001, one field worker, OMP_NUM_THREADS=1,
metrics enabled. Wall time includes process startup; memory is GNU time peak RSS.
These are local smoke measurements, not large-snapshot throughput claims.

| Codecs (position/velocity) | Runtime | Median wall (s) | Peak RSS (MiB) | Median package bytes |
| --- | --- | ---: | ---: | ---: |

| pcodec/pcodec | Python | 0.1551 | 50.54 | 422,985 |
| pcodec/pcodec | C++ | 0.0193 | 14.27 | 422,809 |
| szo/szo | Python | 0.1497 | 54.54 | 194,724 |
| szo/szo | C++ | 0.0199 | 16.23 | 194,549 |
| lcp/xnyzip | Python | 0.2390 | 50.61 | 161,249 |
| lcp/xnyzip | C++ | 0.0713 | 19.11 | 159,424 |

Individual runs and stage timings are saved in `build/benchmark/benchmark.json`.
All submodule revisions and the pre-existing dirty state in `origin/` are unchanged.
