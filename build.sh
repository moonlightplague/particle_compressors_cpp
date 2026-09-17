#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -x "$PWD/.deps/cargo/bin/cargo" ]]; then
  export CARGO_HOME="$PWD/.deps/cargo" RUSTUP_HOME="$PWD/.deps/rustup"
  export PATH="$CARGO_HOME/bin:$PATH"
fi
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT="${HDF5_ROOT:-$HOME/.local/hdf5}" "$@"
cmake --build build --parallel "${BUILD_JOBS:-2}"
ctest --test-dir build --output-on-failure
