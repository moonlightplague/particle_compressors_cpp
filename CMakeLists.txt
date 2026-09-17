cmake_minimum_required(VERSION 3.20)
project(
  particle_compressors
  VERSION 0.1.0
  LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(CTest)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
# Discover system installations or a pre-existing Conan package cache.
file(GLOB PARTICLE_CONAN_PACKAGES "$ENV{HOME}/.conan2/p/*/p"
     "$ENV{HOME}/.conan2/p/b/*/p")
find_path(ZSTD_INCLUDE_DIR zstd.h
          HINTS ${PARTICLE_CONAN_PACKAGES} PATH_SUFFIXES include)
find_library(ZSTD_LIBRARY zstd
             HINTS ${PARTICLE_CONAN_PACKAGES} PATH_SUFFIXES lib)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
  message(FATAL_ERROR
    "Zstandard development files were not found. Install libzstd-dev on "
    "Debian/Ubuntu, or set ZSTD_INCLUDE_DIR and ZSTD_LIBRARY to an existing "
    "installation.")
endif()
add_library(particle_zstd UNKNOWN IMPORTED)
set_target_properties(particle_zstd PROPERTIES
  IMPORTED_LOCATION "${ZSTD_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${ZSTD_INCLUDE_DIR}")
find_program(CARGO_EXECUTABLE cargo HINTS "${CMAKE_SOURCE_DIR}/.deps/cargo/bin"
                                          REQUIRED)
set(YAML_CPP_BUILD_TESTS
    OFF
    CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS
    OFF
    CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL
    OFF
    CACHE BOOL "" FORCE)
add_subdirectory(tools/yaml-cpp EXCLUDE_FROM_ALL)
set(PCO_TARGET_DIR "${CMAKE_BINARY_DIR}/cargo")
add_custom_command(
  OUTPUT "${PCO_TARGET_DIR}/release/libcpcodec.a"
  COMMAND
    ${CMAKE_COMMAND} -E env CARGO_TARGET_DIR=${PCO_TARGET_DIR}
    ${CARGO_EXECUTABLE} build --manifest-path
    ${CMAKE_SOURCE_DIR}/tools/pcodec/Cargo.toml --release --locked -p cpcodec
  DEPENDS tools/pcodec/pco_c/src/lib.rs tools/pcodec/Cargo.lock
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools/pcodec
  VERBATIM)
add_custom_target(pcodec_build DEPENDS "${PCO_TARGET_DIR}/release/libcpcodec.a")
add_library(pcodec STATIC IMPORTED GLOBAL)
set_target_properties(
  pcodec PROPERTIES IMPORTED_LOCATION "${PCO_TARGET_DIR}/release/libcpcodec.a")
add_dependencies(pcodec pcodec_build)
# Both upstream C APIs export SZ_compress_args/SZ_decompress. Rename them and
# their template entry points in separate translation units to avoid collisions.
set(SZo_VERSION 1.0.0)
set(SZo_DATA_VERSION 1.0.0)
include(FetchContent)
# pysz 1.0.3 uses the SZ3 3.3.2 stream format, unlike LCP's SZ3 3.1.8.
FetchContent_Declare(
  sz3_compat
  URL https://codeload.github.com/szcompressor/SZ3/tar.gz/0ebe6fa26b0aa2e274d13b78df81a81e2db1b268
      DOWNLOAD_EXTRACT_TIMESTAMP
      TRUE
  URL_HASH
    SHA256=5c225162df7b8e65f79661eae17efa32a44164c4fa1479e6d09f5348f3a3d1b4)
FetchContent_GetProperties(sz3_compat)
if(NOT sz3_compat_POPULATED)
  FetchContent_Populate(sz3_compat)
endif()
set(SZ3_DATA_VERSION 3.3.2)
function(particle_codec_version name major minor patch template)
  set(PROJECT_NAME ${name})
  set(PROJECT_VERSION "${major}.${minor}.${patch}")
  set(PROJECT_VERSION_MAJOR ${major})
  set(PROJECT_VERSION_MINOR ${minor})
  set(PROJECT_VERSION_PATCH ${patch})
  set(PROJECT_VERSION_TWEAK 0)
  configure_file(${template} generated/${name}/version.hpp)
endfunction()
particle_codec_version(SZo 1 0 0 tools/SZo/include/SZo/version.hpp.in)
particle_codec_version(SZ3 3 3 2
                       ${sz3_compat_SOURCE_DIR}/include/SZ3/version.hpp.in)
file(GLOB SZO_FSE CONFIGURE_DEPENDS tools/SZo/include/SZo/encoder/zfse/*.c)
add_library(particle_szo STATIC src/szo_profile.cpp
                                tools/SZo/tools/szoc/src/szo.cpp ${SZO_FSE})
target_include_directories(
  particle_szo PRIVATE tools/SZo/include tools/SZo/tools/szoc/include
                       "${CMAKE_BINARY_DIR}/generated")
target_compile_definitions(
  particle_szo
  PRIVATE SZ_compress_args=particle_szo_compress
          SZ_decompress=particle_szo_decompress
          SZ_compress=particle_szo_template_compress
          SZo_DEBUG_TIMINGS=0
          SZo_USE_SKA_HASH=1)
# SZo bundles entropy functions also exported by static libzstd. Keep its
# implementation private so the two codec versions cannot collide at link time.
foreach(symbol
    ERR_getErrorString
    FSE_NCountWriteBound
    FSE_buildCTable
    FSE_buildCTable_raw
    FSE_buildCTable_rle
    FSE_buildCTable_wksp
    FSE_buildDTable
    FSE_buildDTable_raw
    FSE_buildDTable_rle
    FSE_compress
    FSE_compress2
    FSE_compressBound
    FSE_compress_usingCTable
    FSE_compress_wksp
    FSE_createCTable
    FSE_createDTable
    FSE_decompress
    FSE_decompress_usingDTable
    FSE_decompress_wksp
    FSE_freeCTable
    FSE_freeDTable
    FSE_getErrorName
    FSE_isError
    FSE_normalizeCount
    FSE_optimalTableLog
    FSE_optimalTableLog_internal
    FSE_readNCount
    FSE_versionNumber
    FSE_writeNCount
    HIST_count
    HIST_countFast
    HIST_countFast_wksp
    HIST_count_simple
    HIST_count_wksp
    HIST_isError
    HUF_getErrorName
    HUF_isError
    HUF_readStats
)
  target_compile_definitions(particle_szo PRIVATE "${symbol}=particle_szo_${symbol}")
endforeach()
target_link_libraries(particle_szo PRIVATE particle_zstd OpenMP::OpenMP_CXX)
add_library(particle_sz3 STATIC tools/LCP/tools/sz3c/src/sz3c.cpp)
target_include_directories(
  particle_sz3
  PRIVATE ${sz3_compat_SOURCE_DIR}/include tools/LCP/tools/sz3c/include
          "${CMAKE_BINARY_DIR}/generated")
target_compile_definitions(
  particle_sz3
  PRIVATE SZ3=ParticleSZ3 SZ_compress_args=particle_sz3_compress
          SZ_decompress=particle_sz3_decompress
          SZ_compress=particle_sz3_template_compress SZ3_DEBUG_TIMINGS=0)
target_link_libraries(particle_sz3 PRIVATE particle_zstd OpenMP::OpenMP_CXX)
add_library(
  particle_pipeline
  src/data.cpp
  src/json_io.cpp
  src/huffman.cpp
  src/layout.cpp
  src/layout_pipeline.cpp
  src/shaped.cpp
  src/native.cpp
  src/metrics.cpp
  src/batch.cpp
  src/hdf5.cpp
  src/codecs.cpp
  src/pipeline.cpp
  src/cli.cpp)
target_include_directories(
  particle_pipeline
  PUBLIC include ${HDF5_INCLUDE_DIRS}
  PRIVATE tools/pcodec/pco_c/include)
target_link_libraries(
  particle_pipeline
  PUBLIC ${HDF5_CXX_LIBRARIES} yaml-cpp::yaml-cpp
  PRIVATE pcodec particle_sz3 particle_szo Threads::Threads ${CMAKE_DL_LIBS})
target_compile_definitions(
  particle_pipeline
  PRIVATE PARTICLE_DEFAULT_CONFIG="${CMAKE_SOURCE_DIR}/config.yaml")
target_compile_options(particle_pipeline PRIVATE -Wall -Wextra -Wpedantic)
add_executable(particle_compressors src/main.cpp)
target_link_libraries(particle_compressors PRIVATE particle_pipeline)
if(BUILD_TESTING)
  add_executable(particle_core_test tests/core.cpp)
  target_link_libraries(particle_core_test PRIVATE particle_pipeline)
  add_test(NAME core COMMAND particle_core_test)
endif()
include(cmake/NativeTriplets.cmake)
set(PARTICLE_REFERENCE_PYTHON
    ""
    CACHE FILEPATH
          "Python with h5py/numpy/pcodec/pysz/pyszo for differential tests")
if(BUILD_TESTING AND PARTICLE_REFERENCE_PYTHON)
  foreach(suite interop triplet_interop layout_interop runtime_interop
                contracts_interop)
    add_test(
      NAME ${suite}
      COMMAND
        ${PARTICLE_REFERENCE_PYTHON} -B ${CMAKE_SOURCE_DIR}/tests/${suite}.py
        --binary $<TARGET_FILE:particle_compressors> --reference
        ${CMAKE_SOURCE_DIR}/origin --work ${CMAKE_BINARY_DIR}/tests/${suite})
    set_tests_properties(
      ${suite} PROPERTIES TIMEOUT 1200 ENVIRONMENT
                          "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1")
  endforeach()
endif()
option(PARTICLE_SANITIZERS
       "Instrument application code with address/undefined sanitizers" OFF)
if(PARTICLE_SANITIZERS)
  target_compile_options(particle_pipeline PRIVATE -fsanitize=address,undefined
                                                   -fno-omit-frame-pointer)
  target_link_options(particle_pipeline PUBLIC -fsanitize=address,undefined)
endif()
