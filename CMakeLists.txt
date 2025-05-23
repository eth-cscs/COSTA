cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

include(FetchContent)
include(CMakeDependentOption)
cmake_policy(SET CMP0048 NEW)

project(costa
  DESCRIPTION "Communication Optimal Shuffle and Transpose Algorithms"
  HOMEPAGE_URL "https://github.com/eth-cscs/COSTA"
  VERSION 2.2.2
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cmake/build_type.cmake)
include(cmake/adjust_mpiexec_flags.cmake)

# Options
set(COSTA_SCALAPACK_IMPLEMENTATION_LIST "OFF" "MKL" "CRAY_LIBSCI" "NVPL" "CUSTOM")

option(COSTA_WITH_TESTS "Generate the test target." OFF)
option(COSTA_WITH_APPS "Generate the miniapp targets." OFF)
option(COSTA_WITH_BENCHMARKS "Generate the benchmark targets." OFF)
option(COSTA_WITH_PROFILING "Enable profiling." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

set(COSTA_SCALAPACK "OFF" CACHE STRING "ScaLAPACK backend. Can be MKL, CRAY_LIBSCI, NVPL, CUSTOM or OFF.")
set_property(CACHE COSTA_SCALAPACK PROPERTY STRINGS ${COSTA_SCALAPACK_IMPLEMENTATION_LIST})

# always generate compile_commands.json for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS "YES")

# check if scalapack backend is valid
message(STATUS "Selected ScaLAPACK backend (or implementation) for COSTA: ${COSTA_SCALAPACK}")

if(NOT COSTA_SCALAPACK IN_LIST COSTA_SCALAPACK_IMPLEMENTATION_LIST)
  message(FATAL_ERROR "Invalid value for COSTA_SCALAPACK! It should be MKL, CRAY_LIBSCI, NVPL, CUSTOM or OFF.")
endif()

# preserve rpaths when installing and make the install folder relocatable
# use `CMAKE_SKIP_INSTALL_RPATH` to skip this
# https://spack.readthedocs.io/en/latest/workflows.html#write-the-cmake-build
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
# skip RPATH if the library is installed to system directories
if(isSystemDir STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    if(APPLE)
        set(basePoint @loader_path)
    else()
        set(basePoint $ORIGIN)
    endif()
    file(RELATIVE_PATH relDir ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir})
endif()

# Dependencies
#
find_package(MPI COMPONENTS CXX REQUIRED)
adjust_mpiexec_flags()

# OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -DDEBUG")

# Bundled dependencies
#
# Do not install cxxopts when including it.

if (COSTA_WITH_TESTS OR COSTA_WITH_APPS)
  FetchContent_Declare(cxxopts
                       GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
                       GIT_TAG        4bf61f08697b110d9e3991864650a405b3dd515d # v3.2.1
                       FIND_PACKAGE_ARGS NAMES cxxopts
  )
  FetchContent_MakeAvailable(cxxopts)
endif()

if (COSTA_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(COSTA_WITH_PROFILING)
  FetchContent_Declare(
    semiprof
    GIT_REPOSITORY  https://github.com/bcumming/semiprof.git
    GIT_TAG         f132142ff2215dfa073e416fa7911d8877d62752
    FIND_PACKAGE_ARGS NAMES semiprof
  )
  FetchContent_MakeAvailable(semiprof)
endif()

if (NOT COSTA_SCALAPACK STREQUAL "OFF")
  find_package(SCALAPACK REQUIRED)
endif()

# COSTA
#

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

add_subdirectory(src/costa)

install(DIRECTORY "${costa_SOURCE_DIR}/src/costa"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING
  PATTERN "*.hpp")

write_basic_package_version_file(
  "${costa_BINARY_DIR}/costaConfigVersion.cmake"
  VERSION ${costa_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_file("${costa_SOURCE_DIR}/cmake/costa.pc.in"
  "${costa_BINARY_DIR}/costa.pc"
  @ONLY)

configure_file("${costa_SOURCE_DIR}/cmake/costaConfig.cmake.in"
  "${costa_BINARY_DIR}/costaConfig.cmake"
  @ONLY)

install(FILES "${costa_BINARY_DIR}/costaConfig.cmake"
  "${costa_BINARY_DIR}/costaConfigVersion.cmake"
  "${costa_SOURCE_DIR}/cmake/FindMKL.cmake"
  "${costa_SOURCE_DIR}/cmake/FindCRAY_LIBSCI.cmake"
  "${costa_SOURCE_DIR}/cmake/Findnvpl.cmake"
  "${costa_SOURCE_DIR}/cmake/FindSCALAPACK.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/costa")

install(FILES "${costa_BINARY_DIR}/costa.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

if(COSTA_WITH_APPS)
  add_subdirectory(miniapps)
  add_subdirectory(examples)
endif()
