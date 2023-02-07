cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

include(CMakeDependentOption)
cmake_policy(SET CMP0048 NEW)

project(costa
  DESCRIPTION "Communication Optimal Shuffle and Transpose Algorithms"
  HOMEPAGE_URL "https://github.com/eth-cscs/COSTA"
  VERSION 2.1.1
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cmake/build_type.cmake)
include(cmake/adjust_mpiexec_flags.cmake)

# Options
set(COSTA_SCALAPACK_IMPLEMENTATION_LIST "OFF" "MKL" "CRAY_LIBSCI" "CUSTOM")

option(COSTA_WITH_TESTS "Generate the test target." OFF)
option(COSTA_WITH_APPS "Generate the miniapp targets." OFF)
option(COSTA_WITH_BENCHMARKS "Generate the benchmark targets." OFF)
option(COSTA_WITH_PROFILING "Enable profiling." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

set(COSTA_SCALAPACK "OFF" CACHE STRING "ScaLAPACK backend. Can be MKL, CRAY_LIBSCI, CUSTOM or OFF.")
set_property(CACHE COSTA_SCALAPACK PROPERTY STRINGS ${COSTA_SCALAPACK_IMPLEMENTATION_LIST})

cmake_dependent_option(COSTA_CXXOPTS_BUILD_EXAMPLES "Enable examples build" OFF "COSTA_WITH_APPS" OFF)
cmake_dependent_option(COSTA_CXXOPTS_BUILD_TESTS "Enable tests suite" OFF "COSTA_WITH_TESTS" OFF)

# check if scalapack backend is valid
message(STATUS "Selected ScaLAPACK backend (or implementation) for COSTA: ${COSTA_SCALAPACK}")

if(NOT COSTA_SCALAPACK IN_LIST COSTA_SCALAPACK_IMPLEMENTATION_LIST)
  message(FATAL_ERROR "Invalid value for COSTA_SCALAPACK! It should be MKL, CRAY_LIBSCI, CUSTOM or OFF.")
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

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops -march=native -DNDEBUG")

# add_compile_options(-fauto-inc-dec -fbranch-count-reg -fcombine-stack-adjustments -fcompare-elim -fcprop-registers -fdce -fdefer-pop -fdelayed-branch -fdse -fforward-propagate -fguess-branch-probability -fif-conversion -fif-conversion2 -finline-functions-called-once -fipa-profile -fipa-pure-const -fipa-reference -fmerge-constants -fmove-loop-invariants -fomit-frame-pointer -freorder-blocks -fshrink-wrap -fshrink-wrap-separate -fsplit-wide-types -fssa-backprop -fssa-phiopt -ftree-bit-ccp -ftree-ccp -ftree-ch -ftree-coalesce-vars -ftree-copy-prop -ftree-dce -ftree-dominator-opts -ftree-dse -ftree-forwprop -ftree-fre -ftree-phiprop -ftree-pta -ftree-scev-cprop -ftree-sink -ftree-slsr -ftree-sra -ftree-ter -funit-at-a-time -falign-functions  -falign-jumps -falign-labels  -falign-loops -fcaller-saves -fcode-hoisting -fcrossjumping -fcse-follow-jumps  -fcse-skip-blocks -fdelete-null-pointer-checks -fdevirtualize  -fdevirtualize-speculatively -fexpensive-optimizations -fsplit-loops -fgcse  -fgcse-lm -fhoist-adjacent-loads -finline-functions -finline-small-functions -findirect-inlining -fipa-bit-cp  -fipa-cp  -fipa-icf -fipa-ra  -fipa-sra  -fipa-vrp -fisolate-erroneous-paths-dereference -flra-remat -foptimize-sibling-calls -foptimize-strlen -fpartial-inlining -fpeephole2 -freorder-blocks-algorithm=stc -freorder-blocks-and-partition  -freorder-functions -frerun-cse-after-loop -fschedule-insns  -fschedule-insns2 -fsched-interblock  -fsched-spec -fstore-merging -fstrict-aliasing -fthread-jumps -ftree-builtin-call-dce -ftree-pre -ftree-switch-conversion  -ftree-tail-merge -ftree-vrp -fgcse-after-reload -fipa-cp-clone -floop-interchange -floop-unroll-and-jam -fpeel-loops -fpredictive-commoning -fsplit-loops -fsplit-paths -ftree-loop-distribution -ftree-loop-vectorize -ftree-partial-pre -ftree-slp-vectorize -funswitch-loops -fvect-cost-model -fvect-cost-model=dynamic)


# Bundled dependencies
#

# Do not install cxxopts when including it.

if (COSTA_WITH_TESTS OR COSTA_WITH_APPS)
  find_package(cxxopts REQUIRED)
endif()

if (COSTA_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(COSTA_WITH_PROFILING)
  find_package(semiprof)
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
  "${costa_SOURCE_DIR}/cmake/FindSCALAPACK.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/costa")

install(FILES "${costa_BINARY_DIR}/costa.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

if(COSTA_WITH_APPS)
  add_subdirectory(miniapps)
  add_subdirectory(examples)
endif()
