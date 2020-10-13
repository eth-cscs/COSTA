cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Check if COSTA is a subproject.
#
set(MASTER_PROJECT OFF)
if(NOT DEFINED PROJECT_NAME)
  set(MASTER_PROJECT ON)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cmake/build_type.cmake)
include(cmake/adjust_mpiexec_flags.cmake)

# Options
#
set(COSTA_SCALAPACK "OFF" CACHE STRING
    "Blas backend. Can be MKL, CRAY_LIBSCI, CUSTOM or OFF.")
set_property(CACHE COSTA_SCALAPACK PROPERTY STRINGS
    "OFF" "MKL" "CRAY_LIBSCI" "CUSTOM")

option(COSTA_WITH_TESTS "Generate the test target." ${MASTER_PROJECT})
option(COSTA_WITH_APPS "Generate the miniapp targets." ${MASTER_PROJECT})
option(COSTA_WITH_BENCHMARKS "Generate the benchmark targets." ${MASTER_PROJECT})
option(COSTA_WITH_INSTALL "Enable installation." ${MASTER_PROJECT})
option(COSTA_WITH_PROFILING "Enable profiling." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

# check if scalapack backend is valid
message(STATUS "Selected ScaLAPACK backend for COSTA: ${COSTA_SCALAPACK}")
unset(BACKEND_LIST)
get_property(BACKEND_LIST CACHE COSTA_SCALAPACK PROPERTY STRINGS)
if(COSTA_SCALAPACK AND NOT COSTA_SCALAPACK IN_LIST BACKEND_LIST)
  message(FATAL_ERROR "Invalid value for COSTA_SCALAPACK!")
endif()

project(costa VERSION 1.0.0 LANGUAGES CXX)

# preserve rpaths when installing and make the install folder relocatable
# use `CMAKE_SKIP_INSTALL_RPATH` to skip this
# https://spack.readthedocs.io/en/latest/workflows.html#write-the-cmake-build
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
          "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
# skip RPATH if SIRIUS is installed to system directories
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
find_package(MPI REQUIRED)
adjust_mpiexec_flags()

# Bundled dependencies
#

# Do not install cxxopts when including it.
option(CXXOPTS_BUILD_EXAMPLES OFF)
option(CXXOPTS_BUILD_TESTS OFF)
option(CXXOPTS_ENABLE_INSTALL OFF)
add_subdirectory(libs/cxxopts)

# ScaLAPACK providers
#
set(SCALAPACK_TARGET "")
if (${COSTA_SCALAPACK} STREQUAL "MKL")
    find_package(MKL REQUIRED COMPONENTS SCALAPACK_32BIT_OMP)
    set(SCALAPACK_TARGET "mkl::scalapack_32bit_omp")
elseif (${COSTA_SCALAPACK} STREQUAL "CRAY_LIBSCI")
    find_package(CRAY_LIBSCI REQUIRED)
    set(SCALAPACK_TARGET "${CRAY_LIBSCI_LIBRARIES}")
elseif (${COSTA_SCALAPACK} STREQUAL "CUSTOM")
    find_package(SCALAPACK REQUIRED)
    set(SCALAPACK_TARGET "${SCALAPACK_LIBRARIES}")
else  ()
    message(WARNING "No ScaLAPACK interface is specified. Building COSTA without SCALAPACK wrappers.")
endif ()

# COSTA
#
add_subdirectory(src/costa)

if(COSTA_WITH_PROFILING)
    option(SEMIPROF_WITH_INSTALL "" ${COSTA_WITH_INSTALL})
    add_subdirectory(libs/semiprof)
endif()

if(COSTA_WITH_INSTALL)
    include(CMakePackageConfigHelpers)
    include(GNUInstallDirs)

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
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/costa")

    install(FILES "${costa_BINARY_DIR}/costa.pc"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endif()

if(COSTA_WITH_APPS)
    add_subdirectory(miniapps)
endif()
