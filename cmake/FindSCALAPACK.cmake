include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)

find_package(MKL)
find_package(CRAY_LIBSCI)

if(COSTA_SCALAPACK_IMPLEMENTATION STREQUAL "MKL")
  get_target_property(COSTA_SCALAPACK costa::BLAS::MKL::scalapack_link
    INTERFACE_LINK_LIBRARIES)
elseif(COSTA_SCALAPACK_IMPLEMENTATION STREQUAL "CRAY_LIBSCI")
  get_target_property(COSTA_SCALAPACK costa::BLAS::SCI::scalapack_link
    INTERFACE_LINK_LIBRARIES)
elseif(COSTA_SCALAPACK_IMPLEMENTATION STREQUAL "CUSTOM")
  pkg_search_module(_COSTA_SCALAPACK scalapack)

  find_library(COSTA_SCALAPACK
    NAMES scalapack
    HINTS
    ${_COSTA_SCALAPACK_LIBRARY_DIRS}
    ENV SCALAPACKROOT
    ENV SCALAPACK_ROOT
    ENV ORNL_SCALAPACK_ROOT
    ENV SCALAPACK_PREFIX
    ENV SCALAPACK_DIR
    ENV SCALAPACKDIR
    /usr/bin
    PATH_SUFFIXES lib
    DOC "Path to the scalapack library.")

elseif()
  message(ERROR "Unknown COSTA_SCALAPACK_IMPLEMENTATION: ${COSTA_SCALAPACK_IMPLEMENTATION}")
endif()

message(INFO "COSTA DEBUG : ${COSTA_SCALAPACK}")

find_package_handle_standard_args(SCALAPACK REQUIRED_VARS COSTA_SCALAPACK)

set(COSTA_SCALAPACK_FOUND "YES")

if (NOT TARGET costa::scalapack)
  add_library(costa::scalapack INTERFACE IMPORTED)
endif()

set_target_properties(
  costa::scalapack PROPERTIES INTERFACE_LINK_LIBRARIES
  "${COSTA_SCALAPACK_LINK_LIBRARIES}")

mark_as_advanced(COSTA_SCALAPACK_LINK_LIBRARIES COSTA_SCALAPACK_FOUND)
