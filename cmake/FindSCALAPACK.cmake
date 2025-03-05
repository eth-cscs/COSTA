include(FindPackageHandleStandardArgs)

if(COSTA_SCALAPACK STREQUAL "MKL")
  find_package(MKL REQUIRED)
  get_target_property(COSTA_SCALAPACK_LINK_LIBRARIES costa::BLAS::MKL::scalapack_link
    INTERFACE_LINK_LIBRARIES)
elseif(COSTA_SCALAPACK STREQUAL "CRAY_LIBSCI")
  find_package(CRAY_LIBSCI REQUIRED)
  get_target_property(COSTA_SCALAPACK_LINK_LIBRARIES costa::BLAS::SCI::scalapack_link
    INTERFACE_LINK_LIBRARIES)
elseif(COSTA_SCALAPACK STREQUAL "NVPL")
  find_package(nvpl REQUIRED)
  get_target_property(COSTA_SCALAPACK_LINK_LIBRARIES costa::BLAS::NVPL::scalapack_link
    INTERFACE_LINK_LIBRARIES)
elseif(COSTA_SCALAPACK STREQUAL "CUSTOM")
  find_library(COSTA_SCALAPACK_LINK_LIBRARIES
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
endif()

find_package_handle_standard_args(SCALAPACK REQUIRED_VARS COSTA_SCALAPACK_LINK_LIBRARIES)

set(COSTA_SCALAPACK_FOUND "YES")

if (NOT TARGET costa::scalapack::scalapack)
  add_library(costa::scalapack::scalapack INTERFACE IMPORTED)
  set_target_properties(costa::scalapack::scalapack PROPERTIES INTERFACE_LINK_LIBRARIES
    "${COSTA_SCALAPACK_LINK_LIBRARIES}")
endif()

mark_as_advanced(COSTA_SCALAPACK_LINK_LIBRARIES COSTA_SCALAPACK_FOUND)
