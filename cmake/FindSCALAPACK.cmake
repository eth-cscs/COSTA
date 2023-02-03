include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)

find_package(MKL)
find_package(CRAY_LIBSCI)

if (TARGET costa::BLAS::MKL::scalapack_link OR TARGET costa::BLAS::SCI::scalapack_link)
  if (TARGET costa::BLAS::MKL::scalapack_link AND COSTA_SCALAPACK_IMPLEMENTATION MATCHES "MKL")
    get_target_property(COSTA_SCALAPACK_LINK_LIBRARIES costa::BLAS::MKL::scalapack_link
      INTERFACE_LINK_LIBRARIES)
  endif()
  
  if (TARGET costa::BLAS::SCI::scalapack_link  AND COSTA_SCALAPACK_IMPLEMENTATION MATCHES "CRAY_LIBSCI")
    get_target_property(COSTA_SCALAPACK_LINK_LIBRARIES costa::BLAS::SCI::scalapack_link
      INTERFACE_LINK_LIBRARIES)
  endif ()
else ()
  pkg_search_module(COSTA_SCALAPACK scalapack)
  
  if (NOT COSTA_SCALAPACK_FOUND)
    find_library(COSTA_SCALAPACK
      NAMES scalapack 
      HINTS
      ENV SCALAPACKROOT
      ENV SCALAPACK_ROOT
      ENV ORNL_SCALAPACK_ROOT
      ENV SCALAPACK_PREFIX
      ENV SCALAPACK_DIR
      ENV SCALAPACKDIR
      ENV /usr/bin
      PATH_SUFFIXES lib
      DOC "Path to the scalapack library.")
  endif()
endif()
message(INFO "COSTA DEBUG : ${COSTA_SCALAPACK_LINK_LIBRARIES}")

find_package_handle_standard_args(SCALAPACK REQUIRED_VARS COSTA_SCALAPACK_LINK_LIBRARIES)

set(COSTA_SCALAPACK_FOUND "YES")

if (NOT TARGET costa::scalapack)
  add_library(costa::scalapack INTERFACE IMPORTED)
endif()

set_target_properties(
  costa::scalapack PROPERTIES INTERFACE_LINK_LIBRARIES
  "${COSTA_SCALAPACK_LINK_LIBRARIES}")

mark_as_advanced(COSTA_SCALAPACK_LINK_LIBRARIES COSTA_SCALAPACK_FOUND)
