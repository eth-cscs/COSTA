include(FindPackageHandleStandardArgs)

# we are using the GNU compiler
set(_sciname "sci_gnu_mpi_mp")
set(_sciname_acc "sci_acc_gnu_nv60")

find_library(COSTA_CRAY_LIBSCI_LIBRARIES
    NAMES ${_sciname_acc} ${_sciname}
    HINTS
    ${_SCALAPACK_LIBRARY_DIRS}
    ENV CRAY_LIBSCI_PREFIX_DIR
    ENV CRAY_LIBSCI_ACC_PREFIX_DIR
    PATH_SUFFIXES lib
    DOC "Path to the Cray-libsci library.")

message("CRAY_LIBSCI: ${CRAY_LIBSCI_LIBRARIES}")

find_package_handle_standard_args(CRAY_LIBSCI DEFAULT_MSG COSTA_CRAY_LIBSCI_LIBRARIES)

if (NOT TARGET costa::BLAS::SCI::scalapack)
  add_library(costa::BLAS::SCI::scalapack INTERFACE IMPORTED)
  set_target_properties(costa::BLAS::SCI::scalapack PROPERTIES INTERFACE_LINK_LIBRARIES "${COSTA_CRAY_LIBSCI_LIBRARIES}")
endif()
