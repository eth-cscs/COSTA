set(INSTALLED_TARGETS_LIST "")
# original grid2grid files
set(costa_src_files 
  grid2grid/grid_layout.hpp
  grid2grid/block.cpp
  grid2grid/grid2D.cpp
  grid2grid/interval.cpp
  grid2grid/scalapack_layout.cpp
  grid2grid/communication_data.cpp
  grid2grid/grid_cover.cpp
  grid2grid/ranks_reordering.cpp
  grid2grid/utils.cpp
  grid2grid/transform.cpp
  grid2grid/transformer.hpp
  layout.cpp
)
set(costa_prefixed_scalapack_src_files 
  scalapack.cpp
  # templatized cpp implementations
  pxtran_op/costa_pxtran_op.cpp
  pxgemr2d/costa_pxgemr2d.cpp
  # scalapck api with costa prefix
  pxtran/prefixed_pxtran.cpp
  pxtranu/prefixed_pxtranu.cpp
  pxtranc/prefixed_pxtranc.cpp
  pxgemr2d/prefixed_pxgemr2d.cpp
)
set(costa_scalapack_src_files 
  scalapack.cpp
  # templatized cpp implementations
  pxtran_op/costa_pxtran_op.cpp
  pxgemr2d/costa_pxgemr2d.cpp
  # scalapck api overwritten
  pxtran/pxtran.cpp
  pxtranu/pxtranu.cpp
  pxtranc/pxtranc.cpp
  pxgemr2d/pxgemr2d.cpp
)

add_library(costa ${costa_src_files})
target_include_directories(costa PUBLIC 
  
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/grid2grid>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/pxgemr2d>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/pxtran_op>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/pxtran>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/pxtranu>
  $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src/pxtranc>
)
target_compile_features(costa PUBLIC cxx_std_14)
target_link_libraries(costa PUBLIC MPI::MPI_CXX
  OpenMP::OpenMP_CXX
  $<$<BOOL:${COSTA_WITH_PROFILING}>:semiprof::semiprof>)

# adds the -fPIC flag. This is automatically added by cmake
# for shared libraries, but for static libraries it has
# to be added manually. It's not needed for COSTA when built in isolation.
# However, if a 3rd-party library that is shared and built with the -fPIC
# flags is using COSTA built as a static library, it might be a problem 
# if COSTA is not compiled with this flag.
# if needed, can be enforced by -DCMAKE_POSITION_INDEPENDENT_CODE=ON

set_property(TARGET costa PROPERTY POSITION_INDEPENDENT_CODE ON)

if(COSTA_WITH_PROFILING)
  target_compile_definitions(costa PRIVATE
    $<BUILD_INTERFACE:COSTA_WITH_PROFILING>
  )
endif()

list(APPEND INSTALLED_TARGETS_LIST "costa")

if(NOT COSTA_SCALAPACK_IMPLEMENTATION STREQUAL "OFF")

  add_library(costa_scalapack ${costa_scalapack_src_files})
  add_library(costa_prefixed_scalapack ${costa_prefixed_scalapack_src_files})

  set_property(TARGET costa_prefixed_scalapack PROPERTY POSITION_INDEPENDENT_CODE ON)
  set_property(TARGET costa_scalapack PROPERTY POSITION_INDEPENDENT_CODE ON)

  target_link_libraries(costa_scalapack PUBLIC costa
    costa::scalapack
    $<$<BOOL:${COSTA_WITH_PROFILING}>:semiprof::semiprof>)

  target_link_libraries(costa_prefixed_scalapack PUBLIC costa 
    costa::scalapack
    $<$<BOOL:${COSTA_WITH_PROFILING}>:semiprof::semiprof>)

  if(COSTA_WITH_PROFILING)
    target_compile_definitions(costa_scalapack PRIVATE
      $<BUILD_INTERFACE:COSTA_WITH_PROFILING>
    )
    target_compile_definitions(costa_prefixed_scalapack PRIVATE
      $<BUILD_INTERFACE:COSTA_WITH_PROFILING>
    )
  endif()
  
  list(APPEND INSTALLED_TARGETS_LIST "costa_scalapack")
  list(APPEND INSTALLED_TARGETS_LIST "costa_prefixed_scalapack")
endif()

install(TARGETS ${INSTALLED_TARGETS_LIST}
  EXPORT costa_targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(EXPORT costa_targets
  FILE costaTargets.cmake
  NAMESPACE costa::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/costa")

