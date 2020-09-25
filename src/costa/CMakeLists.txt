set(INSTALLED_TARGETS_LIST "")
if(NOT TARGET costa)
    set(costa_src_files scalapack.cpp
                        costa_pxtran.cpp
    )

    add_library(costa ${costa_src_files})
    target_include_directories(costa PUBLIC $<BUILD_INTERFACE:${costa_SOURCE_DIR}/src>)
    target_compile_features(costa PUBLIC cxx_std_14)
    target_link_libraries(costa PUBLIC MPI::MPI_CXX
                                       grid2grid
                                       ${ScaLAPACK_TARGET}
    )

    list(APPEND INSTALLED_TARGETS_LIST "costa")
endif()

if(COSTA_WITH_INSTALL AND INSTALLED_TARGETS_LIST)
    install(TARGETS ${INSTALLED_TARGETS_LIST}
            EXPORT costa_targets
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

    install(EXPORT costa_targets
            FILE costaTargets.cmake
            NAMESPACE costa::
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/costa")
endif()

