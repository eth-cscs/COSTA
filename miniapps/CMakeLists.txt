################
#  Build test  #
################
# COSTA-only executables
set(executables "block_cyclic" "comm_volume")
foreach(exec ${executables})
    add_executable(${exec} "${exec}.cpp")
    target_link_libraries(${exec} PRIVATE costa cxxopts::cxxopts)
endforeach()

install(TARGETS ${executables} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

# COSTA without scalapack executable
set(scalapack_executables "pxtran_miniapp" "pxgemr2d_miniapp")

foreach(exec ${scalapack_executables})
    add_executable(${exec} "${exec}.cpp")
    target_link_libraries(${exec} PRIVATE costa_prefixed_scalapack cxxopts::cxxopts)
    if (COSTA_WITH_PROFILING)
        target_link_libraries(${exec} PRIVATE semiprof::semiprof)
        target_compile_definitions(${exec} PRIVATE COSTA_WITH_PROFILING)
    endif()
endforeach()
install(TARGETS ${scalapack_executables} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
