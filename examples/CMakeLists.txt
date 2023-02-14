################
#  Build test  #
################
# COSTA-only executables
set(executables "example0" "example1" "copy_bench")
foreach(exec ${executables})
    add_executable(${exec} "${exec}.cpp")
    target_link_libraries(${exec} PRIVATE costa cxxopts::cxxopts)
endforeach()

# COSTA without scalapack executable
set(scalapack_executables "")

if (SCALAPACK_TARGET) 
    foreach(exec ${scalapack_executables})
        add_executable(${exec} "${exec}.cpp")
        target_link_libraries(${exec} PRIVATE costa_scalapack cxxopts::cxxopts)
        if (COSTA_WITH_PROFILING)
            target_link_libraries(${exec} PRIVATE semiprof::semiprof)
            target_compile_definitions(${exec} PRIVATE COSTA_WITH_PROFILING)
        endif()
    endforeach()
endif()
