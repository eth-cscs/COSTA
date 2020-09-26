################
#  Build test  #
################
# COSTA-only executables
set(executables "")
foreach(exec ${executables})
    add_executable(${exec} "${exec}.cpp")
    target_link_libraries(${exec} PRIVATE costa cxxopts)
    target_include_directories(${exec} PUBLIC "${cosma_SOURCE_DIR}/libs/cxxopts/include/")
endforeach()

# COSTA without scalapack executable
set(scalapack_executables "pxtran_miniapp")

foreach(exec ${scalapack_executables})
    add_executable(${exec} "${exec}.cpp")
    target_link_libraries(${exec} PRIVATE costa_prefixed_scalapack cxxopts)
    target_include_directories(${exec} PUBLIC "${cosma_SOURCE_DIR}/libs/cxxopts/include/")
endforeach()
