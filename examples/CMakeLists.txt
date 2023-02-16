################
#  Build test  #
################
# COSTA-only executables
add_executable(costa_example0 example0.cpp)
target_link_libraries(costa_example0 PRIVATE costa cxxopts::cxxopts)

add_executable(costa_example1 example1.cpp)
target_link_libraries(costa_example1 PRIVATE costa cxxopts::cxxopts)

add_executable(costa_copy_bench copy_bench.cpp)
target_link_libraries(costa_copy_bench PRIVATE costa cxxopts::cxxopts)

install(TARGETS costa_example0 costa_example1 costa_copy_bench RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
