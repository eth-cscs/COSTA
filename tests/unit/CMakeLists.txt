
set(unit_test_files
    # test_utils.cpp
    test.cpp
)

add_executable(unit ${unit_test_files})
target_include_directories(unit PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(unit PRIVATE gtest costa)
add_test(NAME unit COMMAND unit)
add_dependencies(tests unit)
