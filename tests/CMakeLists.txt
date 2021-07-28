# add google test library
if (NOT TARGET gtest)
    add_library(gtest EXCLUDE_FROM_ALL STATIC gtest-all.cpp)
endif()
target_include_directories(gtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# All test targets, should be dependencies of tests target
add_custom_target(tests COMMENT "Builds all tests.")

# non-distributed, unit tests
add_subdirectory(unit)
