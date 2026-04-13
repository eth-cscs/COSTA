include(FetchContent)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG 52eb8108c5bdec04579160ae17225d66034bd723
  FIND_PACKAGE_ARGS NAMES GTest
)
FetchContent_MakeAvailable(googletest)

# All test targets, should be dependencies of tests target
add_custom_target(tests COMMENT "Builds all tests.")

# non-distributed, unit tests
add_subdirectory(unit)
