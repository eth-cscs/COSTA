#pragma once

// The header makes semiprof an optional dependency that needs not be shipped when costa is installed.
//
#ifdef GRID2GRID_WITH_PROFILING
#include <semiprof/semiprof.hpp>
#else
#define PE(name)
#define PL()
#define PP()
#define PC()
#endif
