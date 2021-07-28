#pragma once
#include <memory>
#include <omp.h>

namespace costa {
namespace memory {

template <typename T>
struct threads_workspace {
    threads_workspace() = default;

    threads_workspace(int block_dim, int max_threads=omp_get_max_threads())
        : block_dim(block_dim)
        , max_threads(max_threads) {
        // std::cout << "max_threads = " << max_threads << std::endl;
        buffer = std::unique_ptr<T[]>(new T[block_dim * max_threads]);
    }

    int block_dim = 0;
    int max_threads = 0;
    std::unique_ptr<T[]> buffer;
};
}}
