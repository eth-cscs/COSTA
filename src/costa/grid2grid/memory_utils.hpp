#pragma once
#include <costa/grid2grid/block.hpp>
#include <costa/grid2grid/tiling_manager.hpp>

#include <algorithm>
#include <cstring>
#include <complex>
#include <cmath>
#include <type_traits>
#include <utility>
#include <omp.h>
#include <memory>

namespace costa {
namespace memory {

// copies n entries of elem_type from src_ptr to desc_ptr
// if alpha!=1 or beta != 0, then also performs: 
//     dest[i] = alpha * src[i] + beta*dest[i]
template <typename elem_type>
void copy(const std::size_t n, const elem_type *src_ptr,
          elem_type *dest_ptr,
          const bool should_conjugate = false,
          const elem_type alpha=elem_type{1},
          const elem_type beta=elem_type{0}) {
    static_assert(std::is_trivially_copyable<elem_type>(),
                  "Element type must be trivially copyable!");
    bool perform_operation = std::abs(alpha - elem_type{1}) > 0 || std::abs(beta - elem_type{0}) > 0;
    if (!perform_operation && !should_conjugate) {
        std::memcpy(dest_ptr, src_ptr, sizeof(elem_type) * n);
        assert(dest_ptr[0] == src_ptr[0]);
    } else {
        for (int i = 0; i < n; ++i) {
            auto el = src_ptr[i];
            if (should_conjugate) {
                el = conjugate_f(el);
            }
            dest_ptr[i] = beta * dest_ptr[i] + alpha * el;
        }
    }
}

// copies 2D block of given size from src_ptr with stride ld_src
// to dest_ptr with stride ld_dest
template <class elem_type>
void copy2D(int n_rows, int n_cols,
            const elem_type *src_ptr, const int ld_src,
            elem_type *dest_ptr, const int ld_dest,
            const bool should_conjugate = false,
            const elem_type alpha = elem_type{1},
            const elem_type beta = elem_type{0},
            const bool col_major = true) {
    static_assert(std::is_trivially_copyable<elem_type>(),
                  "Element type must be trivially copyable!");
    auto block_size = n_rows * n_cols;
    // std::cout << "invoking copy2D." << std::endl;
    assert(block_size >= 0);

    // stop if 0-sized
    if (block_size == 0) return;

    if (!col_major) {
        std::swap(n_rows, n_cols);
    }

    // if not strided, copy in a single piece
    if (n_rows == (size_t)ld_src &&
        n_rows == (size_t)ld_dest) {
        copy(block_size, src_ptr, dest_ptr, should_conjugate, alpha, beta);
    } else {
        // if strided, copy column-by-column
        // #pragma omp task firstprivate(dim, src_ptr, ld_src, dest_ptr, ld_dest)
        for (size_t col = 0; col < n_cols; ++col) {
            // #pragma omp task firstprivate(dim, src_ptr, ld_src, col, dest_ptr, ld_dest)
            copy(n_rows,
                 src_ptr + ld_src * col,
                 dest_ptr + ld_dest * col, 
                 should_conjugate, 
                 alpha, beta);
        }
    }
}

// transpose (out of place) data that is in col-major order
template <typename T>
void transpose_col_major(const int n_rows, const int n_cols, 
               const T* src_ptr, const int src_stride, 
               T* dest_ptr, const int dest_stride, 
               const bool should_conjugate, 
               const T alpha, const T beta,
               const tiling_manager<T>& tiling) {
    static_assert(std::is_trivially_copyable<T>(),
            "Element type must be trivially copyable!");
    // n_rows and n_cols before transposing
    // int block_dim = std::max(8, 128/(int)sizeof(T));
    int block_dim = tiling.block_dim;

    int n_blocks_row = (n_rows+block_dim-1)/block_dim;
    int n_blocks_col = (n_cols+block_dim-1)/block_dim;
    int n_blocks = n_blocks_row * n_blocks_col;

    int n_threads = std::min(n_blocks, tiling.max_threads);

    bool perform_operation = std::abs(alpha - T{1}) > 0 || std::abs(beta - T{0}) > 0;

#pragma omp parallel for num_threads(n_threads)
    for (int block = 0; block < n_blocks; ++block) {
        int thread_id = omp_get_thread_num();
        int b_offset = thread_id * block_dim;

        // col-major traversing blocks
        int block_i = (block % n_blocks_row) * block_dim;
        int block_j = (block / n_blocks_row) * block_dim;

        int upper_i = std::min(n_rows, block_i + block_dim);
        int upper_j = std::min(n_cols, block_j + block_dim);

        if (block_i == block_j) {
            for (int i = block_i; i < upper_i; ++i) {
                for (int j = block_j; j < upper_j; ++j) {
                    // (i, j) in the original block, column-major
                    auto el = src_ptr[j * src_stride + i];
                    // auto el = b.local_element(i, j);
                    // (j, i) in the send buffer, column-major
                    if (should_conjugate)
                        el = conjugate_f(el);
                    tiling.buffer[b_offset + j-block_j] = el;
                }
                for (int j = block_j; j < upper_j; ++j) {
                    auto& dst = dest_ptr[i*dest_stride + j];
                    if (perform_operation) {
                        dst = beta * dst + alpha * tiling.buffer[b_offset + j-block_j];
                    } else {
                        dst = tiling.buffer[b_offset + j-block_j];
                    }
                }
            }
        } else {
            // #pragma omp task firstprivate(block_i, block_j, dest_ptr, ptr, stride, conj, n_rows_t)
            for (int i = block_i; i < upper_i; ++i) {
                for (int j = block_j; j < upper_j; ++j) {
                    auto el = src_ptr[j * src_stride + i];
                    // (i, j) in the original block, column-major
                    // auto el = b.local_element(i, j);
                    // (j, i) in the send buffer, column-major
                    if (should_conjugate)
                        el = conjugate_f(el);
                    auto& dst = dest_ptr[i*dest_stride + j];
                    if (perform_operation) {
                        dst = beta * dst + alpha * el;
                    } else {
                        dst = el;
                    }
                }
            }
        }
    }
}

// transpose (out of place) data that is in row-major order
template <typename T>
void transpose_row_major(const int n_rows, const int n_cols, 
               const T* src_ptr, const int src_stride, 
               T* dest_ptr, const int dest_stride, 
               const bool should_conjugate, 
               const T alpha, const T beta,
               const tiling_manager<T>& tiling) {
    static_assert(std::is_trivially_copyable<T>(),
            "Element type must be trivially copyable!");
    // n_rows and n_cols before transposing
    // int block_dim = std::max(8, 128/(int)sizeof(T));
    int block_dim = tiling.block_dim;

    int n_blocks_row = (n_rows+block_dim-1)/block_dim;
    int n_blocks_col = (n_cols+block_dim-1)/block_dim;
    int n_blocks = n_blocks_row * n_blocks_col;

    int n_threads = std::min(n_blocks, tiling.max_threads);

    bool perform_operation = std::abs(alpha - T{1}) > 0 || std::abs(beta - T{0}) > 0;

#pragma omp parallel for num_threads(n_threads)
    for (int block = 0; block < n_blocks; ++block) {
        int thread_id = omp_get_thread_num();
        int b_offset = thread_id * block_dim;

        // row-major traversing blocks
        int block_i = (block / n_blocks_row) * block_dim;
        int block_j = (block % n_blocks_row) * block_dim;

        int upper_i = std::min(n_rows, block_i + block_dim);
        int upper_j = std::min(n_cols, block_j + block_dim);

        if (block_i == block_j) {
            for (int j = block_j; j < upper_j; ++j) {
                for (int i = block_i; i < upper_i; ++i) {
                    // (i, j) in the original block, column-major
                    auto el = src_ptr[i * src_stride + j];
                    // auto el = b.local_element(i, j);
                    // (j, i) in the send buffer, column-major
                    if (should_conjugate)
                        el = conjugate_f(el);
                    tiling.buffer[b_offset + i-block_i] = el;
                }
                for (int i = block_i; i < upper_i; ++i) {
                    auto& dst = dest_ptr[j*dest_stride + i];
                    if (perform_operation) {
                        dst = beta * dst + alpha * tiling.buffer[b_offset + i-block_i];
                    } else {
                        dst = tiling.buffer[b_offset + i-block_i];
                    }
                }
            }
        } else {
            // #pragma omp task firstprivate(block_i, block_j, dest_ptr, ptr, stride, conj, n_rows_t)
            for (int j = block_j; j < upper_j; ++j) {
                for (int i = block_i; i < upper_i; ++i) {
                    auto el = src_ptr[i * src_stride + j];
                    // (i, j) in the original block, column-major
                    // auto el = b.local_element(i, j);
                    // (j, i) in the send buffer, column-major
                    if (should_conjugate)
                        el = conjugate_f(el);
                    auto& dst = dest_ptr[j*dest_stride + i];
                    if (perform_operation) {
                        dst = beta * dst + alpha * el;
                    } else {
                        dst = el;
                    }
                }
            }
        }
    }
}

template <typename T>
void transpose(const int n_rows, const int n_cols, 
               const T* src_ptr, const int src_stride, 
               T* dest_ptr, const int dest_stride, 
               const bool should_conjugate, 
               const T alpha, const T beta,
               const bool col_major,
               const tiling_manager<T>& tiling) {
    if (col_major) {
        transpose_col_major(n_rows, n_cols, 
                            src_ptr, src_stride,
                            dest_ptr, dest_stride,
                            should_conjugate,
                            alpha, beta,
                            tiling);
    } else {
        transpose_row_major(n_rows, n_cols,
                            src_ptr, src_stride,
                            dest_ptr, dest_stride,
                            should_conjugate,
                            alpha, beta,
                            tiling);
    }
}

inline
int default_stride(int n_rows, int n_cols, bool should_transpose, bool col_major) {
    int stride = 0;
    if (should_transpose) {
        stride = col_major ? n_rows : n_cols;
    } else {
        stride = col_major ? n_cols : n_rows;
    }
    return stride;
}


template <typename T>
void copy_and_transform(const int n_rows, const int n_cols,
                        const T* src_ptr, int src_stride,
                        const bool src_col_major,
                        T* dest_ptr, int dest_stride,
                        const bool dest_col_major,
                        const bool should_transpose,
                        const bool should_conjugate,
                        const T alpha, const T beta,
                        const tiling_manager<T>& tiling) {
    // BE CAREFUL: transpose and different src and dest orderings might cancel out
    // ===========
    // Row-major + Transpose + Row-major = Transpose (Row-major)
    // Col-major + Transpose + Col-major = Transpose (Col-major)
    // Row-major + Transpose + Col-major = Copy(Row-major) // cancels out
    // Col-major + Transpose + Row-major = Copy(Col-major) // cancels out
    //
    // Row-major + NoTranspose + Row-major = Copy(Row-major)
    // Col-major + NoTranspose + Col-major = Copy(Col-major)
    // Row-major + NoTranspose + Col-major = Transpose(Row-major)
    // Col-major + NoTranspose + Row-major = Transpose(Col-major)
    bool will_transpose = (should_transpose && src_col_major == dest_col_major)
                           ||
                          (!should_transpose && src_col_major != dest_col_major);

    assert(dest_stride >= 0);

    // if dest_stride == 0, then no stride (i.e. default stride)
    if (dest_stride == 0) {
        dest_stride = default_stride(n_rows, n_cols,
                                     will_transpose, dest_col_major);
    }
    // if src_stride == 0, then no stride (i.e. default stride)
    if (src_stride == 0) {
        // src is not transposed
        src_stride = default_stride(n_rows, n_cols,
                                    false, src_col_major);
    }

    if (will_transpose) {
        // transpose dimensions
        transpose(n_rows, n_cols,
                  src_ptr, src_stride, 
                  dest_ptr, dest_stride, 
                  should_conjugate, 
                  alpha, beta,
                  src_col_major,
                  tiling);
    } else {
        copy2D(n_rows, n_cols,
               src_ptr, src_stride,
               dest_ptr, dest_stride,
               should_conjugate,
               alpha, beta,
               src_col_major);
    }
}

} // namespace memory
} // namespace costa
