#pragma once
#include <costa/grid2grid/grid_layout.hpp>
#include <costa/grid2grid/scalapack_layout.hpp>

namespace costa {
/**
 * A local block of the matrix.
 * data: a pointer to the start of the local matrix A_loc
 * ld: leading dimension or distance between two columns of A_loc
 * row: the global block row index
 * col: the global block colum index
 * ordering: 'R' (row-major) or 'C' (col-major).
 */
struct block_t {
    void *data;
    int ld;
    int row;
    int col;
};

/**
 * Description of a distributed layout of a matrix
 * rowblocks: number of gobal blocks
 * colblocks: number of gobal blocks
 * rowsplit: [rowsplit[i], rowsplit[i+1]) is range of rows of block i
 * colsplit: [colsplit[i], colsplit[i+1]) is range of columns of block i
 * owners: owners[i][j] is the rank owning block (i,j). 
 *         Owners are given in row-major order as assumed by C++.
 * nlocalblocks: number of blocks owned by the current rank
 * localblcoks: an array of block descriptions of the current rank
 * ordering: 'R' or 'C' depending on whether each local block is given in
 *           row- or col-major ordering
 */
template <typename T>
grid_layout<T> custom_layout(const int rowblocks,
                             const int colblocks,
                             const int* rowsplit,
                             const int* colsplit,
                             const int* owners,
                             const int nlocalblocks,
                             const block_t* localblocks,
                             const char ordering);

// contains only the global grid, without local data
assigned_grid2D custom_grid(const int rowblocks,
                            const int colblocks,
                            const int* rowsplit,
                            const int* colsplit,
                            const int* owners);

/**
 * creates a block cyclic layout (scalapack data layout) of some matrix A
 * and represents a submatrix sub(A) of the global matrix. The submatrix 
 * starts at (i, j) and has dimensions (sub_m, sub_n).
 * The detailed arguments are described below:
 * (m, n): global matrix dimensions of matrix A
 * (block_m, block_n): block dimensions
 * (i, j): submatrix start, 1-based. By default should be set to (1, 1).
 * (sub_m, sub_n): size of the submatrix sub(A). By default can be set to (m, n).
 * (proc_m, proc_n): processor, i.e. MPI ranks grid
 * (rsrc, crsrc): coordinates of ranks owning 
 *                the first row/col of the global matrix
 *                By default, should be set to (0, 0)
 * ptr: pointer to local data of the global matrix A.
 * lld: local leading dimension
 * ordering: 'R' or 'C' depending on whether each local block is given in 
 *           row- or col-major ordering
 * rank: MPI rank
 */
template <typename T>
grid_layout<T> block_cyclic_layout(
                   const int m, const int n, // global matrix dimensions
                   const int block_m, const int block_n, // block dimensions
                   const int i, const int j, // submatrix start
                                             // (1-based, scalapack-compatible)
                   const int sub_m, const int sub_n, // submatrix size
                   const int p_m, const int p_n, // processor grid dimension
                   const char order, // rank grid ordering ('R' or 'C')
                   const int rsrc, const int csrc, // coordinates of ranks oweing
                                                   // the first row (0-based)
                   T* ptr, // local data of matrix A
                           // (not the submatrix)
                   const int lld, // local leading dimension
                   const char ordering, // whether local blocks are row/col-major
                   const int rank // processor rank
               );

// same as block_cyclic_layout but without local data,
// so just the global grid
assigned_grid2D block_cyclic_grid(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int rsrc, const int csrc // coordinates of ranks oweing 
                                       // the first row 
);
}
