#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * A local block of the matrix.
 * data: a pointer to the start of the local matrix A_loc
 * ld: leading dimension or distance between two columns of A_loc
 * row: the global block row index
 * col: the global block colum index
 */
struct block_t {
    void *data;
    const int ld;
    const int row;
    const int col;
};

/**
 * Description of a distributed layout of a matrix
 * rowblocks: number of gobal blocks
 * colblocks: number of gobal blocks
 * rowsplit: [rowsplit[i], rowsplit[i+1]) is range of rows of block i
 * colsplit: [colsplit[i], colsplit[i+1]) is range of columns of block i
 */
struct grid_t {
    // global view
    int rowblocks;
    int colblocks;
    const int *rowsplit;
    const int *colsplit;
    const int *owners;
};

/**
 * Description of a distributed layout of a matrix
 * global_view: describes the global matrix grid
 * nlocalblock: number of blocks owned by the current rank
 * localblcoks: an array of block descriptions of the current rank
 */
struct layout_t {
    // global_view
    grid_t* grid;
    // local view (local blocks)
    int nlocalblocks;
    block_t* localblocks;
};

#ifdef __cplusplus
}
#endif
