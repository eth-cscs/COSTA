#include <costa/layout.hpp>

template <typename T>
costa::grid_layout<T> costa::custom_layout(grid_t* grid,
                                     int n_local_blocks,
                                     block_t* local_blocks) {
    // Create the local blocks
    std::vector<costa::block<T>> loc_blks;

    // Create blocks
    for (int i = 0; i < n_local_blocks; ++i) {
        auto &block = local_blocks[i];
        auto row = block.row;
        auto col = block.col;
        auto ptr = reinterpret_cast<T *>(block.data);
        auto stride = block.ld;

        costa::block_coordinates coord{row, col};
        costa::interval rows{grid->rowsplit[row],
                                 grid->rowsplit[row + 1]};
        costa::interval cols{grid->colsplit[col],
                                 grid->colsplit[col + 1]};
        loc_blks.emplace_back(rows, cols, coord, ptr, stride);
    }

    // Grid specification
    std::vector<int> rows_split(grid->rowblocks + 1);
    std::copy_n(grid->rowsplit, rows_split.size(), rows_split.begin());

    std::vector<int> cols_split(grid->colblocks + 1);
    std::copy_n(grid->colsplit, cols_split.size(), cols_split.begin());

    int n_ranks = 1;
    std::vector<std::vector<int>> owners_matrix(grid->rowblocks);
    for (int i = 0; i < grid->rowblocks; ++i) {
        owners_matrix[i].resize(grid->colblocks);
        for (int j = 0; j < grid->colblocks; ++j) {
            owners_matrix[i][j] =
                grid->owners[j * grid->rowblocks + i];
            n_ranks = std::max(n_ranks, owners_matrix[i][j] + 1);
        }
    }

    return {{{std::move(rows_split), std::move(cols_split)},
             std::move(owners_matrix),
             n_ranks},
            {std::move(loc_blks)}};
}

template <typename T>
costa::grid_layout<T> block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int ia, const int ja, // coordinates of ranks oweing 
                                    // the first row 
                                    // (1-based, scalapack-compatible)
        const T* ptr, // local data of matrix A (not the submatrix)
        const int lld, // leading dimension
        const int rank // processor rank
        ) {

    auto scalapack_layout = costa::get_scalapack_grid<T>(
        lld,
        {m, n},
        {ia, ja},
        {sub_m, sub_n},
        {block_m, block_n},
        {proc_m, proc_n},
        rank_grid_ordering,
        {ia, ja},
        ptr,
        rank);

    return scalapack_layout;
}
