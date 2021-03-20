#include <costa/layout.hpp>
#include <complex>

costa::assigned_grid2D costa::custom_grid(const int rowblocks,
                                          const int colblocks,
                                          const int* rowsplit,
                                          const int* colsplit,
                                          const int* owners) {
    // Grid specification
    std::vector<int> rows_split(rowblocks + 1);
    std::copy_n(rowsplit, rows_split.size(), rows_split.begin());

    std::vector<int> cols_split(colblocks + 1);
    std::copy_n(colsplit, cols_split.size(), cols_split.begin());

    int n_ranks = 1;
    std::vector<std::vector<int>> owners_matrix(rowblocks);
    for (int i = 0; i < rowblocks; ++i) {
        owners_matrix[i].resize(colblocks);
        for (int j = 0; j < colblocks; ++j) {
            owners_matrix[i][j] =
                owners[i * colblocks + j];
            n_ranks = std::max(n_ranks, owners_matrix[i][j] + 1);
        }
    }

    return assigned_grid2D{{std::move(rows_split), std::move(cols_split)},
                          std::move(owners_matrix),
                          n_ranks};
}

template <typename T>
costa::grid_layout<T> costa::custom_layout(const int rowblocks,
                                           const int colblocks,
                                           const int* rowsplit,
                                           const int* colsplit,
                                           const int* owners,
                                           const int nlocalblocks,
                                           const block_t* localblocks,
                                           const char ordering) {
    // Create the local blocks
    std::vector<costa::block<T> > loc_blks;

    // Create blocks
    for (int i = 0; i < nlocalblocks; ++i) {
        auto &block = localblocks[i];
        auto row = block.row;
        auto col = block.col;
        auto ptr = reinterpret_cast<T *>(block.data);
        auto stride = block.ld;

        costa::block_coordinates coord{row, col};
        costa::interval rows{rowsplit[row],
                             rowsplit[row + 1]};
        costa::interval cols{colsplit[col],
                             colsplit[col + 1]};
        loc_blks.emplace_back(rows, cols, coord, ptr, stride);
    }

    auto grid = custom_grid(rowblocks, colblocks, rowsplit, colsplit, owners);
    return grid_layout<T>{std::move(grid), std::move(loc_blks), ordering};
}

costa::assigned_grid2D costa::block_cyclic_grid(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start (1-based, due to scalapack)
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int rsrc, const int csrc // coordinates of ranks oweing 
                                   // the first row (0-based)
        ) {
    // get rank grid ordering
    char rank_ordering = std::toupper(rank_grid_ordering);
    assert(rank_ordering == 'R' || rank_ordering == 'C');
    auto ordering = costa::scalapack::ordering::row_major;

    if (rank_ordering == 'C') {
        ordering = costa::scalapack::ordering::column_major;
    }

    auto scalapack_grid = costa::get_scalapack_grid(
        {m, n},
        {i, j},
        {sub_m, sub_n},
        {block_m, block_n},
        {proc_m, proc_n},
        ordering,
        {rsrc, csrc});

    return scalapack_grid;
}

template <typename T>
costa::grid_layout<T> costa::block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start (1-based,
                                  // due to scalapack convention)
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int rsrc, const int csrc, // coordinates of ranks oweing 
                                    // the first row (0-based)
        T* ptr, // local data of matrix A (not the submatrix)
        const int lld, // leading dimension
        const char data_ordering, // whether each local block is row- or col-major
        const int rank // processor rank
        ) {

    // get rank grid ordering
    char r_ordering = std::toupper(rank_grid_ordering);
    assert(r_ordering == 'R' || r_ordering == 'C');
    auto rank_ordering = costa::scalapack::ordering::row_major;
    if (r_ordering == 'C') {
        rank_ordering = costa::scalapack::ordering::column_major;
    }

    // get local blocks data ordering
    char d_ordering = std::toupper(data_ordering);
    assert(d_ordering == 'R' || d_ordering == 'C');

    auto scalapack_layout = costa::get_scalapack_layout<T>(
        lld,
        {m, n},
        {i, j},
        {sub_m, sub_n},
        {block_m, block_n},
        {proc_m, proc_n},
        rank_ordering,
        {rsrc, csrc},
        ptr,
        d_ordering,
        rank);

    return scalapack_layout;
}

// template instantiation for custom_layout
template 
costa::grid_layout<float> costa::custom_layout(
                             const int rowblocks,
                             const int colblocks,
                             const int* rowsplit,
                             const int* colsplit,
                             const int* owners,
                             const int nlocalblocks,
                             const block_t* localblocks,
                             const char data_ordering);
template 
costa::grid_layout<double> costa::custom_layout(
                             const int rowblocks,
                             const int colblocks,
                             const int* rowsplit,
                             const int* colsplit,
                             const int* owners,
                             const int nlocalblocks,
                             const block_t* localblocks,
                             const char data_ordering);
template 
costa::grid_layout<std::complex<float>> costa::custom_layout(
                             const int rowblocks,
                             const int colblocks,
                             const int* rowsplit,
                             const int* colsplit,
                             const int* owners,
                             const int nlocalblocks,
                             const block_t* localblocks,
                             const char data_ordering);
template 
costa::grid_layout<std::complex<double>> costa::custom_layout(
                             const int rowblocks,
                             const int colblocks,
                             const int* rowsplit,
                             const int* colsplit,
                             const int* owners,
                             const int nlocalblocks,
                             const block_t* localblocks,
                             const char data_ordering);

template
costa::grid_layout<float> costa::block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int ia, const int ja, // coordinates of ranks oweing 
                                    // the first row 
                                    // (1-based, scalapack-compatible)
        float* ptr, // local data of matrix A (not the submatrix)
        const int lld, // local leading dimension
        const char data_ordering, // whether each local block is row- or col-major
        const int rank // processor rank
);

template
costa::grid_layout<double> costa::block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int ia, const int ja, // coordinates of ranks oweing 
                                    // the first row 
                                    // (1-based, scalapack-compatible)
        double* ptr, // local data of matrix A (not the submatrix)
        const int lld, // local leading dimension
        const char data_ordering, // whether each local block is row- or col-major
        const int rank // processor rank
);

template
costa::grid_layout<std::complex<float>> costa::block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int ia, const int ja, // coordinates of ranks oweing 
                                    // the first row 
                                    // (1-based, scalapack-compatible)
        std::complex<float>* ptr, // local data of matrix A (not the submatrix)
        const int lld, // local leading dimension
        const char data_ordering, // whether each local block is row- or col-major
        const int rank // processor rank
);

template
costa::grid_layout<std::complex<double>> costa::block_cyclic_layout(
        const int m, const int n, // global matrix dimensions
        const int block_m, const int block_n, // block dimensions
        const int i, const int j, // submatrix start
        const int sub_m, const int sub_n, // submatrix size
        const int proc_m, const int proc_n, // processor grid dimension
        const char rank_grid_ordering, // rank grid ordering ('R' or 'C')
        const int ia, const int ja, // coordinates of ranks oweing 
                                    // the first row 
                                    // (1-based, scalapack-compatible)
        std::complex<double>* ptr, // local data of matrix A (not the submatrix)
        const int lld, // local leading dimension
        const char data_ordering, // whether each local block is row- or col-major
        const int rank // processor rank
);

