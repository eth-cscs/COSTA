#pragma once
#include <costa/grid2grid/block.hpp>
#include <costa/grid2grid/grid2D.hpp>
#include <costa/grid2grid/mpi_type_wrapper.hpp>

namespace costa {
template <typename T>
class grid_layout {
  public:
    grid_layout() = default;

    grid_layout(assigned_grid2D &&g, local_blocks<T> &&b)
        : grid(std::forward<assigned_grid2D>(g))
        , blocks(std::forward<local_blocks<T>>(b)) {}

    int num_ranks() const { return grid.num_ranks(); }

    void transpose_or_conjugate(char flag) {
        flag = std::toupper(flag);
        assert(flag == 'N' || flag == 'T' || flag == 'C');
        if (flag == 'T' || flag == 'C') {
            grid.transpose();
            blocks.transpose_or_conjugate(flag);
        }
    }

    void reorder_ranks(std::vector<int>& reordering) {
        grid.reorder_ranks(reordering);
    }

    int reordered_rank(int rank) const {
        return grid.reordered_rank(rank);
    }

    bool ranks_reordered() {
        return grid.ranks_reordered();
    }

    // returns the number of rows and column in the matrix
    // that this grid represents, and not the number of blocks
    // in a row or column of the grid
    int num_cols() const noexcept { return grid.num_cols(); }
    int num_rows() const noexcept { return grid.num_rows(); }

    // returns the number of blocks in a row/column of the grid
    int num_blocks_col() const noexcept { return grid.num_blocks_col(); }
    int num_blocks_row() const noexcept { return grid.num_blocks_row(); }

    // scales all local blocks by beta
    void scale_by(const T beta) {
        if (beta == T{1}) return;
        // iterate over all local blocks
        // local_blocks contains only blocks within the specified submatrix
        for (unsigned i = 0u; i < blocks.num_blocks(); ++i) {
            auto& block = blocks.get_block(i);
            block.scale_by(beta);
        }
    }

    assigned_grid2D grid;
    local_blocks<T> blocks;
    // std::optional<T> scalar = std::nullopt;
};

template <typename T>
using layout_ref = std::reference_wrapper<grid_layout<T>>;

} // namespace costa
