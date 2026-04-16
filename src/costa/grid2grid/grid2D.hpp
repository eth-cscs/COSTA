#pragma once
#include <costa/grid2grid/interval.hpp>

#include <vector>

namespace costa {
/*
A class describing a matrix n_rows x n_cols split into an arbitrary grid.
A grid is defined by rows_split and cols_split. More precisely, a matrix is
split into blocks such that a block with coordinates (i, j) is defined by:
    - row interval [rows_split[i], rows_split[i+1]) and similarly
    - column interval [cols_split[i], cols_split[i+1])
For each i, the following must hold:
    - 0 <= rows_split[i] < n_rows
    - 0 <= cols_split[i] < n_cols
*/
struct grid2D {
    // number of blocks in a row
    int n_rows = 0;
    // number of blocks in a column
    int n_cols = 0;
    // defines how rows are split
    std::vector<int> rows_split;
    // defines how cols are split
    std::vector<int> cols_split;

    grid2D() = default;
    grid2D(std::vector<int> &&r_split, std::vector<int> &&c_split);

    // returns index-th row interval, i.e. [rows_split[index], rows_split[index
    // + 1])
    interval row_interval(int index) const;

    // returns index-th column interval, i.e. [cols_split[index],
    // cols_split[index + 1])
    interval col_interval(int index) const;

    void transpose();
};

/*
A class describing a matrix split into an arbitrary grid (grid2D) where each
block is assigned to an arbitrary MPI rank. More precisely, a block with
coordinates (i, j) is assigned to an MPI rank defined by ranks[i][j].
*/
class assigned_grid2D {
  public:
    assigned_grid2D() = default;
    assigned_grid2D(grid2D &&g,
                    std::vector<std::vector<int>> &&proc,
                    int n_ranks);

    // returns the rank owning block (i, j)
    int owner(int i, int j) const;

    // returns a grid
    const grid2D &grid() const noexcept;

    // returns a total number of MPI ranks in the communicator
    // which owns the whole matrix. However, not all ranks have
    // to own a block of the matrix. Thus, it might happen that
    // maximum(ranks[i][j]) < n_ranks - 1 over all i, j
    int num_ranks() const noexcept;

    // returns index-th row interval
    interval rows_interval(int index) const;

    // returns index-th column interval
    interval cols_interval(int index) const;

    int block_size(int row_index, int col_index);

    void transpose();

    void reorder_ranks(std::vector<int>& reordering);

    int reordered_rank(int rank) const;

    bool ranks_reordered() const;

    friend std::ostream &operator<<(std::ostream &os, const assigned_grid2D &other) {
        for (int i = 0; i < other.num_blocks_row(); ++i) {
            for (int j = 0; j < other.num_blocks_col(); ++j) {
                os << "block " << other.rows_interval(i) << " x " << other.cols_interval(j) <<  " is owned by rank " 
                   << other.owner(i, j) << std::endl;
            }
        }
        return os;
    }

    // returns the number of rows/cols in the grid
    // i.e. the number of blocks in a row/col of the grid
    int num_blocks_row() const noexcept { 
        return g.n_rows;
    }

    int num_blocks_col() const noexcept {
        return g.n_cols;
    }

    // returns the number of rows/cols in the matrix
    int num_rows() const noexcept { 
        return g.rows_split.back();
    }

    int num_cols() const noexcept { 
        return g.cols_split.back();
    }

  private:
    friend bool operator==(assigned_grid2D const &,
                           assigned_grid2D const &) noexcept;

    bool transposed = false;

    grid2D g;
    std::vector<std::vector<int>> ranks;
    int n_ranks = 0;

    std::vector<int> ranks_reordering;
};

inline bool operator==(assigned_grid2D const &lhs,
                assigned_grid2D const &rhs) noexcept {
    return lhs.g.rows_split == rhs.g.rows_split &&
           lhs.g.cols_split == rhs.g.cols_split && lhs.ranks == rhs.ranks;
}

/*
    A class describing a matrix n_rows x n_cols split into an arbitrary grid.
    A grid is defined by rows_split and cols_split. More precisely, a matrix is
   split into blocks such that a block with coordinates (i, j) is defined by:
        - row interval [rows_split[i], rows_split[i+1]) and similarly
        - column interval [cols_split[i], cols_split[i+1])
    For each i, the following must hold:
        - 0 <= rows_split[i] < n_rows
        - 0 <= cols_split[i] < n_cols
*/
inline grid2D::grid2D(std::vector<int> &&r_split, std::vector<int> &&c_split)
    : n_rows(r_split.size() ? r_split.size() - 1 : 0)
    , n_cols(c_split.size() ? c_split.size() - 1 : 0)
    , rows_split(std::forward<std::vector<int>>(r_split))
    , cols_split(std::forward<std::vector<int>>(c_split)) {}

// returns index-th row interval, i.e. [rows_split[index], rows_split[index +
// 1])
inline interval grid2D::row_interval(int index) const {
    if ((size_t)index >= rows_split.size() - 1) {
        throw std::runtime_error(
            "ERROR: in class grid2D, row index out of range.");
    }
    return {rows_split[index], rows_split[index + 1]};
}

// returns index-th column interval, i.e. [cols_split[index], cols_split[index +
// 1])
inline interval grid2D::col_interval(int index) const {
    if ((size_t)index >= cols_split.size() - 1) {
        throw std::runtime_error(
            "ERROR: in class grid2D, col index out of range.");
    }
    return {cols_split[index], cols_split[index + 1]};
}

inline void grid2D::transpose() {
    std::swap(rows_split, cols_split);
    std::swap(n_rows, n_cols);
}

/*
A class describing a matrix split into an arbitrary grid (grid2D) where each
block is assigned to an arbitrary MPI rank. More precisely, a block with
coordinates (i, j) is assigned to an MPI rank defined by ranks[i][j].
*/
inline assigned_grid2D::assigned_grid2D(grid2D &&g,
                                 std::vector<std::vector<int>> &&proc,
                                 int n_ranks)
    : g(std::forward<grid2D>(g))
    , ranks(std::forward<std::vector<std::vector<int>>>(proc))
    , n_ranks(n_ranks) {}

// returns the rank owning block (i, j)
inline int assigned_grid2D::owner(int i, int j) const {
    int b_row = transposed ? j : i;
    int b_col = transposed ? i : j;
    return reordered_rank(ranks[b_row][b_col]);
}

// returns a grid
inline const grid2D &assigned_grid2D::grid() const noexcept { return g; }

// returns a total number of MPI ranks in the communicator
// which owns the whole matrix. However, not all ranks have
// to own a block of the matrix. Thus, it might happen that
// maximum(ranks[i][j]) < n_ranks - 1 over all i, j
inline int assigned_grid2D::num_ranks() const noexcept { return n_ranks; }

// returns index-th row interval
inline interval assigned_grid2D::rows_interval(int index) const {
    return g.row_interval(index);
}

// returns index-th column interval
inline interval assigned_grid2D::cols_interval(int index) const {
    return g.col_interval(index);
}

inline int assigned_grid2D::block_size(int row_index, int col_index) {
    return rows_interval(row_index).length() *
           cols_interval(col_index).length();
}

inline void assigned_grid2D::transpose() {
    g.transpose();
    // ranks = transpose(ranks);
    transposed = !transposed;
}

inline void assigned_grid2D::reorder_ranks(std::vector<int>& reordering) {
    ranks_reordering = reordering;
}

inline int assigned_grid2D::reordered_rank(int rank) const {
    assert(rank < std::max((int) ranks_reordering.size(), n_ranks));
    if (ranks_reordered())
        return ranks_reordering[rank];
    else
        return rank;
}

inline bool assigned_grid2D::ranks_reordered() const {
    return ranks_reordering.size() > 0;
}

} // namespace costa
