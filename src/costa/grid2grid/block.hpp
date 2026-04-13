#pragma once
#include <costa/grid2grid/grid2D.hpp>
#include <costa/grid2grid/interval.hpp>

#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

namespace costa {

[[gnu::always_inline]] inline int conjugate_f(int el) { return el; }

[[gnu::always_inline]] inline double conjugate_f(double el) { return el; }

[[gnu::always_inline]] inline float conjugate_f(float el) { return el; }

[[gnu::always_inline]] inline std::complex<float> conjugate_f(std::complex<float> el) {
    return std::conj(el);
}

[[gnu::always_inline]] inline std::complex<double> conjugate_f(std::complex<double> el) {
    return std::conj(el);
}

struct block_coordinates {
    int row = 0;
    int col = 0;
    block_coordinates() = default;
    block_coordinates(int r, int c);

    [[gnu::always_inline]] void transpose() { std::swap(row, col); }
};

struct block_range {
    interval rows_interval;
    interval cols_interval;

    block_range() = default;

    [[gnu::always_inline]] block_range(interval r, interval c)
        : rows_interval(r), cols_interval(c) {}

    [[gnu::always_inline]] bool outside_of(const block_range &range) const {
        return (rows_interval.end <= range.rows_interval.start ||
                rows_interval.start >= range.rows_interval.end) &&
               (cols_interval.end <= range.cols_interval.start ||
                cols_interval.end <= range.cols_interval.start);
    }

    [[gnu::always_inline]] bool inside(const block_range &range) const {
        return range.rows_interval.start < rows_interval.start &&
               range.rows_interval.end > rows_interval.end &&
               range.cols_interval.start < cols_interval.start &&
               range.cols_interval.end > cols_interval.end;
    }

    [[gnu::always_inline]] bool intersects(const block_range &range) const {
        return !outside_of(range) && !inside(range);
    }

    [[gnu::always_inline]] block_range intersection(const block_range &other) const {
        return {rows_interval.intersection(other.rows_interval),
                cols_interval.intersection(other.cols_interval)};
    }

    [[gnu::always_inline]] bool non_empty() const {
        return rows_interval.non_empty() && cols_interval.non_empty();
    }

    [[gnu::always_inline]] bool empty() const {
        return rows_interval.empty() || cols_interval.empty();
    }

    [[gnu::always_inline]] bool operator==(const block_range &other) const {
        if (empty()) return other.empty();
        return rows_interval == other.rows_interval &&
               cols_interval == other.cols_interval;
    }

    [[gnu::always_inline]] bool operator!=(const block_range &other) const {
        return !(*this == other);
    }
};

inline std::ostream &operator<<(std::ostream &os, const block_range &other) {
    os << "rows:" << other.rows_interval << ", cols:" << other.cols_interval
       << std::endl;
    return os;
}

// assumes column-major ordering inside block
template <typename T>
struct block {
    // blocks from different matrices
    // should have different tags
    int tag = 0;
    // start and end index of the block
    interval rows_interval;
    interval cols_interval;

    // std::optional<T> scalar = std::nullopt;

    block_coordinates coordinates;

    T *data = nullptr;
    int stride = 0;

    char _ordering = 'C';

    bool transposed = false;

    block() = default;

    block(const assigned_grid2D &grid,
          block_coordinates coord,
          T *ptr,
          const int stride);

    block(const assigned_grid2D &grid, block_range &range, T *ptr,
          const int stride);

    block(const assigned_grid2D &grid,
          interval r_inter,
          interval c_inter,
          T *ptr,
          const int stride);

    block(interval r_inter,
          interval c_inter,
          block_coordinates coord,
          T *ptr,
          const int stride);

    block(block_range &range, block_coordinates coord, T *ptr, 
          const int stride);;

    // finds the index of the interval inter in splits
    int interval_index(const std::vector<int> &splits, interval inter);

    block<T> subblock(interval r_range, interval c_range) const;

    bool non_empty() const;

    // implementing comparator
    bool operator<(const block &other) const;

    int n_rows() const { return rows_interval.length(); }

    int n_cols() const { return cols_interval.length(); }

    std::pair<int, int> size() const { return {n_rows(), n_cols()}; }

    size_t total_size() const { return n_rows() * n_cols(); }

    T local_element(int li, int lj) const;
    T& local_element(int li, int lj);

    std::pair<int, int> local_to_global(int li, int lj) const;
    std::pair<int, int> global_to_local(int gi, int gj) const;

    void transpose();
    void set_ordering(const char ordering);

    // scales the local block by beta
    void scale_by(T beta);
    // fill local block with beta
    void fill(T beta);
};

template <typename T>
bool operator==(block<T> const &lhs, block<T> const &rhs) noexcept {
    return lhs.rows_interval.start == rhs.rows_interval.start &&
           lhs.rows_interval.end == rhs.rows_interval.end &&
           lhs.cols_interval.start == rhs.cols_interval.start &&
           lhs.cols_interval.end == rhs.cols_interval.end &&
           lhs.coordinates.row == rhs.coordinates.row &&
           lhs.coordinates.col == rhs.coordinates.col &&
           lhs.stride == rhs.stride && lhs.data == rhs.data &&
           lhs._ordering == rhs._ordering &&
           lhs.tag == rhs.tag;
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const block<T> &other) {
    return os << "rows: " << other.rows_interval
              << "cols: " << other.cols_interval << std::endl;
}

template <typename T>
class local_blocks {
  public:
    local_blocks() = default;
    local_blocks(std::vector<block<T>> &&blocks);

    block<T> &get_block(int i);

    const block<T> &get_block(int i) const { return blocks[i]; }

    int num_blocks() const;

    size_t size() const;

    void transpose();

  private:
    template <typename T_>
    friend bool operator==(local_blocks<T_> const &,
                           local_blocks<T_> const &) noexcept;

    std::vector<block<T>> blocks;
    size_t total_size = 0;
};

template <typename T>
bool operator==(local_blocks<T> const &lhs,
                local_blocks<T> const &rhs) noexcept {
    return lhs.blocks == rhs.blocks && lhs.total_size == rhs.total_size;
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os,
                                const local_blocks<T> &other) {
    for (unsigned i = 0; i < (unsigned)other.num_blocks(); ++i) {
        os << "block " << i << ":\n" << other.get_block(i) << std::endl;
    }
    return os;
}
} // namespace costa
