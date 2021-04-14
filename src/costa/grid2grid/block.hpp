#pragma once
#include <costa/grid2grid/grid2D.hpp>
#include <costa/grid2grid/interval.hpp>

#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

namespace costa {

int conjugate_f(int el);

double conjugate_f(double el);

float conjugate_f(float el);

std::complex<float> conjugate_f(std::complex<float> el);

std::complex<double> conjugate_f(std::complex<double> el);

struct block_coordinates {
    int row = 0;
    int col = 0;
    block_coordinates() = default;
    block_coordinates(int r, int c);

    void transpose();
};

struct block_range {
    interval rows_interval;
    interval cols_interval;

    block_range() = default;
    block_range(interval r, interval c);

    bool outside_of(const block_range &range) const;

    bool inside(const block_range &range) const;

    bool intersects(const block_range &range) const;

    block_range intersection(const block_range &other) const;

    bool non_empty() const;

    bool empty() const;

    bool operator==(const block_range &other) const;

    bool operator!=(const block_range &other) const;
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
