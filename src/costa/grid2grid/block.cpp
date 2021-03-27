#include <costa/grid2grid/block.hpp>

#include <complex>

namespace costa {
int conjugate_f(int el) {
    return el; 
}

double conjugate_f(double el) {
    return el; 
}

float conjugate_f(float el) {
    return el; 
}

std::complex<float> conjugate_f(std::complex<float> el) {
    return std::conj(el); 
}

std::complex<double> conjugate_f(std::complex<double> el) {
    return std::conj(el); 
}

block_coordinates::block_coordinates(int r, int c)
    : row(r)
    , col(c) {}

void block_coordinates::transpose() {
    std::swap(row, col);
}

block_range::block_range(interval r, interval c)
    : rows_interval(r)
    , cols_interval(c) {}

bool block_range::outside_of(const block_range &range) const {
    return (rows_interval.end <= range.rows_interval.start ||
            rows_interval.start >= range.rows_interval.end) &&
           (cols_interval.end <= range.cols_interval.start ||
            cols_interval.end <= range.cols_interval.start);
}

bool block_range::inside(const block_range &range) const {
    return range.rows_interval.start < rows_interval.start &&
           range.rows_interval.end > rows_interval.end &&
           range.cols_interval.start < cols_interval.start &&
           range.cols_interval.end > cols_interval.end;
}

bool block_range::intersects(const block_range &range) const {
    return !outside_of(range) && !inside(range);
}

block_range block_range::intersection(const block_range &other) const {
    interval rows_intersection =
        rows_interval.intersection(other.rows_interval);
    interval cols_intersection =
        cols_interval.intersection(other.cols_interval);
    return {rows_intersection, cols_intersection};
}

bool block_range::non_empty() const {
    return rows_interval.non_empty() && cols_interval.non_empty();
}

bool block_range::empty() const {
    return rows_interval.empty() || cols_interval.empty();
}

bool block_range::operator==(const block_range &other) const {
    if (empty()) {
        return other.empty();
    }
    return rows_interval == other.rows_interval &&
           cols_interval == other.cols_interval;
}

bool block_range::operator!=(const block_range &other) const {
    return !(*this == other);
}

template <typename T>
block<T>::block(const assigned_grid2D &grid,
                block_coordinates coord,
                T *ptr,
                const int stride)
    : rows_interval(grid.rows_interval(coord.row))
    , cols_interval(grid.cols_interval(coord.col))
    , coordinates(coord)
    , data(ptr)
    , stride(stride)
{}

template <typename T>
block<T>::block(const assigned_grid2D &grid,
                block_range &range,
                T *ptr,
                const int stride)
    : block(grid, range.rows_interval, range.cols_interval, ptr, stride) {}

template <typename T>
block<T>::block(const assigned_grid2D &grid,
                interval r_inter,
                interval c_inter,
                T *ptr,
                const int stride)
    : rows_interval(r_inter)
    , cols_interval(c_inter)
    , data(ptr)
    , stride(stride)
{
    // compute the coordinates based on the grid and intervals
    int row_coord = interval_index(grid.grid().rows_split, rows_interval);
    int col_coord = interval_index(grid.grid().cols_split, cols_interval);
    coordinates = block_coordinates(row_coord, col_coord);
}

template <typename T>
block<T>::block(interval r_inter,
                interval c_inter,
                block_coordinates coord,
                T *ptr,
                const int stride)
    : rows_interval(r_inter)
    , cols_interval(c_inter)
    , coordinates(coord)
    , data(ptr)
    , stride(stride)
{}

template <typename T>
block<T>::block(block_range &range, block_coordinates coord, T *ptr, 
                const int stride)
    : block(range.rows_interval, range.cols_interval, coord, ptr, stride) {}


// finds the index of the interval inter in splits
template <typename T>
int block<T>::interval_index(const std::vector<int> &splits, interval inter) {
    auto ptr = std::lower_bound(splits.begin(), splits.end(), inter.start);
    int index = std::distance(splits.begin(), ptr);
    return index;
}

template <typename T>
block<T> block<T>::subblock(interval r_range, interval c_range) const {
    if (!rows_interval.contains(r_range) || !cols_interval.contains(c_range)) {
        std::cout << "BLOCK: row_interval = " << rows_interval
                  << ", column_interval = " << cols_interval << std::endl;
        std::cout << "SUBBLOCK: row_interval = " << r_range
                  << ", column_interval = " << c_range << std::endl;
        throw std::runtime_error(
            "ERROR: current block does not contain requested subblock.");
    }
    // column-major ordering inside block assumed here
    auto r_interval = rows_interval;
    auto c_interval = cols_interval;
    auto coord = coordinates;

    T *ptr = data + (c_range.start - c_interval.start) * stride +
             (r_range.start - r_interval.start);
    // std::cout << "stride = " << stride << std::endl;
    // std::cout << "ptr offset = " << (ptr - data) << std::endl;
    block<T> b(r_range, c_range, coord, ptr, stride); // correct
    b._ordering = _ordering;
    b.tag = tag;
    return b;
}

template <typename T>
bool block<T>::non_empty() const {
    bool non_empty_intervals =
        cols_interval.non_empty() && rows_interval.non_empty();
    assert(!non_empty_intervals || data);
    return non_empty_intervals;
}

template <typename T>
bool block<T>::operator<(const block &other) const {
    bool tags_less = tag < other.tag;
    bool tags_equal = tag == other.tag;
    bool cols_less = cols_interval < other.cols_interval;
    bool cols_equal = cols_interval == other.cols_interval;
    bool rows_less = rows_interval < other.rows_interval;
    bool rows_equal = rows_interval == other.rows_interval;

    bool blocks_less = cols_less || 
                       (cols_equal && rows_less);
    bool blocks_equal = cols_equal && rows_equal;

    return blocks_less || (blocks_equal && tags_less);

    // return cols_interval.start < other.cols_interval.start ||
    //        (cols_interval.start == other.cols_interval.start &&
    //         rows_interval.start < other.rows_interval.start);
}

template <typename T>
T block<T>::local_element(int li, int lj) const {
    assert(li >= 0 && li < n_rows());
    assert(lj >= 0 && lj < n_cols());
    assert(_ordering == 'C' || _ordering == 'R');

    int offset = stride * lj + li;
    if (_ordering == 'R') {
        offset = stride * li + lj;
    }
    return data[offset];
}

template <typename T>
T& block<T>::local_element(int li, int lj) {
    assert(li >= 0 && li < n_rows());
    assert(lj >= 0 && lj < n_cols());
    assert(_ordering == 'C' || _ordering == 'R');

    int offset = stride * lj + li;
    if (_ordering == 'R') {
        offset = stride * li + lj;
    }
    return data[offset];
}

template <typename T>
std::pair<int, int> block<T>::local_to_global(int li, int lj) const {
    assert(li >= 0 && li < n_rows());
    assert(lj >= 0 && lj < n_cols());

    int gi = rows_interval.start + li;
    int gj = cols_interval.start + lj;

    return std::pair<int, int>{gi, gj};
}

template <typename T>
std::pair<int, int> block<T>::global_to_local(int gi, int gj) const {
    int li = -1;
    int lj = -1;

    if (rows_interval.contains(gi)) {
        li = gi - rows_interval.start;
    }
    if (cols_interval.contains(gj)) {
        lj = gj - cols_interval.start;
    }

    return std::pair<int, int>{li, lj};
}

// transpose local block
template <typename T>
void block<T>::transpose() {
    std::swap(rows_interval, cols_interval);
    coordinates.transpose();
}

template <typename T>
void block<T>::set_ordering(const char ordering) {
    _ordering = std::toupper(ordering);
    assert(_ordering == 'R' || _ordering == 'C');
}

template <typename T>
void block<T>::scale_by(T beta) {
    if (beta == T{1}) return;

    int num_rows = n_rows();
    int num_cols = n_cols();

    for (int lj = 0; lj < num_cols; ++lj) {
        for (int li = 0; li < num_rows; ++li) {
            int offset = stride * lj + li;
            data[offset] *= beta;
        }
    }
}


template <typename T>
local_blocks<T>::local_blocks(std::vector<block<T>> &&blocks)
    : blocks(std::forward<std::vector<block<T>>>(blocks)) {
    for (const auto &b : blocks) {
        this->total_size += b.total_size();
    }
}

template <typename T>
block<T> &local_blocks<T>::get_block(int i) {
    return blocks[i];
}

template <typename T>
int local_blocks<T>::num_blocks() const {
    return blocks.size();
}

template <typename T>
size_t local_blocks<T>::size() const {
    return total_size;
}

template <typename T>
void local_blocks<T>::transpose() {
    for (auto& b: blocks) {
        b.transpose();
    }
}

template struct block<double>;
template struct block<std::complex<double>>;
template struct block<float>;
template struct block<std::complex<float>>;

template class local_blocks<double>;
template class local_blocks<std::complex<double>>;
template class local_blocks<float>;
template class local_blocks<std::complex<float>>;

} // end namespace costa
