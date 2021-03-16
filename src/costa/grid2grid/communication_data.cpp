#include <costa/grid2grid/communication_data.hpp>

#include <complex>
#include <omp.h>

namespace costa {
// *********************
//     MESSAGE
// *********************
template <typename T>
message<T>::message(block<T> b, int rank,
        T alpha, T beta,
        char trans, char ordering)
    : b(b)
    , rank(rank)
    , alpha(alpha)
    , beta(beta) {

    // ignore upper/lower-case
    trans = std::toupper(trans);
    ordering = std::toupper(ordering);

    // check if the values are valid
    assert(trans == 'N' || trans == 'T' || trans == 'C');
    assert(ordering == 'R' || ordering == 'C');

    // set boolean variables based on given chars
    transpose = trans != 'N';
    col_major = ordering == 'C';

    // check if the data type is complex
    bool is_complex = std::is_same<T, std::complex<double>>::value ||
                      std::is_same<T, std::complex<float>>::value;
    conjugate = trans == 'C' && is_complex;
}

template <typename T>
block<T> message<T>::get_block() const {
    return b;
}

template <typename T>
int message<T>::get_rank() const {
    return rank;
}

// implementing comparator
template <typename T>
bool message<T>::operator<(const message<T> &other) const {
    return get_rank() < other.get_rank() ||
           (get_rank() == other.get_rank() && b < other.get_block()); 
}

template <typename T>
void communication_data<T>::partition_messages() {
    if (mpi_messages.size() == 0) 
        return;

    int pivot = -1; 
    for (int i = 0; i < mpi_messages.size(); ++i) {
        int rank = mpi_messages[i].get_rank();
        if (pivot != rank) {
            pivot = rank;
            package_ticks.push_back(i);
        }
    }
    package_ticks.push_back(mpi_messages.size());
}

// ************************
//   COMMUNICATION DATA
// ************************
template <typename T>
communication_data<T>::communication_data(std::vector<message<T>> &messages,
                                          int rank, int n_ranks)
    : n_ranks(n_ranks)
    , my_rank(rank) {
    // std::cout << "constructor of communciation data invoked" << std::endl;
    dspls = std::vector<int>(n_ranks);
    counts = std::vector<int>(n_ranks);
    mpi_messages.reserve(messages.size());
    offset_per_message.reserve(messages.size());

    int offset = 0;

    int prev_rank = -1;

    for (unsigned i = 0; i < messages.size(); ++i) {
        const auto &m = messages[i];
        int target_rank = m.get_rank();
        block<T> b = m.get_block();
        assert(b.non_empty());

        // if the message should be communicated to 
        // a different rank
        if (target_rank != my_rank) {
            mpi_messages.push_back(m);
            offset_per_message.push_back(offset);
            offset += b.total_size();
            counts[target_rank] += b.total_size();
            total_size += b.total_size();
            prev_rank = target_rank;
        } else {
            local_messages.push_back(m);
        }
    }

    buffer = std::unique_ptr<T[]>(new T[total_size]);
    for (unsigned i = 1; i < (unsigned)n_ranks; ++i) {
        dspls[i] = dspls[i - 1] + counts[i - 1];
    }

    n_packed_messages = 0;
    for (unsigned i = 0; i < (unsigned) n_ranks; ++i) {
        if (counts[i]) {
            ++n_packed_messages;
        }
    }

    partition_messages();
}

template <typename T>
void communication_data<T>::copy_to_buffer() {
    if (mpi_messages.size()) {
#pragma omp parallel for schedule(dynamic, 1)
        for (unsigned i = 0; i < mpi_messages.size(); ++i) {
            const auto &m = mpi_messages[i];
            block<T> b = m.get_block();

            copy_and_transform(b.n_rows(), b.n_cols(),
                               b.data, b.stride, b.ordering, 
                               // dest_ptr
                               data() + offset_per_message[i],
                               0,
                               m.col_major,
                               false, // no transpose on copy to buffer
                               false, // no conjugate on copy to buffer
                               tiling,
                               T{1}, T{0} // no scaling on copy to buffer
                               );
        }
    }
}

template <typename T>
void communication_data<T>::copy_from_buffer(int idx) {
    assert(idx >= 0 && idx+1 < package_ticks.size());
    if (package_ticks[idx+1] - package_ticks[idx]) {
#pragma omp parallel for schedule(dynamic, 1)
        for (unsigned i = package_ticks[idx]; i < package_ticks[idx+1]; ++i) {
            const auto &m = mpi_messages[i];
            block<T> b = m.get_block();
            copy_and_transform(b.n_rows(), b.n_cols(),
                               b.data, 0, m.col_major, 
                               data() + offset_per_message[i],
                               b.stride, b.ordering,
                               m.transpose,
                               m.conjugate,
                               tiling,
                               m.alpha, m.beta);
        }
    }
}

template <typename T>
T *communication_data<T>::data() {
    return buffer.get();
}

template <typename T>
void copy_local_blocks(std::vector<message<T>>& from,
                       std::vector<message<T>>& to) {
    assert(from.size() == to.size());
    if (from.size() > 0) {
    memory::tiling_manager<T> tiling;
#pragma omp parallel for 
        for (unsigned i = 0u; i < from.size(); ++i) {
            assert(from[i].alpha == to[i].alpha);
            assert(from[i].beta == to[i].beta);
            assert(from[i].transpose == to[i].transpose);
            assert(from[i].conjugate == to[i].conjugate);
            assert(from[i].get_rank() == to[i].get_rank());

            auto b_src = from[i].get_block();
            auto b_dest = to[i].get_block();
            assert(b_src.non_empty());
            assert(b_dest.non_empty());
            assert(b_src.total_size() == b_dest.total_size());

            copy_and_transform(b_src.n_rows(), b_src.n_cols(),
                               b_src.data, b_src.stride, b_src.ordering,
                               b_dest.data,
                               b_dest.stride, b_dest.ordering,
                               from[i].transpose,
                               from[i].conjugate,
                               tiling,
                               from[i].alpha, from[i].beta);
    }
    }
}

// template instantiation for communication_data
template class communication_data<double>;
template class communication_data<std::complex<double>>;
template class communication_data<float>;
template class communication_data<std::complex<float>>;

// template instantiation for message
template class message<double>;
template class message<std::complex<double>>;
template class message<float>;
template class message<std::complex<float>>;

// template instantiation for copy_local_blocks
template void
copy_local_blocks(std::vector<message<double>>& from, 
                  std::vector<message<double>>& to);
template void
copy_local_blocks(std::vector<message<float>>& from, 
                  std::vector<message<float>>& to);
template void
copy_local_blocks(std::vector<message<std::complex<float>>>& from, 
                  std::vector<message<std::complex<float>>>& to);
template void
copy_local_blocks(std::vector<message<std::complex<double>>>& from, 
                  std::vector<message<std::complex<double>>>& to);
} // namespace costa
