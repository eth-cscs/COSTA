#include <costa/grid2grid/communication_data.hpp>
#include <costa/grid2grid/workspace.hpp>

#include <complex>
#include <omp.h>

namespace costa {
// *********************
//     MESSAGE
// *********************
template <typename T>
message<T>::message(block<T> b, int rank,
                    char ordering,
                    T alpha, T beta,
                    bool trans, bool conj)
    : b(b)
    , rank(rank)
    , alpha(alpha)
    , beta(beta)
    , transpose(trans) {

    assert(b.non_empty());

    // check if the values are valid
    assert(ordering == 'R' || ordering == 'C');

    // set boolean variables based on given chars
    col_major = ordering == 'C';

    // check if the data type is complex
    bool is_complex = std::is_same<T, std::complex<double>>::value ||
                      std::is_same<T, std::complex<float>>::value;
    conjugate = conj && is_complex;
}

template <typename T>
std::string message<T>::to_string() const {
    std::string transposed = transpose ? "True" : "False";
    std::string conjugated = conjugate ? "True" : "False";
    std::string col_majored = col_major ? "True" : "False";
    std::string block_col_majored = b._ordering=='C' ? "Col-major" : "Row-major";


    std::string res = "";
    res += "Message: \n";
    res += "rank = " + std::to_string(rank) + "\n";
    res += "transpose = " + transposed + "\n";
    res += "conjugate = " + conjugated + "\n";
    res += "col_major = " + col_majored + "\n";
    res += "block: " + std::to_string(b.n_rows()) + " x " + std::to_string(b.n_cols()) + "\n";
    res += "block tag: " + std::to_string(b.tag) + "\n";
    res += "block ordering: " + block_col_majored + "\n";
    return res;
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
    int rank = get_rank();
    int other_rank = other.get_rank();

    const double alpha = std::abs(this->alpha);
    const double beta = std::abs(this->beta);

    const double other_alpha = std::abs(other.alpha);
    const double other_beta = std::abs(other.beta);

    return std::tie(rank, b, alpha, beta, transpose, conjugate)
           <
           std::tie(other_rank, other.b, other_alpha, other_beta, other.transpose, other.conjugate)
           ;
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
                                          int rank, int n_ranks,
					  CommType type)
    : n_ranks(n_ranks)
    , my_rank(rank)
    , type(type) {
    // std::cout << "constructor of communication data invoked" << std::endl;
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
            /*
            if (rank == 0) {
                std::cout << "message " << i << ":\n";
                std::cout << m.to_string() << std::endl;
            }
            */
            local_messages.push_back(m);
        }
    }
    /*
    if(rank == 0) {
        std::cout << "==========================" <<std::endl;
    }
    */
    memory::get_costa_context_instance<T>()->resize_buffer(type, total_size);  

    for (unsigned i = 1; i < (unsigned)n_ranks; ++i) {
        dspls[i] = dspls[i - 1] + counts[i - 1];
    }

    n_packed_messages = 0;
    for (unsigned i = 0; i < (unsigned) n_ranks; ++i) {
        if (counts[i] > 0) {
            ++n_packed_messages;
        }
    }

    partition_messages();
}

template <typename T>
void communication_data<T>::copy_from_buffer() {
    if (mpi_messages.size() > 0) {
	auto& workspace = *memory::get_costa_context_instance<T>();
#pragma omp parallel for shared(mpi_messages, offset_per_message, workspace)
        for (unsigned i = 0; i < mpi_messages.size(); ++i) {
            const auto &m = mpi_messages[i];
            block<T> b = m.get_block();
            bool b_col_major = b._ordering == 'C';
            // std::cout <<"To buffer: Stride = " << b.stride << " -> 0" << std::endl;
            int num_rows = b.n_rows();
            int num_cols = b.n_cols();
            if (b.transposed) std::swap(num_rows, num_cols);
            copy_and_transform(num_rows, num_cols,
                               data() + offset_per_message[i],
                               0, m.col_major,
                               b.data, b.stride, b_col_major,
                               m.transpose,
                               m.conjugate,
                               m.alpha, m.beta,
                               workspace);
        }
    }
}

template <typename T>
void communication_data<T>::copy_to_buffer() {
    if (mpi_messages.size() > 0) {
	auto& workspace = *memory::get_costa_context_instance<T>();
#pragma omp parallel for shared(mpi_messages, offset_per_message, workspace)
        for (unsigned i = 0; i < mpi_messages.size(); ++i) {
            const auto &m = mpi_messages[i];
            block<T> b = m.get_block();
            bool b_col_major = b._ordering == 'C';
            // std::cout <<"To buffer: Stride = " << b.stride << " -> 0" << std::endl;
            int num_rows = b.n_rows();
            int num_cols = b.n_cols();
            if (b.transposed) std::swap(num_rows, num_cols);
            copy_and_transform(num_rows, num_cols,
                               b.data, b.stride, b_col_major, 
                               // dest_ptr
                               data() + offset_per_message[i],
                               0,
                               b_col_major,
                               false, // no transpose on copy to buffer
                               false, // no conjugate on copy to buffer
                               T{1}, T{0},// no scaling on copy to buffer
                               workspace 
                               );
        }
    }
}

template <typename T>
void communication_data<T>::copy_from_buffer(int idx) {
    assert(idx >= 0 && idx+1 < package_ticks.size());
    if (package_ticks[idx+1] - package_ticks[idx] > 0) {
	auto& workspace = *memory::get_costa_context_instance<T>();
#pragma omp parallel for shared(idx, package_ticks, mpi_messages, offset_per_message, workspace)
        for (unsigned i = package_ticks[idx]; i < package_ticks[idx+1]; ++i) {
            const auto &m = mpi_messages[i];
            block<T> b = m.get_block();
            bool b_col_major = b._ordering == 'C';
            int num_rows = b.n_rows();
            int num_cols = b.n_cols();
            if (m.transpose) std::swap(num_rows, num_cols);
            // std::cout <<"From buffer: Stride = 0 -> " << b.stride << std::endl;
            // std::cout <<"Block ordering = " << b._ordering << std::endl;
            copy_and_transform(num_rows, num_cols,
                               data() + offset_per_message[i],
                               0, m.col_major,
                               b.data, b.stride, b_col_major,
                               m.transpose,
                               m.conjugate,
                               m.alpha, m.beta,
                               workspace);
        }
    }
}

template <typename T>
T *communication_data<T>::data() {
    return memory::get_costa_context_instance<T>()->buffer_ptr(type);
}

template <typename T>
void copy_local_blocks(std::vector<message<T>>& from,
                       std::vector<message<T>>& to) {
    assert(from.size() == to.size());
    if (from.size() > 0) {
	auto& workspace = *memory::get_costa_context_instance<T>();
#pragma omp parallel for shared(from, to, workspace)
        for (int i = 0; i < from.size(); ++i) {
            /*
            if (from[i].transpose != to[i].transpose) {
                std::cout << "from = " << from[i].to_string() <<", to = " << to[i].to_string() <<std::endl;
            }
            */
            assert(from[i].alpha == to[i].alpha);
            assert(from[i].beta == to[i].beta);
            assert(from[i].transpose == to[i].transpose);
            assert(from[i].conjugate == to[i].conjugate);
            assert(from[i].get_rank() == to[i].get_rank());

            auto b_src = from[i].get_block();
            auto b_dest = to[i].get_block();
            assert(b_src.non_empty());
            assert(b_dest.non_empty());
            /*
            if (b_src.total_size() != b_dest.total_size()) {
                std::cout << "bsrc size = " << b_src.total_size() <<", bdest size" << b_dest.total_size() <<std::endl;
            }
            */
            assert(b_src.total_size() == b_dest.total_size());
            assert(b_src.tag == b_dest.tag);
            assert(from[i].transpose == b_src.transposed);
            assert(!b_dest.transposed);

            bool b_src_col_major = b_src._ordering == 'C';
            bool b_dest_col_major = b_dest._ordering == 'C';
            // std::cout <<"src = " << b_src._ordering <<", dest = " << b_dest._ordering << std::endl;
            // std::cout <<"LOCAL: src stride = " << b_src.stride <<", dest stride = " << b_dest.stride << std::endl;
            int num_rows = b_src.n_rows();
            int num_cols = b_src.n_cols();
            if (b_src.transposed) std::swap(num_rows, num_cols);

            copy_and_transform(num_rows, num_cols,
                               b_src.data, b_src.stride, b_src_col_major,
                               b_dest.data,
                               b_dest.stride, b_dest_col_major,
                               from[i].transpose,
                               from[i].conjugate,
                               from[i].alpha, from[i].beta,
                               workspace);
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
