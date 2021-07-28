#include <costa/grid2grid/transform.hpp>
#include <costa/grid2grid/profiler.hpp>
#include <costa/grid2grid/utils.hpp>
#include <assert.h>
#include <complex>
#include <thread>

namespace costa {

comm_volume communication_volume(assigned_grid2D& g_init,
                                 assigned_grid2D& g_final,
                                 char trans) {
    // transpose
    if (trans != 'N') g_init.transpose();
    grid_cover g_cover(g_init.grid(), g_final.grid());

    int n_blocks_row = g_init.grid().n_rows;
    int n_blocks_col = g_init.grid().n_cols;

    std::unordered_map<edge_t, size_t> weights;

    for (int i = 0; i < n_blocks_row; ++i) {
        for (int j = 0; j < n_blocks_col; ++j) {
            auto rank_to_comm_vol = utils::rank_to_comm_vol_for_block(
                g_init, block_coordinates{i, j}, g_cover, g_final);
            int rank = g_init.owner(i, j);

            for (const auto& comm_vol : rank_to_comm_vol) {
                int target_rank = comm_vol.first;
                int weight = comm_vol.second;

                int smaller_rank = std::min(rank, target_rank);
                int larger_rank = std::max(rank, target_rank);

                edge_t edge_between_ranks{smaller_rank, larger_rank};
                // edge_t edge_between_ranks{rank, target_rank};
                weights[edge_between_ranks] += weight;
            }
        }
    }

    // transpose back
    if (trans != 'N') g_init.transpose();
    return comm_volume(std::move(weights));
}

template <typename T>
void exchange_async(communication_data<T>& send_data, 
                    communication_data<T>& recv_data,
                    MPI_Comm comm) {
    memory::threads_workspace<T> workspace(64);

    PE(transform_irecv);
    MPI_Request* recv_reqs;
    // protect from empty data
    if (recv_data.n_packed_messages > 0) {
        recv_reqs = new MPI_Request[recv_data.n_packed_messages];
    }
    int request_idx = 0;
    // initiate all receives
    for (unsigned i = 0u; i < recv_data.n_ranks; ++i) {
        if (recv_data.counts[i] > 0) {
            MPI_Irecv(recv_data.data() + recv_data.dspls[i],
                      recv_data.counts[i],
                      mpi_type_wrapper<T>::type(),
                      i, 0, comm,
                      &recv_reqs[request_idx]);
            ++request_idx;
        }
    }
    // MPI_Startall(recv_data.n_packed_messages, recv_reqs);
    PL();

    PE(transform_packing);
    // copy blocks to temporary send buffers
    send_data.copy_to_buffer(workspace);
    PL();

    PE(transform_isend);
    MPI_Request* send_reqs;
    if (send_data.n_packed_messages > 0) {
        send_reqs = new MPI_Request[send_data.n_packed_messages];
    }
    request_idx = 0;
    // initiate all sends
    for (unsigned i = 0u; i < send_data.n_ranks; ++i) {
        if (send_data.counts[i] > 0) {
            MPI_Isend(send_data.data() + send_data.dspls[i], 
                      send_data.counts[i],
                      mpi_type_wrapper<T>::type(),
                      i, 0, comm,
                      &send_reqs[request_idx]);
            ++request_idx;
        }
    }
    // MPI_Startall(send_data.n_packed_messages, send_reqs);
    PL();

    PE(transform_localblocks);
    // copy local data (that are on the same rank in both initial and final layout)
    // this is independent of MPI and can be executed in parallel
    copy_local_blocks(send_data.local_messages,
                      recv_data.local_messages,
                      workspace);
    PL();

    // wait for any package and immediately unpack it
    for (unsigned i = 0u; i < recv_data.n_packed_messages; ++i) {
        int idx;
        PE(transform_waitany);
        MPI_Waitany(recv_data.n_packed_messages,
                    recv_reqs,
                    &idx,
                    MPI_STATUS_IGNORE);
        PL();
        PE(transform_unpacking);
        // unpack the package that arrived
        recv_data.copy_from_buffer(idx, workspace);
        PL();
    }

    if (recv_data.n_packed_messages > 0) {
        delete[] recv_reqs;
    }

    PE(transform_waitall);
    // finish up the send requests since all the receive requests are finished
    if (send_data.n_packed_messages > 0) {
        MPI_Waitall(send_data.n_packed_messages, send_reqs, MPI_STATUSES_IGNORE);
        delete[] send_reqs;
    }
    PL();
}

template <typename T>
void exchange(communication_data<T>& send_data, communication_data<T>& recv_data, MPI_Comm comm) {
    memory::threads_workspace<T> workspace(64);

    // copy blocks to temporary send buffers
    PE(transformation_copy2buffer);
    send_data.copy_to_buffer(workspace);

    PE(transformation_all2all);
    MPI_Alltoallv(send_data.data(),
                  send_data.counts.data(),
                  send_data.dspls.data(),
                  mpi_type_wrapper<T>::type(),
                  recv_data.data(),
                  recv_data.counts.data(),
                  recv_data.dspls.data(),
                  mpi_type_wrapper<T>::type(),
                  comm);
    PL();

    PE(transformation_copyfrombuffer);
    recv_data.copy_from_buffer(workspace);
    PL();

    // copy local data (that are on the same rank in both initial and final layout)
    // this is independent of MPI and can be executed in parallel
    PE(transformation_localblocks);
    // copy local data (that are on the same rank in both initial and final layout)
    // this is independent of MPI and can be executed in parallel
    copy_local_blocks(send_data.local_messages,
                      recv_data.local_messages,
                      workspace);
    PL();
}

template <typename T>
void transform(grid_layout<T> &initial_layout,
               grid_layout<T> &final_layout,
               MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    // transposing depends also on the ordering
    // i.e. whether the initial and the final layouts
    // have blocks in row- or col-major ordering
    bool transpose = utils::if_should_transpose(
                                         initial_layout.ordering, 
                                         final_layout.ordering, 
                                         'N');
    // transpose if needed to make layouts compatible 
    // (i.e. same dimensions) before comparing the grids
    if (transpose) initial_layout.transpose();

    // no transpose, no conjugate (false, false)
    costa::communication_data<T> send_data;
    costa::communication_data<T> recv_data;
    costa::memory::memory_buffer<message<T>> messages_buffer;
#pragma omp parallel sections
    {
#pragma omp section
    send_data = 
        utils::prepare_to_send(messages_buffer, initial_layout, final_layout, rank,
                               T{1}, T{0}, transpose, false);
#pragma omp section
    recv_data = 
        utils::prepare_to_recv(messages_buffer, final_layout, initial_layout, rank,
                               T{1}, T{0}, transpose, false);
    }
    // undo the transpose
    if (transpose) initial_layout.transpose();

    // perform the communication
    // exchange_async(send_data, recv_data, comm);
    exchange(send_data, recv_data, comm);
}

template <typename T>
void transform(grid_layout<T> &initial_layout,
               grid_layout<T> &final_layout,
               const char trans,
               const T alpha, const T beta,
               MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    char flag = std::toupper(trans);
    assert(flag == 'N' || flag == 'T' || flag == 'C');

    // transposing depends also on the ordering
    // i.e. whether the initial and the final layouts
    // have blocks in row- or col-major ordering
    bool transpose = utils::if_should_transpose(
                                         initial_layout.ordering, 
                                         final_layout.ordering, 
                                         flag);
    bool conjugate = flag == 'C';

    // transpose if needed to make layouts compatible 
    // (i.e. same dimensions) before comparing the grids
    if (transpose) initial_layout.transpose();

    costa::communication_data<T> send_data;
    costa::communication_data<T> recv_data;
    costa::memory::memory_buffer<message<T>> messages_buffer;
#pragma omp parallel sections
    {
#pragma omp section
    send_data =
        utils::prepare_to_send(messages_buffer, 
                               initial_layout, final_layout, rank, 
                               alpha, beta, transpose, conjugate);

#pragma omp section
    recv_data =
        utils::prepare_to_recv(messages_buffer,
                               final_layout, initial_layout, rank, 
                               alpha, beta, transpose, conjugate);
    }

    // undo the transpose
    if (transpose) initial_layout.transpose();

    // perform the communication
    // exchange_async(send_data, recv_data, comm);
    exchange(send_data, recv_data, comm);
}

template <typename T>
void transform(std::vector<layout_ref<T>>& from,
               std::vector<layout_ref<T>>& to,
               MPI_Comm comm) {
    assert(from.size() == to.size());

    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<T> alpha(from.size(), T{1});
    std::vector<T> beta(from.size(), T{0});

    // transpose
    bool transpose[from.size()];
    std::fill_n(transpose, from.size(), false);
    // conjugate
    bool conjugate[from.size()];
    std::fill_n(conjugate, from.size(), false);

    costa::communication_data<T> send_data;
    costa::communication_data<T> recv_data;
    costa::memory::memory_buffer<message<T>> messages_buffer;
#pragma omp parallel sections
    {
#pragma omp section
    send_data = utils::prepare_to_send(messages_buffer,
                                       from, to, rank, 
                                       &alpha[0], &beta[0], 
                                       transpose, conjugate);
#pragma omp section
    recv_data = utils::prepare_to_recv(messages_buffer,
                                       to, from, rank, 
                                       &alpha[0], &beta[0], 
                                       transpose, conjugate);
    }

    // exchange_async(send_data, recv_data, comm);
    exchange(send_data, recv_data, comm);
}

template <typename T>
void transform(std::vector<layout_ref<T>>& from,
               std::vector<layout_ref<T>>& to,
               const char* trans,
               const T* alpha, const T* beta,
               MPI_Comm comm) {
    // THIS IS THE ONE
    assert(from.size() == to.size());

    int rank;
    MPI_Comm_rank(comm, &rank);

    // transpose
    bool transpose[from.size()];
    // conjugate
    bool conjugate[from.size()];

    // transpose each initial layout as specified by flags
    for (unsigned i = 0u; i < from.size(); ++i) {
        char flag = std::toupper(trans[i]);
        assert(flag == 'N' || flag == 'T' || flag == 'C');

        // transposing depends also on the ordering
        // i.e. whether the initial and the final layouts
        // have blocks in row- or col-major ordering
        transpose[i] = utils::if_should_transpose(from[i].get().ordering,
                                                  to[i].get().ordering,
                                                  flag);

        // check if should conjugate
        conjugate[i] = flag == 'C';

        if (transpose[i]) {
            from[i].get().transpose();
        }
    }

    costa::memory::memory_buffer<message<T>> messages_buffer;

    costa::communication_data<T> send_data;
    costa::communication_data<T> recv_data;
//#pragma omp parallel sections num_threads(2)
    {
// #pragma omp section
        {
            auto start = std::chrono::steady_clock::now();
            send_data = utils::prepare_to_send(messages_buffer,
                    from, to, rank, 
                    alpha, beta, 
                    transpose, conjugate);
            auto end = std::chrono::steady_clock::now();
            auto timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            /*
            if (rank ==0)
                std::cout << "Preparation Time for send [ms] = " << timing << std::endl;
                */
        }
// #pragma omp section
        {
            auto start = std::chrono::steady_clock::now();
            recv_data = utils::prepare_to_recv(messages_buffer,
                    to, from, rank, 
                    alpha, beta, 
                    transpose, conjugate);
            auto end = std::chrono::steady_clock::now();
            auto timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            /*
            if (rank ==0)
                std::cout << "Preparation Time for recv [ms] = " << timing << std::endl;
                */
        }
    }

    // undo the transpose
    for (unsigned i = 0u; i < from.size(); ++i) {
        if (transpose[i]) {
            from[i].get().transpose();
        }
    }

    auto start = std::chrono::steady_clock::now();
    exchange_async(send_data, recv_data, comm);
    // exchange(send_data, recv_data, comm);
    auto end = std::chrono::steady_clock::now();
    auto timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    /*
    if (rank ==0)
        std::cout << "Exchange async [ms] = " << timing << std::endl;
        */
}

// explicit instantiation of transforming a single pair of layouts
template void transform<float>(grid_layout<float> &initial_layout,
                               grid_layout<float> &final_layout,
                               MPI_Comm comm);

template void transform<double>(grid_layout<double> &initial_layout,
                                grid_layout<double> &final_layout,
                                MPI_Comm comm);

template void
transform<std::complex<float>>(grid_layout<std::complex<float>> &initial_layout,
                               grid_layout<std::complex<float>> &final_layout,
                               MPI_Comm comm);

template void transform<std::complex<double>>(
    grid_layout<std::complex<double>> &initial_layout,
    grid_layout<std::complex<double>> &final_layout,
    MPI_Comm comm);

// explicit instantiation of transforming a single pair of layouts
// (with scaling parameters alpha and beta)
template void transform<float>(grid_layout<float> &initial_layout,
                               grid_layout<float> &final_layout,
                               const char trans,
                               const float alpha, const float beta,
                               MPI_Comm comm);

template void transform<double>(grid_layout<double> &initial_layout,
                                grid_layout<double> &final_layout,
                                const char trans,
                                const double alpha, const double beta,
                                MPI_Comm comm);

template void
transform<std::complex<float>>(grid_layout<std::complex<float>> &initial_layout,
                               grid_layout<std::complex<float>> &final_layout,
                               const char trans,
                               const std::complex<float> alpha, const std::complex<float> beta,
                               MPI_Comm comm);

template void transform<std::complex<double>>(
    grid_layout<std::complex<double>> &initial_layout,
    grid_layout<std::complex<double>> &final_layout,
    const char trans,
    const std::complex<double> alpha, const std::complex<double> beta,
    MPI_Comm comm);

// explicit instantiation of transform with vectors
template void transform<float>(std::vector<layout_ref<float>>& initial_layouts,
                               std::vector<layout_ref<float>>& final_layouts,
                               MPI_Comm comm);

template void transform<double>(std::vector<layout_ref<double>>& initial_layouts,
                                std::vector<layout_ref<double>>& final_layouts,
                                MPI_Comm comm);

template void transform<std::complex<float>>(
                               std::vector<layout_ref<std::complex<float>>>& initial_layouts,
                               std::vector<layout_ref<std::complex<float>>>& final_layouts,
                               MPI_Comm comm);

template void transform<std::complex<double>>(
                               std::vector<layout_ref<std::complex<double>>>& initial_layouts,
                               std::vector<layout_ref<std::complex<double>>>& final_layouts,
                               MPI_Comm comm);

// explicit instantiation of transform with vectors
// (with scaling parameters alpha and beta)
template void transform<float>(std::vector<layout_ref<float>>& initial_layouts,
                               std::vector<layout_ref<float>>& final_layouts,
                               const char* trans,
                               const float* alpha, const float* beta,
                               MPI_Comm comm);

template void transform<double>(std::vector<layout_ref<double>>& initial_layouts,
                                std::vector<layout_ref<double>>& final_layouts,
                                const char* trans,
                                const double* alpha, const double* beta,
                                MPI_Comm comm);

template void transform<std::complex<float>>(
                               std::vector<layout_ref<std::complex<float>>>& initial_layouts,
                               std::vector<layout_ref<std::complex<float>>>& final_layouts,
                               const char* trans,
                               const std::complex<float>* alpha, const std::complex<float>* beta,
                               MPI_Comm comm);

template void transform<std::complex<double>>(
                               std::vector<layout_ref<std::complex<double>>>& initial_layouts,
                               std::vector<layout_ref<std::complex<double>>>& final_layouts,
                               const char* trans,
                               const std::complex<double>* alpha, const std::complex<double>* beta,
                               MPI_Comm comm);

} // namespace costa
