#include <costa/grid2grid/utils.hpp>

#include <umpire/ResourceManager.hpp>
#include <umpire/strategy/QuickPool.hpp>
#include <umpire/strategy/ThreadSafeAllocator.hpp>
#include <umpire/TypedAllocator.hpp>

template <typename T>
using messages_vector = std::vector<T, umpire::TypedAllocator<T>>;

template <typename T>
messages_vector<T> get_messages_vector() {
    auto& rm = umpire::ResourceManager::getInstance();
    auto alloc = rm.getAllocator("HOST");

    // _sphinx_tag_tut_typed_alloc_start
    umpire::TypedAllocator<T> allocator{alloc};
    messages_vector<T> messages{allocator};

    return messages;
}

bool costa::utils::if_should_transpose(const char src_ordering, 
                                       const char dest_ordering,
                                       const char trans) {
    assert(src_ordering == 'R' || src_ordering == 'C');
    assert(dest_ordering == 'R' || dest_ordering == 'C');
    assert(trans == 'N' || trans == 'T' || trans == 'C');
    return trans != 'N';
    /*
    // BE CAREFUL: transpose and different src and dest orderings might cancel out
    // ===========
    // Row-major + Transpose + Row-major = Transpose (Row-major)
    // Col-major + Transpose + Col-major = Transpose (Col-major)
    // Row-major + Transpose + Col-major = Copy(Row-major) // cancels out
    // Col-major + Transpose + Row-major = Copy(Col-major) // cancels out
    //
    // Row-major + NoTranspose + Row-major = Copy(Row-major)
    // Col-major + NoTranspose + Col-major = Copy(Col-major)
    // Row-major + NoTranspose + Col-major = Transpose(Row-major)
    // Col-major + NoTranspose + Row-major = Transpose(Col-major)
    bool transpose = trans != 'N';

    bool should_transpose = (transpose && src_ordering == dest_ordering)
                             ||
                            (!transpose && src_ordering != dest_ordering);
    return should_transpose;
    */
}
std::vector<std::vector<int>> costa::topology_cost(MPI_Comm comm) {
    int P;
    MPI_Comm_size(comm, &P);
    // default cost factor 2
    // ranks sharing the same node will have the cost factor = 1
    std::vector<std::vector<int>> cost(P, std::vector<int>(P, 1));

    // split global comm into local_comms where each local comm contains
    // all ranks which share the same node
    int rank, local_rank, local_size;
    MPI_Comm local_comm;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank,  MPI_INFO_NULL, &local_comm);

    // size of the local subcomm
    int small_P;
    MPI_Comm_size(local_comm, &small_P);

    // global group
    MPI_Group group;
    MPI_Comm_group(comm, &group);

    // local group
    MPI_Group local_group;
    MPI_Comm_group(local_comm, &local_group);

    // we want to translate all local ranks (sharing the current node) 
    // to their global identifiers
    std::vector<int> local_ranks(small_P);
    for (int i = 0; i < local_ranks.size(); ++i) {
        local_ranks[i] = i;
    }

    std::vector<int> local2global(small_P);

    // translate local-group ranks to global ranks
    MPI_Group_translate_ranks(local_group, small_P, &local_ranks[0], 
                              group, &local2global[0]);

    // each node is labeled by the smallest rank residing on that node
    int node_idx = *std::min_element(local2global.begin(), local2global.end());

    std::vector<int> rank2node(P);

    // each rank reports its local node
    MPI_Allgather(&node_idx, 1, MPI_INT, &rank2node[0], 1, MPI_INT, comm);

    for (int i = 0; i < P; ++i) {
        for (int j = i; j < P; ++j) {
            // if ranks share the same node, halve the cost
            if (rank2node[i] == rank2node[j]) {
                cost[i][j] = 2;
                cost[j][i] = 2;
            }
        }
    }

    return cost;
}

std::unordered_map<int, int> costa::utils::rank_to_comm_vol_for_block(
        const assigned_grid2D& g_init,
        const block_coordinates &b_coord,
        grid_cover &g_cover,
        const assigned_grid2D& g_final) {
    // std::cout << "decomposing block " << b << std::endl;
    block_cover b_cover = g_cover.decompose_block(b_coord);

    int row_first = b_cover.rows_cover.start_index;
    int row_last = b_cover.rows_cover.end_index;

    int col_first = b_cover.cols_cover.start_index;
    int col_last = b_cover.cols_cover.end_index;

    auto rows_interval = g_init.rows_interval(b_coord.row);
    auto cols_interval = g_init.cols_interval(b_coord.col);

    std::unordered_map<int, int> comm_vol;

    int row_start = rows_interval.start;
    // use start of the interval to get the rank and the end of the interval
    // to get the block which has to be sent
    // skip the last element
    for (int i = row_first; i < row_last; ++i) {
        int row_end = std::min(g_final.grid().rows_split[i + 1], rows_interval.end);

        int col_start = cols_interval.start;
        for (int j = col_first; j < col_last; ++j) {
            // use i, j to find out the rank
            int rank = g_final.owner(i, j);
            // std::cout << "owner of block " << i << ", " << j << " is " <<
            // rank << std::endl;

            // use i+1 and j+1 to find out the block
            int col_end =
                std::min(g_final.grid().cols_split[j + 1], cols_interval.end);

            int size = (row_end - row_start) * (col_end - col_start);
            // if non empty, add this block
            if (size > 0) {
                comm_vol[rank] += size;
            }

            col_start = col_end;
        }
        row_start = row_end;
    }
    return comm_vol;
}

template <typename T>
void decompose_block(std::vector<costa::message<T>>& messages,
                     const costa::block<T> &b,
                     costa::grid_cover &g_cover,
                     const costa::assigned_grid2D &g,
                     const char final_ordering,
                     const T alpha, const T beta,
                     bool transpose, bool conjugate
                     ) {
    // std::cout << "decomposing block " << b << std::endl;
    costa::block_cover b_cover = g_cover.decompose_block(b);

    int row_first = b_cover.rows_cover.start_index;
    int row_last = b_cover.rows_cover.end_index;

    int col_first = b_cover.cols_cover.start_index;
    int col_last = b_cover.cols_cover.end_index;

    int n_blocks = (col_last - col_first) * (row_last - row_first);

    messages.reserve(n_blocks);

    // use start of the interval to get the rank and the end of the interval
    // to get the block which has to be sent
    // skip the last element
    int col_start = b.cols_interval.start;
    for (int j = col_first; j < col_last; ++j) {
        int row_start = b.rows_interval.start;
        // use i+1 and j+1 to find out the block
        int col_end =
            std::min(g.grid().cols_split[j + 1], b.cols_interval.end);
        for (int i = row_first; i < row_last; ++i) {
            int row_end = std::min(g.grid().rows_split[i + 1], b.rows_interval.end);
            // use i, j to find out the rank
            int rank = g.owner(i, j);
            // std::cout << "owner of block " << i << ", " << j << " is " <<
            // rank << std::endl;

            // get pointer to this block of data based on the internal local
            // layout
            costa::block<T> subblock =
                b.subblock({row_start, row_end}, {col_start, col_end});

            assert(subblock.non_empty());
            // if non empty, add this block
            if (subblock.non_empty()) {
                // std::cout << "for rank " << rank << ", adding subblock: " <<
                // subblock << std::endl; std::cout << "owner of " << subblock
                // << " is " << rank << std::endl;
                messages.push_back({subblock, rank,
                                    final_ordering,
                                    alpha, beta,
                                    transpose, conjugate});
            }
            row_start = row_end;
        }
        col_start = col_end;
    }
}

template <typename T>
std::vector<T> k_way_merge(std::vector<std::vector<T>>& arr) {
    if (arr.size() == 0) return {};
    if (arr.size() == 1) return arr[0];
    // A pair of pairs, first element is going to
    // store value, second element index of array
    // and third element index in the array.
    using ppi = std::pair<T, std::pair<std::size_t, std::size_t>>;

    std::vector<T> output;
    output.reserve(arr.size() * arr[0].size());

    // Create a min heap with k heap nodes. Every
    // heap node has first element of an array
    std::priority_queue<ppi, std::vector<ppi>, std::greater<ppi> > pq;

    for (std::size_t i = 0; i < arr.size(); i++) {
        pq.push({ arr[i][0], { i, 0 } });
    }

    // Now one by one get the minimum element
    // from min heap and replace it with next
    // element of its array
    while (!pq.empty()) {
        ppi curr = pq.top();
        pq.pop();

        // i ==> Array Number
        // j ==> Index in the array number
        auto i = curr.second.first;
        auto j = curr.second.second;

        output.push_back(curr.first);

        // The next element belongs to same array as
        // current.
        if (j + 1 < arr[i].size())
            pq.push({ arr[i][j + 1], { i, j + 1 } });
    }

    return output;
}

template <typename T>
void merge_sort(std::vector<T>& v, std::size_t left, std::size_t right) {
    if (left < right) {
        if (right-left >= 128) {
            std::size_t mid = (left + right) / 2;
#pragma omp taskgroup
            {
#pragma omp task shared(v) untied if(right - left >= (1<<14))
                merge_sort(v, left, mid);
#pragma omp task shared(v) untied if(right - left >= (1<<14))
                merge_sort(v, mid+1, right);
#pragma omp taskyield
            }
            std::inplace_merge(v.begin() + left, v.begin() + mid + 1, v.begin() + right + 1);
        } else{
            std::sort(v.begin() + left, v.begin() + right + 1);
        }
    }
}

template <typename T>
void parallel_sort(std::vector<T>& v) {
    std::size_t n_threads = (std::size_t) std::max(1, omp_get_max_threads());
    n_threads = std::min(n_threads, v.size());

    #pragma omp parallel num_threads(n_threads)
    #pragma omp single
    merge_sort(v, 0, v.size()-1);
}

template <typename T>
std::vector<costa::message<T>> decompose_blocks(
                                         costa::memory::memory_buffer<costa::message<T>>& messages_buffer,
                                         costa::grid_layout<T> &init_layout,
                                         costa::grid_layout<T> &final_layout,
                                         const T alpha, const T beta,
                                         bool transpose,
                                         bool conjugate,
                                         int tag = 0) {
    auto start = std::chrono::steady_clock::now();
    PE(transform_decompose);
    costa::grid_cover g_overlap(init_layout.grid.grid(), final_layout.grid.grid());

    int half_threads = std::max(1, omp_get_max_threads());
    int n_threads = std::min(half_threads, init_layout.blocks.num_blocks());

    // messages_vector<T> messages = get_messages_vector<T>();

    if (init_layout.blocks.num_blocks() == 1) {
        std::vector<costa::message<T>> result;
        // result.reserve(init_layout.blocks.num_blocks());
        for (int i = 0; i < init_layout.blocks.num_blocks(); ++i) {
            auto blk = init_layout.blocks.get_block(i);
            blk.tag = tag;
            assert(blk.non_empty());
            decompose_block(result,
                            blk, g_overlap, 
                            final_layout.grid,
                            final_layout.ordering,
                            alpha, beta, transpose, conjugate);
            // result.insert(result.end(), decomposed_blocks.begin(), decomposed_blocks.end());
        }

        if (result.size() > 1000)
            parallel_sort(result);
        else 
#pragma omp single
            std::sort(result.begin(), result.end());
        return result;
    }

    std::vector<std::vector<costa::message<T>>> thread_messages(n_threads);

#pragma omp parallel num_threads(n_threads) shared(g_overlap, final_layout, alpha, beta, transpose, conjugate)
    {
    int thread_id = omp_get_thread_num();
    thread_messages[thread_id].reserve(init_layout.blocks.num_blocks()/n_threads);
#pragma omp for
    for (int i = 0; i < init_layout.blocks.num_blocks(); ++i) {
        int thread_id = omp_get_thread_num();
        // std::cout << "decomposing block " << i << " out of " <<
        // init_layout.blocks.num_blocks() << std::endl;
        auto blk = init_layout.blocks.get_block(i);
        blk.tag = tag;
        assert(blk.non_empty());
        decompose_block(thread_messages[thread_id],
                        blk, g_overlap, 
                        final_layout.grid,
                        final_layout.ordering,
                        alpha, beta, transpose, conjugate);
    }

    std::sort(thread_messages[thread_id].begin(), thread_messages[thread_id].end());
    }

    auto result = k_way_merge(thread_messages);

    PL();
    auto end = std::chrono::steady_clock::now();
    auto timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //std::cout << "Decompose blocks [ms] = " << timing << std::endl;

    return result;
}


template <typename T>
costa::communication_data<T> costa::utils::prepare_to_send(
                                      costa::memory::memory_buffer<costa::message<T>>& messages_buffer,
                                      costa::grid_layout<T> &init_layout,
                                      costa::grid_layout<T> &final_layout,
                                      int rank,
                                      const T alpha, const T beta,
                                      bool transpose, bool conjugate) {
    // in case ranks were reordered to minimize the communication
    // this might not be the identity function
    // if (rank == 0) {
    //     std::cout << "prepare to send: changing rank to " << init_layout.reordered_rank(rank) << std::endl;
    // }
    // rank = init_layout.reordered_rank(rank);
    std::vector<costa::message<T>> messages =
        decompose_blocks(messages_buffer, init_layout, final_layout, 
                         alpha, beta, transpose, conjugate);

    return costa::communication_data<T>(messages, rank, std::max(final_layout.num_ranks(), init_layout.num_ranks()));
}

template <typename T>
costa::communication_data<T> costa::utils::prepare_to_send(
                                      costa::memory::memory_buffer<costa::message<T>>& messages_buffer,
                                      std::vector<layout_ref<T>>& from,
                                      std::vector<layout_ref<T>>& to,
                                      int rank,
                                      const T* alpha, const T* beta,
                                      bool* transpose,
                                      bool* conjugate) {
    std::vector<std::vector<costa::message<T>>> messages(from.size());
    int n_ranks = 0;

    for (unsigned i = 0u; i < from.size(); ++i) {
        auto& init_layout = from[i].get();
        auto& final_layout = to[i].get();

        messages[i] = decompose_blocks(messages_buffer,
                                       init_layout, final_layout, 
                                       alpha[i], beta[i], 
                                       transpose[i], conjugate[i], 
                                       i);
        n_ranks = std::max(n_ranks, std::max(final_layout.num_ranks(), init_layout.num_ranks()));
    }
    if (from.size() == 0) {
        return costa::communication_data<T>(messages[0], rank, n_ranks);
    } else {
        auto sorted_messages = k_way_merge(messages);
        return costa::communication_data<T>(sorted_messages, rank, n_ranks);
    }
}
template <typename T> 
costa::communication_data<T> costa::utils::prepare_to_recv(
                                      costa::memory::memory_buffer<costa::message<T>>& messages_buffer,
                                      costa::grid_layout<T> &final_layout,
                                      costa::grid_layout<T> &init_layout,
                                      int rank,
                                      const T alpha, const T beta,
                                      const bool transpose, const bool conjugate) {
    std::vector<costa::message<T>> messages =
        decompose_blocks(messages_buffer, final_layout, init_layout, 
                         alpha, beta, transpose, conjugate);

    return costa::communication_data<T>(messages, rank, std::max(init_layout.num_ranks(), final_layout.num_ranks()));
}

template <typename T>
costa::communication_data<T> costa::utils::prepare_to_recv(costa::memory::memory_buffer<costa::message<T>>& messages_buffer,
                                      std::vector<layout_ref<T>>& to,
                                      std::vector<layout_ref<T>>& from,
                                      int rank,
                                      const T* alpha, const T* beta,
                                      bool* transpose,
                                      bool* conjugate) {
    std::vector<std::vector<costa::message<T>>> messages(from.size());
    int n_ranks = 0;

    for (unsigned i = 0u; i < from.size(); ++i) {
        auto& init_layout = from[i].get();
        auto& final_layout = to[i].get();

        messages[i] = decompose_blocks(
                                       messages_buffer,
                                       final_layout, init_layout, 
                                       alpha[i], beta[i], 
                                       transpose[i], conjugate[i], 
                                       i);
        n_ranks = std::max(n_ranks, std::max(init_layout.num_ranks(), final_layout.num_ranks()));
    }
    if (from.size() == 0) {
        return costa::communication_data<T>(messages[0], rank, n_ranks);
    } else {
        auto sorted_messages = k_way_merge(messages);
        return costa::communication_data<T>(sorted_messages, rank, n_ranks);
    }
}

// template instantiation for prepare_to_send<T>
template costa::communication_data<float> costa::utils::prepare_to_send<float>(
                                      costa::memory::memory_buffer<message<float>>& messages_buffer,
                                      costa::grid_layout<float>& from,
                                      costa::grid_layout<float>& to,
                                      int rank,
                                      const float alpha, const float beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<double> costa::utils::prepare_to_send<double>(
                                      costa::memory::memory_buffer<message<double>>& messages_buffer,
                                      costa::grid_layout<double>& from,
                                      costa::grid_layout<double>& to,
                                      int rank,
                                      const double alpha, const double beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<std::complex<float>> costa::utils::prepare_to_send<std::complex<float>>(
                                      costa::memory::memory_buffer<message<std::complex<float>>>& messages_buffer,
                                      costa::grid_layout<std::complex<float>>& from,
                                      costa::grid_layout<std::complex<float>>& to,
                                      int rank,
                                      const std::complex<float> alpha, const std::complex<float> beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<std::complex<double>> costa::utils::prepare_to_send<std::complex<double>>(
                                      costa::memory::memory_buffer<message<std::complex<double>>>& messages_buffer,
                                      costa::grid_layout<std::complex<double>>& from,
                                      costa::grid_layout<std::complex<double>>& to,
                                      int rank,
                                      const std::complex<double> alpha, const std::complex<double> beta,
                                      bool transpose, bool conjugate);

// template instantiation for the vector version of prepare_to_send<T>
template costa::communication_data<float> costa::utils::prepare_to_send<float>(
                                      costa::memory::memory_buffer<message<float>>& messages_buffer,
                                      std::vector<layout_ref<float>>& from,
                                      std::vector<layout_ref<float>>& to,
                                      int rank,
                                      const float* alpha, const float* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<double> costa::utils::prepare_to_send<double>(
                                      costa::memory::memory_buffer<message<double>>& messages_buffer,
                                      std::vector<layout_ref<double>>& from,
                                      std::vector<layout_ref<double>>& to,
                                      int rank,
                                      const double* alpha, const double* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<std::complex<float>> costa::utils::prepare_to_send<std::complex<float>>(
                                      costa::memory::memory_buffer<message<std::complex<float>>>& messages_buffer,
                                      std::vector<layout_ref<std::complex<float>>>& from,
                                      std::vector<layout_ref<std::complex<float>>>& to,
                                      int rank,
                                      const std::complex<float>* alpha, const std::complex<float>* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<std::complex<double>> costa::utils::prepare_to_send<std::complex<double>>(
                                      costa::memory::memory_buffer<message<std::complex<double>>>& messages_buffer,
                                      std::vector<layout_ref<std::complex<double>>>& from,
                                      std::vector<layout_ref<std::complex<double>>>& to,
                                      int rank,
                                      const std::complex<double>* alpha, const std::complex<double>* beta,
                                      bool* transpose,
                                      bool* conjugate);

// template instantiation for prepare_to_recv<T>
template costa::communication_data<float> costa::utils::prepare_to_recv<float>(
                                      costa::memory::memory_buffer<message<float>>& messages_buffer,
                                      costa::grid_layout<float>& to,
                                      costa::grid_layout<float>& from,
                                      int rank,
                                      const float alpha, const float beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<double> costa::utils::prepare_to_recv<double>(
                                      costa::memory::memory_buffer<message<double>>& messages_buffer,
                                      costa::grid_layout<double>& to,
                                      costa::grid_layout<double>& from,
                                      int rank,
                                      const double alpha, const double beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<std::complex<float>> costa::utils::prepare_to_recv<std::complex<float>>(
                                      costa::memory::memory_buffer<message<std::complex<float>>>& messages_buffer,
                                      costa::grid_layout<std::complex<float>>& to,
                                      costa::grid_layout<std::complex<float>>& from,
                                      int rank,
                                      const std::complex<float> alpha, const std::complex<float> beta,
                                      bool transpose, bool conjugate);
template costa::communication_data<std::complex<double>> costa::utils::prepare_to_recv<std::complex<double>>(
                                      costa::memory::memory_buffer<message<std::complex<double>>>& messages_buffer,
                                      costa::grid_layout<std::complex<double>>& to,
                                      costa::grid_layout<std::complex<double>>& from,
                                      int rank,
                                      const std::complex<double> alpha, const std::complex<double> beta,
                                      bool transpose, bool conjugate);

// template instantiation for the vector version of prepare_to_recv<T>
template costa::communication_data<float> costa::utils::prepare_to_recv<float>(
                                      costa::memory::memory_buffer<message<float>>& messages_buffer,
                                      std::vector<layout_ref<float>>& to,
                                      std::vector<layout_ref<float>>& from,
                                      int rank,
                                      const float* alpha, const float* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<double> costa::utils::prepare_to_recv<double>(
                                      costa::memory::memory_buffer<message<double>>& messages_buffer,
                                      std::vector<layout_ref<double>>& to,
                                      std::vector<layout_ref<double>>& from,
                                      int rank,
                                      const double* alpha, const double* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<std::complex<float>> costa::utils::prepare_to_recv<std::complex<float>>(
                                      costa::memory::memory_buffer<message<std::complex<float>>>& messages_buffer,
                                      std::vector<layout_ref<std::complex<float>>>& to,
                                      std::vector<layout_ref<std::complex<float>>>& from,
                                      int rank,
                                      const std::complex<float>* alpha, const std::complex<float>* beta,
                                      bool* transpose,
                                      bool* conjugate);
template costa::communication_data<std::complex<double>> costa::utils::prepare_to_recv<std::complex<double>>(
                                      costa::memory::memory_buffer<message<std::complex<double>>>& messages_buffer,
                                      std::vector<layout_ref<std::complex<double>>>& to,
                                      std::vector<layout_ref<std::complex<double>>>& from,
                                      int rank,
                                      const std::complex<double>* alpha, const std::complex<double>* beta,
                                      bool* transpose,
                                      bool* conjugate);


