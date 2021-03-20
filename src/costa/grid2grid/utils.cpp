#include <costa/grid2grid/utils.hpp>

bool costa::utils::if_should_transpose(const char src_ordering, 
                                       const char dest_ordering,
                                       const char trans) {
    assert(src_ordering == 'R' || src_ordering == 'C');
    assert(dest_ordering == 'R' || dest_ordering == 'C');
    assert(trans == 'N' || trans == 'T' || trans == 'C');
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

