#include <cassert>
#include <mpi.h>

#include <costa/blacs.hpp>
#include <costa/pxgemr2d_params.hpp>
#include <costa/costa_pxgemr2d.hpp>

#include <grid2grid/ranks_reordering.hpp>
#include <grid2grid/transformer.hpp>

namespace costa {
// returns true if subcomm is a subcommunicator of comm
// i.e. checks if intersection(comm, subcomm) == subcomm
bool is_subcommunicator(MPI_Comm comm, MPI_Comm subcomm) {
    // get the groups from the given communicators
    MPI_Group group;
    MPI_Group subgroup;
    MPI_Comm_group(comm, &group);
    MPI_Comm_group(subcomm, &subgroup);

    // get the intersection of the two groups
    MPI_Group intersection;
    MPI_Group_intersection(group, subgroup, &intersection);

    // check if intersection == subcomm (meaning that subcomm is a subset of comm)
    int comp;
    MPI_Group_compare(intersection, subgroup, &comp);

    return comp == MPI_SIMILAR;
}

template <typename T>
void pxgemr2d(
           const int m,
           const int n,
           const T *a,
           const int ia,
           const int ja,
           const int *desca,
           T *c,
           const int ic,
           const int jc,
           const int *descc,
           int ctxt) {
    // **********************************
    //           CORNER CASES
    // **********************************
    // edge cases, which are allowed by the standard
    if (m == 0 || n == 0) return;

    // **********************************
    //           MAIN CODE
    // **********************************
    // blas context
    int ctxt_a = scalapack::get_grid_context(desca);
    int ctxt_c = scalapack::get_grid_context(descc);

    // scalapack rank grid decomposition
    int procrows, proccols;
    int myrow, mycol;
    blacs::Cblacs_gridinfo(ctxt, &procrows, &proccols, &myrow, &mycol);

    // get MPI communicators
    MPI_Comm comm_a = scalapack::get_communicator(ctxt_a);
    MPI_Comm comm_c = scalapack::get_communicator(ctxt_c);
    MPI_Comm comm = scalapack::get_communicator(ctxt);
    // check if comm is at least the union of comm_a and comm_c
    assert(is_subcommunicator(comm, comm_a));
    assert(is_subcommunicator(comm, comm_c));

    // communicator size and rank
    int rank, P;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &rank);

    // block sizes
    scalapack::block_size b_dim_a(desca);
    scalapack::block_size b_dim_c(descc);

    // global matrix sizes
    scalapack::global_matrix_size mat_dim_a(desca);
    scalapack::global_matrix_size mat_dim_c(descc);

    // sumatrix size to multiply
    int a_subm = m;
    int a_subn = n;

    int c_subm = m;
    int c_subn = n;

    // rank sources (rank coordinates that own first row and column of a matrix)
    scalapack::rank_src rank_src_a(desca);
    scalapack::rank_src rank_src_c(descc);

    // leading dimensions
    int lld_a = scalapack::leading_dimension(desca);
    int lld_c = scalapack::leading_dimension(descc);

    // check whether rank grid is row-major or col-major
    auto ordering = scalapack::rank_ordering(ctxt, P);
    char grid_order =
        ordering == grid2grid::scalapack::ordering::column_major ? 'C' : 'R';

#ifdef DEBUG
    if (rank == 0) {
        pxgemr2d_params<T> params(
                             // global dimensions
                             mat_dim_a.rows, mat_dim_a.cols,
                             mat_dim_c.rows, mat_dim_c.cols,
                             // block dimensions
                             b_dim_a.rows, b_dim_a.cols,
                             b_dim_c.rows, b_dim_c.cols,
                             // submatrix start
                             ia, ja,
                             ic, jc,
                             // problem size
                             m, n,
                             // leading dimensinons
                             lld_a, lld_c,
                             // processor grid
                             procrows, proccols,
                             // processor grid ordering
                             grid_order,
                             // ranks containing first rows
                             rank_src_a.row_src, rank_src_a.col_src,
                             rank_src_c.row_src, rank_src_c.col_src
                         );
        std::cout << params << std::endl;
    }
    MPI_Barrier(comm);
#endif

#ifdef DEBUG
    if (rank == 0) {
        std::cout << strategy << std::endl;
        std::cout << "============================================" << std::endl;
    }
    MPI_Barrier(comm);
#endif

    // get abstract layout descriptions for ScaLAPACK layout
    auto scalapack_layout_a = grid2grid::get_scalapack_grid<T>(
        lld_a,
        {mat_dim_a.rows, mat_dim_a.cols},
        {ia, ja},
        {a_subm, a_subn},
        {b_dim_a.rows, b_dim_a.cols},
        {procrows, proccols},
        ordering,
        'N',
        {rank_src_a.row_src, rank_src_a.col_src},
        a,
        rank);

    auto scalapack_layout_c = grid2grid::get_scalapack_grid<T>(
        lld_c,
        {mat_dim_c.rows, mat_dim_c.cols},
        {ic, jc},
        {c_subm, c_subn},
        {b_dim_c.rows, b_dim_c.cols},
        {procrows, proccols},
        ordering,
        'N',
        {rank_src_c.row_src, rank_src_c.col_src},
        c,
        rank);

    // total communication volume for transformation of layouts
    auto comm_vol = grid2grid::communication_volume(scalapack_layout_a.grid, scalapack_layout_c.grid);

    // compute the optimal rank reordering that minimizes the communication volume
    bool reordered = false;
    std::vector<int> rank_permutation = grid2grid::optimal_reordering(comm_vol, P, reordered);

    // create reordered communicator, which has same ranks
    // but relabelled as given by the rank_permutation
    // (to avoid the communication during layout transformation)
    MPI_Comm reordered_comm = comm;
    if (reordered) {
        MPI_Comm_split(comm, 0, rank_permutation[rank], &reordered_comm);
    }

#ifdef DEBUG
    if (rank == 0) {
        std::cout << "Optimal rank relabeling:" << std::endl;
        for (int i = 0; i < P; ++i) {
            std::cout << i << "->" << rank_permutation[i] << std::endl;
        }
    }
#endif

    // transform A to C
    grid2grid::transform<T>(scalapack_layout_a, scalapack_layout_c, comm);

#ifdef DEBUG
    if (rank == 0) {
        auto reordered_vol = grid2grid::communication_volume(scalapack_layout_a.grid, scalapack_layout_c.grid);

        // std::cout << "Detailed comm volume: " << comm_vol << std::endl;
        // std::cout << "Detailed comm volume reordered: " << reordered_vol << std::endl;

        auto comm_vol_total = comm_vol.total_volume();
        auto reordered_vol_total = reordered_vol.total_volume();
        std::cout << "Initial comm volume = " << comm_vol_total << std::endl;
        std::cout << "Reduced comm volume = " << reordered_vol_total << std::endl;
        auto diff = (long long) comm_vol_total - (long long) reordered_vol_total;
        std::cout << "Comm volume reduction [%] = " << 100.0 * diff / comm_vol_total << std::endl;

    }
#endif
    if (reordered) {
        MPI_Comm_free(&reordered_comm);
    }
}

// explicit instantiation for pxtran
template void pxgemr2d<double>(
                            const int m,
                            const int n,
                            const double *a,
                            const int ia,
                            const int ja,
                            const int *desca,
                            double *c,
                            const int ic,
                            const int jc,
                            const int *descc,
                            int ctxt);

template void pxgemr2d<float>(
                           const int m,
                           const int n,
                           const float *a,
                           const int ia,
                           const int ja,
                           const int *desca,
                           float *c,
                           const int ic,
                           const int jc,
                           const int *descc,
                           int ctxt);

template void pxgemr2d<zdouble_t>(
                               const int m,
                               const int n,
                               const zdouble_t *a,
                               const int ia,
                               const int ja,
                               const int *desca,
                               zdouble_t *c,
                               const int ic,
                               const int jc,
                               const int *descc,
                               int ctxt);

template void pxgemr2d<zfloat_t>(
                              const int m,
                              const int n,
                              const zfloat_t *a,
                              const int ia,
                              const int ja,
                              const int *desca,
                              zfloat_t *c,
                              const int ic,
                              const int jc,
                              const int *descc,
                              int ctxt);
} // namespace costa
