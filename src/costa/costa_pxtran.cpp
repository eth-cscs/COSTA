#include <cassert>
#include <mpi.h>

#include <costa/blacs.hpp>
#include <costa/costa_pxtran.hpp>

#include <grid2grid/ranks_reordering.hpp>
#include <grid2grid/transformer.hpp>

namespace costa {
template <typename T>
void pxtran(
           const int m,
           const int n,
           const T alpha,
           const T *a,
           const int ia,
           const int ja,
           const int *desca,
           const T beta,
           T *c, // result
           const int ic,
           const int jc,
           const int *descc) {
    // **********************************
    //           CORNER CASES
    // **********************************
    // edge cases, which are allowed by the standard
    if (m == 0 || n == 0) return;

    // **********************************
    //           MAIN CODE
    // **********************************
    // blas context
    int ctxt = scalapack::get_grid_context(desca, descc);

    // scalapack rank grid decomposition
    int procrows, proccols;
    int myrow, mycol;
    blacs::Cblacs_gridinfo(ctxt, &procrows, &proccols, &myrow, &mycol);

    // get MPI communicator
    MPI_Comm comm = scalapack::get_communicator(ctxt);

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
    int a_subm = n;
    int a_subn = m;

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
        pxtran_params<T> params(
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
                             // alpha, beta
                             alpha, beta,
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
        'T',
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

    /*
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
    */

    // transform A to C
    if (alpha != T{1} || beta != T{0}) {
        // scale while transforming
        grid2grid::transform<T>(scalapack_layout_a, scalapack_layout_c, alpha, beta, comm);
    } else {
        // use a more efficient copy since C can be overwritten
        grid2grid::transform<T>(scalapack_layout_a, scalapack_layout_c, comm);
    }

    /*
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
    */
}

// explicit instantiation for pxtran
template void pxtran<double>(
                            const int m,
                            const int n,
                            const double alpha,
                            const double *a,
                            const int ia,
                            const int ja,
                            const int *desca,
                            const double beta,
                            double *c,
                            const int ic,
                            const int jc,
                            const int *descc);

template void pxtran<float>(
                           const int m,
                           const int n,
                           const float alpha,
                           const float *a,
                           const int ia,
                           const int ja,
                           const int *desca,
                           const float beta,
                           float *c,
                           const int ic,
                           const int jc,
                           const int *descc);

template void pxtran<zdouble_t>(
                               const int m,
                               const int n,
                               const zdouble_t alpha,
                               const zdouble_t *a,
                               const int ia,
                               const int ja,
                               const int *desca,
                               const zdouble_t beta,
                               zdouble_t *c,
                               const int ic,
                               const int jc,
                               const int *descc);

template void pxtran<zfloat_t>(
                              const int m,
                              const int n,
                              const zfloat_t alpha,
                              const zfloat_t *a,
                              const int ia,
                              const int ja,
                              const int *desca,
                              const zfloat_t beta,
                              zfloat_t *c,
                              const int ic,
                              const int jc,
                              const int *descc);
} // namespace costa
