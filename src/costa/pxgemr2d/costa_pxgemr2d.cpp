#include <cassert>
#include <mpi.h>

#include <costa/blacs.hpp>
#include <costa/pxgemr2d/pxgemr2d_params.hpp>
#include <costa/pxgemr2d/costa_pxgemr2d.hpp>
#include <costa/grid2grid/ranks_reordering.hpp>

#include <costa/grid2grid/transform.hpp>

#include <costa/grid2grid/profiler.hpp>

namespace costa {
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
           const int ctxt) {
    // clear the profiler
    // empty if compiled without profiler
    PC();

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
    /*
    MPI_Comm comm = scalapack::get_communicator(ctxt);
    // MPI_Comm comm = blacs::Cblacs2sys_handle(ctxt);
    // check if comm is at least the union of comm_a and comm_c
    assert(is_subcommunicator(comm, comm_a));
    assert(is_subcommunicator(comm, comm_c));
    */

    MPI_Comm comm = scalapack::comm_union(comm_a, comm_c);

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
    auto ordering_a = scalapack::rank_ordering(ctxt_a, P);
    char grid_order_a =
        ordering_a == costa::scalapack::ordering::column_major ? 'C' : 'R';
    // check whether rank grid is row-major or col-major
    auto ordering_c = scalapack::rank_ordering(ctxt_c, P);
    char grid_order_c =
        ordering_c == costa::scalapack::ordering::column_major ? 'C' : 'R';

#ifdef DEBUG
    if (rank == 0) {
        // scalapack rank grid decomposition
        int procrows_a, proccols_a;
        int procrows_c, proccols_c;
        int myrow, mycol;
        blacs::Cblacs_gridinfo(ctxt_a, &procrows_a, &proccols_a, &myrow, &mycol);
        blacs::Cblacs_gridinfo(ctxt_c, &procrows_c, &proccols_c, &myrow, &mycol);

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
                             // processor grid of A
                             procrows_a, proccols_a,
                             // processor grid of C
                             procrows_c, proccols_c,
                             // processor grid ordering
                             grid_order_a, grid_order_c,
                             // ranks containing first rows
                             rank_src_a.row_src, rank_src_a.col_src,
                             rank_src_c.row_src, rank_src_c.col_src
                         );
        std::cout << params << std::endl;
    }
    MPI_Barrier(comm);
#endif

    // get abstract layout descriptions for ScaLAPACK layout
    auto scalapack_layout_a = costa::get_scalapack_layout<T>(
        lld_a,
        {mat_dim_a.rows, mat_dim_a.cols},
        {ia, ja},
        {a_subm, a_subn},
        {b_dim_a.rows, b_dim_a.cols},
        {procrows, proccols},
        ordering_a,
        {rank_src_a.row_src, rank_src_a.col_src},
        a,
        'C',
        rank);

    auto scalapack_layout_c = costa::get_scalapack_layout<T>(
        lld_c,
        {mat_dim_c.rows, mat_dim_c.cols},
        {ic, jc},
        {c_subm, c_subn},
        {b_dim_c.rows, b_dim_c.cols},
        {procrows, proccols},
        ordering_c,
        {rank_src_c.row_src, rank_src_c.col_src},
        c,
        'C',
        rank);

    // transform A to C
    costa::transform<T>(scalapack_layout_a, scalapack_layout_c, comm);

    // print the profiling data
    if (rank == 0) {
        PP();
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
                            const int ctxt);

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
                           const int ctxt);

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
                               const int ctxt);

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
                              const int ctxt);
} // namespace costa
