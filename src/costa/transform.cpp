#include <costa/transform.hpp>
#include <grid2grid/grid_layout.hpp>
#include <grid2grid/transform.hpp>
#include <grid2grid/transformer.hpp>

namespace costa {
template <typename T>
grid2grid::grid_layout<T> get_layout(int n_ranks,
                                     const ::layout_t *layout) {

    // Create the local blocks
    std::vector<grid2grid::block<T>> loc_blks;

    // Create blocks
    for (int i = 0; i < layout->nlocalblocks; ++i) {
        auto &block = layout->localblocks[i];
        auto row = block.row;
        auto col = block.col;
        auto ptr = reinterpret_cast<T *>(block.data);
        auto stride = block.ld;

        grid2grid::block_coordinates coord{row, col};
        grid2grid::interval rows{layout->grid->rowsplit[row],
                                 layout->grid->rowsplit[row + 1]};
        grid2grid::interval cols{layout->grid->colsplit[col],
                                 layout->grid->colsplit[col + 1]};
        loc_blks.emplace_back(rows, cols, coord, ptr, stride);
    }

    // Grid specification
    std::vector<int> rows_split(layout->grid->rowblocks + 1);
    std::copy_n(layout->grid->rowsplit, rows_split.size(), rows_split.begin());

    std::vector<int> cols_split(layout->grid->colblocks + 1);
    std::copy_n(layout->grid->colsplit, cols_split.size(), cols_split.begin());

    std::vector<std::vector<int>> owners_matrix(layout->grid->rowblocks);
    for (int i = 0; i < layout->grid->rowblocks; ++i) {
        owners_matrix[i].resize(layout->grid->colblocks);
        for (int j = 0; j < layout->grid->colblocks; ++j)
            owners_matrix[i][j] = layout->grid->owners[j * layout->grid->rowblocks + i];
    }

    return {{{std::move(rows_split), std::move(cols_split)},
             std::move(owners_matrix),
             n_ranks},
            {std::move(loc_blks)}};
}

template <typename T>
void transform(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const T alpha, const T beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              ) {
    // communicator info
    int P;
    MPI_Comm_size(comm, &P);

    // create grid2grid::grid_layout object from the frontend description
    auto in_layout = get_layout<T>(P, A);
    auto out_layout = get_layout<T>(P, B);

    // transform A to B
    grid2grid::transform<T>(in_layout, out_layout, 
            transpose_or_conjugate, alpha, beta, comm);
}

template <typename T>
void transform_multiple(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const T* alpha, const T* beta,
               // transpose flags
               const char* trans,
               const MPI_Comm comm
              ) {

    // communicator info
    int P;
    MPI_Comm_size(comm, &P);

    // transformer
    grid2grid::transformer<T> transf(comm);

    // schedule all transforms
    for (int i = 0; i < nlayouts; ++i) {
        // create grid2grid::grid_layout object from the frontend description
        auto in_layout = get_layout<T>(P, &A[i]);
        auto out_layout = get_layout<T>(P, &B[i]);

        // schedule the transformation
        transf.schedule(in_layout, out_layout, trans[i], alpha[i], beta[i]);
    }

    // perform the full transformation
    transf.transform();
}

// ***********************************
// template instantiation transform
// ***********************************
template
void transform<int>(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const int alpha, const int beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform<float>(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const float alpha, const float beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform<double>(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const double alpha, const double beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform<std::complex<float>>(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const std::complex<float> alpha, 
               const std::complex<float> beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform<std::complex<double>>(
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const std::complex<double> alpha, 
               const std::complex<double> beta,
               // transpose flags
               const char transpose_or_conjugate,
               const MPI_Comm comm
              );

// ***********************************
// template instantiation transform
// ***********************************
template
void transform_multiple<int>(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const int* alpha, const int* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform_multiple<float>(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const float* alpha, const float* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform_multiple<double>(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const double* alpha, const double* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform_multiple<std::complex<float>>(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const std::complex<float>* alpha, 
               const std::complex<float>* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const MPI_Comm comm
              );

template
void transform_multiple<std::complex<double>>(
               const int nlayouts,
               const layout_t* A,
               const layout_t* B,
               // scaling parameters
               const std::complex<double>* alpha, 
               const std::complex<double>* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const MPI_Comm comm
              );
}
