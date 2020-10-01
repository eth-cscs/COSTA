#pragma once
#include <complex>
#include <costa/layout.hpp>

/**
 * The general transformation looks as follows (A = input, B = output):
 * B = alpha * (A^T) + beta * B
 * The parameters are as follows:
 * alpha: a scalar to be multiplied by the input matrix
 * beta: a scalar to be multiplied by the output matrix
 * transpose_or_conjugate ('N', 'T' or 'C'): describes whether the input
 * matrix should be left unchanged ('N'), transposed ('T') or conjugated ('C')
 */
namespace costa {

using zdouble_t = std::complex<double>;
using zfloat_t = std::complex<float>;

// applies the transformation: B = alpha * (A^T) + beta * B
// if relabel = true, the ranks are relabeled to minimize the communication
// and the new communicator is returned
// if relabel = false, the ranks are not relabelled and 
// the input communicator is returned without relabelling
template <typename T>
MPI_Comm transform(
               const layout* A,
               const layout* B,
               // scaling parameters
               const T alpha, const T beta,
               // transpose flags
               const char transpose_or_conjugate,
               // if true, ranks will be relabelled to minimize the communication
               // and the new communicator will be returned
               const bool relabel,
               const MPI_Comm comm
              );

// transforming multiple layouts at once (in a single communication round), 
// minimizing the latency. In addition, the relabelling will take into account
// the transformations of all layouts at once.
// Due to the minimized communication latency and also 
// due to the optimal relabelling accross all transformations, 
// this is potentially more efficient than invoking
// single transform for each layout separately.
template <typename T>
MPI_Comm transform_multiple(
               const layout* A,
               const layout* B,
               // scaling parameters
               const T* alpha, const T* beta,
               // transpose flags
               const char* transpose_or_conjugate,
               const bool relabel,
               const MPI_Comm comm
              );
} // namespace costa
