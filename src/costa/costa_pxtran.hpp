#pragma once
#include <complex>
#include <costa/scalapack.hpp>
/*
 * Transposing the matrices
 */
namespace costa {

using zdouble_t = std::complex<double>;
using zfloat_t = std::complex<float>;

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
           T *c,
           const int ic,
           const int jc,
           const int *descc);

// scales the submatrix of C by beta
// The submatrix is defined by (ic-1, jc-1) and (ic-1+m, jc-1+n)
template <typename T>
void scale_matrix(const int* descc, T* c,
                  const int ic, const int jc,
                  const int m, const int n,
                  const T beta);
} // namespace costa

