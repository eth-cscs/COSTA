#include <costa/pxtran_op/costa_pxtran_op.hpp>

extern "C" {
#include <costa/pxtran/pxtran.h>

// Reimplement ScaLAPACK signatures functions
void pdtran(const int *m , const int *n , 
            double *alpha , const double *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const double *beta , double *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    costa::pxtran_op<double>(
                  *m,
                  *n,
                  *alpha,
                  a,
                  *ia,
                  *ja,
                  desca,
                  *beta,
                  c,
                  *ic,
                  *jc,
                  descc,
                  'T');
}

void pstran(const int *m , const int *n , 
            float *alpha , const float *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const float *beta , float *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    costa::pxtran_op<float>(
                  *m,
                  *n,
                  *alpha,
                  a,
                  *ia,
                  *ja,
                  desca,
                  *beta,
                  c,
                  *ic,
                  *jc,
                  descc,
                  'T');
}

// *********************************************************************************
// Same as previously, but with added underscore at the end.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void pstran_(const int *m , const int *n , 
            float *alpha , const float *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const float *beta , float *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    pstran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

void pdtran_(const int *m , const int *n , 
            double *alpha , const double *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const double *beta , double *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

// *********************************************************************************
// Same as previously, but with added double underscores at the end.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void pstran__(const int *m , const int *n , 
              float *alpha , const float *a , 
              const int *ia , const int *ja , 
              const int *desca , 
              const float *beta , float *c , 
              const int *ic , const int *jc ,
              const int *descc ) {
    pstran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

void pdtran__(const int *m , const int *n , 
              double *alpha , const double *a , 
              const int *ia , const int *ja , 
              const int *desca , 
              const double *beta , double *c , 
              const int *ic , const int *jc ,
              const int *descc ) {
    pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

// *********************************************************************************
// Same as previously, but CAPITALIZED.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void PSTRAN(const int *m , const int *n , 
            float *alpha , const float *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const float *beta , float *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    pstran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

void PDTRAN(const int *m , const int *n , 
            double *alpha , const double *a , 
            const int *ia , const int *ja , 
            const int *desca , 
            const double *beta , double *c , 
            const int *ic , const int *jc ,
            const int *descc ) {
    pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

}
