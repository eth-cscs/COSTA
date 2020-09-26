#include <costa/costa_pxgemr2d.hpp>

extern "C" {
#include <costa/pxgemr2d.h>

void pigemr2d(int *m, int *n,
              int *a,
              int *ia, int *ja,
              int *desca,
              int *c,
              int *ic, int *jc,
              int *descc,
              int *ictxt) {
    costa::pxgemr2d<int>(
                  *m,
                  *n,
                  a,
                  *ia,
                  *ja,
                  desca,
                  c,
                  *ic,
                  *jc,
                  descc,
                  ictxt);
}

void psgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *c,
              int *ic, int *jc,
              int *descc,
              int *ictxt) {
    costa::pxgemr2d<float>(
                  *m,
                  *n,
                  a,
                  *ia,
                  *ja,
                  desca,
                  c,
                  *ic,
                  *jc,
                  descc,
                  ictxt);
}

void pdgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *c,
              int *ic, int *jc,
              int *descc,
              int *ictxt) {
    costa::pxgemr2d<double>(
                  *m,
                  *n,
                  a,
                  *ia,
                  *ja,
                  desca,
                  c,
                  *ic,
                  *jc,
                  descc,
                  ictxt);
}

void pcgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *c,
              int *ic, int *jc,
              int *descc,
              int *ictxt) {
    costa::pxgemr2d<std::complex<float>>(
                  *m,
                  *n,
                  reinterpret_cast<const std::complex<float>*>(a),
                  *ia,
                  *ja,
                  desca,
                  reinterpret_cast<std::complex<float>*>(c),
                  *ic,
                  *jc,
                  descc,
                  ictxt);
}

void pzgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *c,
              int *ic, int *jc,
              int *descc,
              int *ictxt) {
    costa::pxgemr2d<std::complex<double>>(
                  *m,
                  *n,
                  reinterpret_cast<const std::complex<double>*>(a),
                  *ia,
                  *ja,
                  desca,
                  reinterpret_cast<std::complex<double>*>(c),
                  *ic,
                  *jc,
                  descc,
                  ictxt);
}

// *********************************************************************************
// Same as previously, but with added underscore at the end.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void psgemr2d_(int *m, int *n,
               float *a,
               int *ia, int *ja,
               int *desca,
               float *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt) {
    psgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pdgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt) {
    pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pcgemr2d_(int *m, int *n,
               float *a,
               int *ia, int *ja,
               int *desca,
               float *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt) {
    pcgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pzgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt) {
    pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pigemr2d_(int *m, int *n,
               int *a,
               int *ia, int *ja,
               int *desca,
               int *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt) {
    pigemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

// *********************************************************************************
// Same as previously, but with added double underscores at the end.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void psgemr2d__(int *m, int *n,
                float *a,
                int *ia, int *ja,
                int *desca,
                float *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt) {
    psgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pdgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt) {
    pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pcgemr2d__(int *m, int *n,
                float *a,
                int *ia, int *ja,
                int *desca,
                float *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt) {
    pcgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pzgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt) {
    pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void pigemr2d__(int *m, int *n,
                int *a,
                int *ia, int *ja,
                int *desca,
                int *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt) {
    pigemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

// *********************************************************************************
// Same as previously, but CAPITALIZED.
// This is used for fortran interfaces, in case fortran expects these symbols
// *********************************************************************************
void PSGEMR2D(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt) {
    psgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void PDGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt) {
    pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void PCGEMR2D(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt) {
    pcgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void PZGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt) {
    pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}

void PIGEMR2D(int *m, int *n,
              int *a,
              int *ia, int *ja,
              int *desca,
              int *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt) {
    pigemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
}
}
