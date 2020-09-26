#pragma once
#ifdef __cplusplus
extern "C" {
#endif

// scalapack api
void psgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void pdgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void pcgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void pzgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void pigemr2d(int *m, int *n,
              int *a,
              int *ia, int *ja,
              int *desca,
              int *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

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
               int *ictxt);

void pdgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void pcgemr2d_(int *m, int *n,
               float *a,
               int *ia, int *ja,
               int *desca,
               float *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void pzgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void pigemr2d_(int *m, int *n,
               int *a,
               int *ia, int *ja,
               int *desca,
               int *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

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
                int *ictxt);

void pdgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void pcgemr2d__(int *m, int *n,
                float *a,
                int *ia, int *ja,
                int *desca,
                float *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void pzgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void pigemr2d__(int *m, int *n,
                int *a,
                int *ia, int *ja,
                int *desca,
                int *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

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
              int *ictxt);

void PDGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void PCGEMR2D(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void PZGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void PIGEMR2D(int *m, int *n,
              int *a,
              int *ia, int *ja,
              int *desca,
              int *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

#ifdef __cplusplus
}
#endif
