#pragma once
#ifdef __cplusplus
extern "C" {
#endif

// scalapack api
void costa_psgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void costa_pdgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void costa_pcgemr2d(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void costa_pzgemr2d(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void costa_pigemr2d(int *m, int *n,
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
void costa_psgemr2d_(int *m, int *n,
               float *a,
               int *ia, int *ja,
               int *desca,
               float *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void costa_pdgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void costa_pcgemr2d_(int *m, int *n,
               float *a,
               int *ia, int *ja,
               int *desca,
               float *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void costa_pzgemr2d_(int *m, int *n,
               double *a,
               int *ia, int *ja,
               int *desca,
               double *b,
               int *ib, int *jb,
               int *descb,
               int *ictxt);

void costa_pigemr2d_(int *m, int *n,
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
void costa_psgemr2d__(int *m, int *n,
                float *a,
                int *ia, int *ja,
                int *desca,
                float *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void costa_pdgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void costa_pcgemr2d__(int *m, int *n,
                float *a,
                int *ia, int *ja,
                int *desca,
                float *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void costa_pzgemr2d__(int *m, int *n,
                double *a,
                int *ia, int *ja,
                int *desca,
                double *b,
                int *ib, int *jb,
                int *descb,
                int *ictxt);

void costa_pigemr2d__(int *m, int *n,
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
void COSTA_PSGEMR2D(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void COSTA_PDGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void COSTA_PCGEMR2D(int *m, int *n,
              float *a,
              int *ia, int *ja,
              int *desca,
              float *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void COSTA_PZGEMR2D(int *m, int *n,
              double *a,
              int *ia, int *ja,
              int *desca,
              double *b,
              int *ib, int *jb,
              int *descb,
              int *ictxt);

void COSTA_PIGEMR2D(int *m, int *n,
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
