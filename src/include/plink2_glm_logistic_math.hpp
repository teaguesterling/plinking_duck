#pragma once

// Extracted logistic/Firth regression math functions from plink2_glm_logistic.cc.
// These are the core numerical solvers from plink2's --glm implementation.

// cmath must be included before plink2_matrix.h because the header uses fabs()
// in inline functions (HalfInvertedDet) before the .cc file gets to include math.h.
#include <cmath>
#include "plink2_matrix.h"

namespace plink2 {

// Convert 2-bit genotype array to float, removing missing (genotype==3) entries.
// Returns count of non-missing samples written to dst.
uint32_t GenoarrToFloatsRemoveMissing(const uintptr_t *genoarr, const float *__restrict table, uint32_t sample_ct,
                                      float *__restrict dst);

// IRLS logistic regression (float precision).
// Returns 1 on convergence failure, 0 otherwise.
// coef is input/output (starting point, overwritten with betas).
// xx is covariate-major, rows vector-aligned, trailing elements zeroed.
// yy is case/control phenotype, trailing elements zeroed.
BoolErr LogisticRegressionF(const float *yy, const float *xx, const float *sample_offsets, uint32_t sample_ct,
                            uint32_t predictor_ct, float *__restrict coef, uint32_t *is_unfinished_ptr,
                            float *__restrict ll, float *__restrict pp, float *__restrict vv, float *__restrict hh,
                            float *__restrict grad, float *__restrict dcoef);

// Firth penalized logistic regression (float precision).
// Returns 1 on convergence failure, 0 otherwise.
// beta is input/output. hh is output (inverted variance-covariance after completion).
BoolErr FirthRegressionF(const float *yy, const float *xx, const float *sample_offsets, uint32_t sample_ct,
                         uint32_t predictor_ct, float *beta, uint32_t *is_unfinished_ptr, float *hh,
                         double *half_inverted_buf, MatrixInvertBuf1 *inv_1d_buf, double *dbl_2d_buf, float *pp,
                         float *vv, float *ustar, float *delta, float *hdiag, float *ww, float *hh0, float *tmpnxk_buf);

// Hessian computation: hh = mm * diag(vv) * mm^T (float, covariate-major layout).
void ComputeHessianF(const float *mm, const float *vv, uint32_t col_ct, uint32_t row_ct, float *__restrict dest);

// Cholesky decomposition of float matrix (covariate-major, vector-aligned rows).
void CholeskyDecompositionF(const float *aa, uint32_t predictor_ct, float *__restrict ll);

// Solve LL^T x = y for float matrix (covariate-major, vector-aligned rows).
void SolveLinearSystemF(const float *ll, const float *yy, uint32_t predictor_ct, float *__restrict xx);

// Check log-likelihood for convergence issues (returns 1 if any p==0 or p==1).
BoolErr ComputeLoglikCheckedF(const float *yy, const float *pp, uint32_t sample_ct, double *loglik_ptr);

} // namespace plink2
