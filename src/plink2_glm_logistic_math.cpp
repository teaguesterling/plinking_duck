// Extracted logistic/Firth regression math functions from plink2_glm_logistic.cc.
//
// Original source: third_party/plink-ng/2.0/plink2_glm_logistic.cc lines 147-1012
// Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang.
// Licensed under GPL v3+.
//
// Contains the core numerical solvers used by plink2 --glm:
//   - LogisticRegressionF (IRLS logistic regression, float precision)
//   - FirthRegressionF (penalized IRLS, logistf port)
//   - Supporting math functions (Hessian, Cholesky, etc.)

#include "plink2_glm_logistic_math.hpp"

#include <pgenlib_misc.h>  // PAIR_TABLE16
#include <plink2_bits.h>
#include <plink2_fmath.h>

#include <assert.h>
#include <math.h>
#include <string.h>

namespace plink2 {

// ---------------------------------------------------------------------------
// Lookup tables for genotype-to-float conversion
// ---------------------------------------------------------------------------

static const float kSmallFloatPairs[32] = PAIR_TABLE16(0.0, 1.0, 2.0, 3.0);

static const float kSmallInvFloatPairs[32] = PAIR_TABLE16(2.0, 1.0, 0.0, 3.0);

static const float kSmallInvFloats[4] = {2.0, 1.0, 0.0, 3.0};

uint32_t GenoarrToFloatsRemoveMissing(const uintptr_t* genoarr, const float* __restrict table, uint32_t sample_ct, float* __restrict dst) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  float* dst_iter = dst;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return S_CAST(uint32_t, dst_iter - dst);
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
        *dst_iter++ = table[cur_geno];
      }
      geno_word >>= 2;
    }
  }
}

// ---------------------------------------------------------------------------
// Architecture-specific SIMD helpers
// ---------------------------------------------------------------------------

#ifdef __LP64__
// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.

#  ifdef FVEC_32
static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  uint32_t row_idx = 0;
  __m256 s1;
  __m256 s2;
  __m256 s3;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (; row_idx < row_ctm3; row_idx += 4) {
      s1 = _mm256_setzero_ps();
      s2 = _mm256_setzero_ps();
      s3 = _mm256_setzero_ps();
      __m256 s4 = _mm256_setzero_ps();
      for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
        const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
        const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
        __m256 a1 = _mm256_load_ps(mm_ptr);
        __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
        __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
        __m256 a4 = _mm256_load_ps(&(mm_ptr[3 * col_ctav]));
        s1 = _mm256_fmadd_ps(a1, vv, s1);
        s2 = _mm256_fmadd_ps(a2, vv, s2);
        s3 = _mm256_fmadd_ps(a3, vv, s3);
        s4 = _mm256_fmadd_ps(a4, vv, s4);
      }
      *dest++ = VecFHsum(R_CAST(VecF, s1));
      *dest++ = VecFHsum(R_CAST(VecF, s2));
      *dest++ = VecFHsum(R_CAST(VecF, s3));
      *dest++ = VecFHsum(R_CAST(VecF, s4));
    }
  }
  s1 = _mm256_setzero_ps();
  s2 = _mm256_setzero_ps();
  s3 = _mm256_setzero_ps();
  switch (row_ct % 4) {
  case 3:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
      s3 = _mm256_fmadd_ps(a3, vv, s3);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest++ = VecFHsum(R_CAST(VecF, s2));
    *dest = VecFHsum(R_CAST(VecF, s3));
    break;
  case 2:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest = VecFHsum(R_CAST(VecF, s2));
    break;
  case 1:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(&(mm[row_idx * col_ctav + col_idx]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
    }
    *dest = VecFHsum(R_CAST(VecF, s1));
    break;
  }
}

#  else  // !FVEC_32
static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

#  endif  // !FVEC_32

static inline void LogisticSseF(uint32_t nn, float* vect) {
  const VecF zero = vecf_setzero();
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(VecF*, &(vect[uii]));
    aa = zero - aa;
    aa = fmath_exp_ps(aa);
    aa = aa + one;
    aa = one / aa;
    *R_CAST(VecF*, &(vect[uii])) = aa;
  }
}

static inline void ComputeVAndPMinusYF(const float* yy, uint32_t nn, float* __restrict pp, float* __restrict vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
    VecF ytmp = *R_CAST(const VecF*, &(yy[uii]));
    *R_CAST(VecF*, &(pp[uii])) = ptmp - ytmp;
  }
}

static inline void ComputeVF(const float* pp, uint32_t nn, float* __restrict vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(const VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
  }
}

static inline float TripleProductF(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  VecF sum = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(const VecF*, &(v1[uii]));
    VecF bb = *R_CAST(const VecF*, &(v2[uii]));
    VecF cc = *R_CAST(const VecF*, &(v3[uii]));
    sum = sum + aa * bb * cc;
  }
  return VecFHsum(sum);
}

static inline void ComputeTwoDiagTripleProductF(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
  VecF saa = vecf_setzero();
  VecF sab = vecf_setzero();
  VecF sbb = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF atmp = *R_CAST(const VecF*, &(aa[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF av = atmp * vtmp;
    const VecF bv = btmp * vtmp;
    saa = saa + atmp * av;
    sab = sab + atmp * bv;
    sbb = sbb + btmp * bv;
  }
  *raa_ptr = VecFHsum(saa);
  *rab_ptr = VecFHsum(sab);
  *rbb_ptr = VecFHsum(sbb);
}

static inline void ComputeThreeTripleProductF(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF a3tmp = *R_CAST(const VecF*, &(a3[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    btmp = btmp * vtmp;
    s1 = s1 + a1tmp * btmp;
    s2 = s2 + a2tmp * btmp;
    s3 = s3 + a3tmp * btmp;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}

static inline void ComputeTwoPlusOneTripleProductF(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF bv = btmp * vtmp;
    s1 = s1 + btmp * bv;
    s2 = s2 + a1tmp * bv;
    s3 = s3 + a2tmp * bv;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}

static void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* __restrict dst) {
  const uintptr_t fullvec_ct = ct / kFloatPerFVec;
  const VecF* src_alias = R_CAST(const VecF*, src);
  VecF vsum = vecf_setzero();
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    vsum += src_alias[vidx];
  }
  float sum = VecFHsum(vsum);
  const uintptr_t trailing_start_idx = fullvec_ct * kFloatPerFVec;
  for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
    sum += src[ulii];
  }

  const float neg_mean = -sum / S_CAST(float, ct);
  const VecF neg_vmean = VCONST_F(neg_mean);
  VecF* dst_alias = R_CAST(VecF*, dst);
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    dst_alias[vidx] = src_alias[vidx] + neg_vmean;
  }
  if (trailing_start_idx != ct) {
    for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
      dst[ulii] = src[ulii] + neg_mean;
    }
    const uintptr_t trailing_stop_idx = trailing_start_idx + kFloatPerFVec;
    for (uintptr_t ulii = ct; ulii != trailing_stop_idx; ++ulii) {
      dst[ulii] = S_CAST(float, 0.0);
    }
  }
}

#else  // no __LP64__
static inline void LogisticSseF(uint32_t nn, float* vect) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vect[uii] = S_CAST(float, 1.0) / (1 + expf(-vect[uii]));
  }
}

static inline void ComputeVAndPMinusYF(const float* yy, uint32_t nn, float* __restrict pp, float* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void ComputeVF(const float* pp, uint32_t nn, float* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
  }
}

static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

static inline float TripleProductF(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    fxx += v1[uii] * v2[uii] * v3[uii];
  }
  return fxx;
}

static inline void ComputeTwoDiagTripleProductF(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = aa[uii];
    const float fyy = bb[uii];
    float fzz = vv[uii];
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void ComputeThreeTripleProductF(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = bb[uii] * vv[uii];
    r1 += a1[uii] * fxx;
    r2 += a2[uii] * fxx;
    r3 += a3[uii] * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void ComputeTwoPlusOneTripleProductF(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = bb[uii];
    const float fyy = fxx * vv[uii];
    r1 += fxx * fyy;
    r2 += a1[uii] * fyy;
    r3 += a2[uii] * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* __restrict dst) {
  float sum = 0.0;
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    sum += src[ulii];
  }
  const float mean = sum / u31tof(ct);
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = src[ulii] - mean;
  }
}
#endif

// ---------------------------------------------------------------------------
// Core math functions (architecture-independent)
// ---------------------------------------------------------------------------

BoolErr ComputeLoglikCheckedF(const float* yy, const float* pp, uint32_t sample_ct, double* loglik_ptr) {
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = S_CAST(double, pp[sample_idx]);
    if ((new_pi == 0.0) || (new_pi == 1.0)) {
      return 1;
    }
    loglik += (yy[sample_idx] != S_CAST(float, 0.0))? log(new_pi) : log1p(-new_pi);
  }
  *loglik_ptr = loglik;
  return 0;
}

void ComputeHessianF(const float* mm, const float* vv, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  const uintptr_t row_ctav = RoundUpPow2(row_ct, kFloatPerFVec);
  const uintptr_t row_ctavp1 = row_ctav + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const float* mm_cur = &(mm[row_idx * col_ctav]);
      ComputeTwoDiagTripleProductF(mm_cur, &(mm_cur[col_ctav]), vv, col_ct, &(dest[row_idx * row_ctavp1]), &(dest[(row_idx + 1) * row_ctavp1 - 1]), &(dest[(row_idx + 1) * row_ctavp1]));
      ComputeTwoPlusOneTripleProductF(&(mm_cur[2 * col_ctav]), &(mm_cur[col_ctav]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_ctavp1]), &(dest[(row_idx + 2) * row_ctavp1 - 1]), &(dest[(row_idx + 2) * row_ctavp1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 != row_ct; ++row_idx2) {
        ComputeThreeTripleProductF(&(mm[row_idx2 * col_ctav]), mm_cur, &(mm_cur[col_ctav]), &(mm_cur[2 * col_ctav]), vv, col_ct, &(dest[row_idx2 * row_ctav + row_idx]), &(dest[row_idx2 * row_ctav + row_idx + 1]), &(dest[row_idx2 * row_ctav + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    ComputeTwoPlusOneTripleProductF(&(mm[(row_ct - 3) * col_ctav]), &(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 3) * row_ctavp1]), &(dest[(row_ct - 2) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1 - 2]));
    // fall through
  case 2:
    ComputeTwoDiagTripleProductF(&(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 2) * row_ctavp1]), &(dest[(row_ct - 1) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_ctavp1] = TripleProductF(&(mm[(row_ct - 1) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct);
  }
}

void CholeskyDecompositionF(const float* aa, uint32_t predictor_ct, float* __restrict ll) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t predictor_ctavp1 = predictor_ctav + 1;
  const float* aa_diag_elem_ptr = aa;
  float* cur_ll_row = ll;
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = *aa_diag_elem_ptr;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const float fyy = cur_ll_row[col_idx];
      fxx -= fyy * fyy;
    }
    float fyy;
    if (fxx >= S_CAST(float, 0.0)) {
      fyy = sqrtf(fxx);
    } else {
      fyy = S_CAST(float, 1e-6);
    }
    cur_ll_row[row_idx] = fyy;
    fyy = S_CAST(float, 1.0) / fyy;
    const float* aa_col_iter = aa_diag_elem_ptr;
    float* cur_ll_row2 = cur_ll_row;
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != predictor_ct; ++row_idx2) {
      aa_col_iter = &(aa_col_iter[predictor_ctav]);
      float fxx2 = *aa_col_iter;
      cur_ll_row2 = &(cur_ll_row2[predictor_ctav]);
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        fxx2 -= cur_ll_row[col_idx] * cur_ll_row2[col_idx];
      }
      cur_ll_row2[row_idx] = fxx2 * fyy;
    }
    aa_diag_elem_ptr = &(aa_diag_elem_ptr[predictor_ctavp1]);
    cur_ll_row = &(cur_ll_row[predictor_ctav]);
  }
}

void SolveLinearSystemF(const float* ll, const float* yy, uint32_t predictor_ct, float* __restrict xx) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  {
    const float* cur_ll_row = ll;
    for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
      float fxx = yy[row_idx];
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        fxx -= cur_ll_row[col_idx] * xx[col_idx];
      }
      xx[row_idx] = fxx / cur_ll_row[row_idx];
      cur_ll_row = &(cur_ll_row[predictor_ctav]);
    }
  }
  for (uint32_t col_idx = predictor_ct; col_idx; ) {
    float* xx_stop = &(xx[--col_idx]);
    float fxx = *xx_stop;
    const float* ll_col_iter = &(ll[(predictor_ct - 1) * predictor_ctav + col_idx]);
    for (float* xx_iter = &(xx[predictor_ct - 1]); xx_iter != xx_stop; --xx_iter) {
      fxx -= (*ll_col_iter) * (*xx_iter);
      ll_col_iter -= predictor_ctav;
    }
    *xx_stop = fxx / (*ll_col_iter);
  }
}

// ---------------------------------------------------------------------------
// LogisticRegressionF — IRLS logistic regression
// ---------------------------------------------------------------------------

BoolErr LogisticRegressionF(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* __restrict coef, uint32_t* is_unfinished_ptr, float* __restrict ll, float* __restrict pp, float* __restrict vv, float* __restrict hh, float* __restrict grad, float* __restrict dcoef) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  float min_delta_coef = 1e9;

  ZeroFArr(predictor_ct * predictor_ctav, ll);
  ZeroFArr(sample_ctav - sample_ct, &(pp[sample_ct]));
  ZeroFArr(sample_ctav - sample_ct, &(vv[sample_ct]));
  for (uint32_t iteration = 0; ; ++iteration) {
    ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddFVec(sample_offsets, sample_ctav, pp);
    }

    LogisticSseF(sample_ct, pp);

    ComputeVAndPMinusYF(yy, sample_ct, pp, vv);

    ComputeHessianF(xx, vv, sample_ct, predictor_ct, hh);

    MultMatrixDxnVectNF(xx, pp, sample_ct, predictor_ct, grad);

    CholeskyDecompositionF(hh, predictor_ct, ll);

    SolveLinearSystemF(ll, grad, predictor_ct, dcoef);

    float delta_coef = 0.0;
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      const float cur_dcoef = dcoef[pred_idx];
      delta_coef += fabsf(cur_dcoef);
      coef[pred_idx] -= cur_dcoef;
    }
    if (delta_coef < min_delta_coef) {
      min_delta_coef = delta_coef;
    }
    if (delta_coef != delta_coef) {
      return 1;
    }
    if (iteration > 3) {
      if (((delta_coef > S_CAST(float, 20.0)) && (delta_coef > 2 * min_delta_coef)) || ((iteration > 6) && (fabsf(S_CAST(float, 1.0) - delta_coef) < S_CAST(float, 1e-3)))) {
        return 1;
      }
      if (iteration > 13) {
        for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
          if (fabsf(coef[pred_idx]) > S_CAST(float, 8e3)) {
            return 1;
          }
        }
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    if (delta_coef < S_CAST(float, 1e-4)) {
      for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
        if (fabsf(coef[pred_idx]) > S_CAST(float, 6e4)) {
          return 1;
        }
      }
      return 0;
    }
  }
}

// ---------------------------------------------------------------------------
// Firth helpers
// ---------------------------------------------------------------------------

static void FirthComputeHdiagWeightsF(const float* yy, const float* xx, const float* pp, const float* hh, const float* vv, uint32_t predictor_ct, uint32_t predictor_ctav, uint32_t sample_ct, uint32_t sample_ctav, float* hdiag, float* ww, float* tmpnxk_buf) {
  ColMajorFmatrixMultiplyStrided(xx, hh, sample_ct, sample_ctav, predictor_ct, predictor_ctav, predictor_ct, sample_ctav, tmpnxk_buf);
#ifdef __LP64__
  const VecF half = VCONST_F(0.5);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    VecF dotprods = vecf_setzero();
    const float* xx_row = &(xx[sample_offset]);
    const float* tmpnxk_row = &(tmpnxk_buf[sample_offset]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const VecF cur_xx = *R_CAST(const VecF*, &(xx_row[pred_uidx * sample_ctav]));
      const VecF cur_tmpnxk = *R_CAST(const VecF*, &(tmpnxk_row[pred_uidx * sample_ctav]));
      dotprods = dotprods + cur_xx * cur_tmpnxk;
    }
    const VecF cur_vv = *R_CAST(const VecF*, &(vv[sample_offset]));
    const VecF cur_pi = *R_CAST(const VecF*, &(pp[sample_offset]));
    const VecF cur_yy = *R_CAST(const VecF*, &(yy[sample_offset]));
    const VecF cur_hdiag = cur_vv * dotprods;
    *R_CAST(VecF*, &(hdiag[sample_offset])) = cur_hdiag;
    const VecF half_minus_cur_pis = half - cur_pi;
    const VecF yy_minus_cur_pis = cur_yy - cur_pi;
    *R_CAST(VecF*, &(ww[sample_offset])) = yy_minus_cur_pis + cur_hdiag * half_minus_cur_pis;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    float dotprod = 0.0;
    const float* xx_row = &(xx[sample_idx]);
    const float* tmpnxk_row = &(tmpnxk_buf[sample_idx]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      dotprod += xx_row[pred_uidx * sample_ctav] * tmpnxk_row[pred_uidx * sample_ctav];
    }
    const float cur_hdiag = vv[sample_idx] * dotprod;
    hdiag[sample_idx] = cur_hdiag;
    const float cur_pi = pp[sample_idx];
    ww[sample_idx] = (yy[sample_idx] - cur_pi) + cur_hdiag * (S_CAST(float, 0.5) - cur_pi);
  }
#endif
}

static void FirthComputeSecondWeightsF(const float* hdiag, const float* vv, __maybe_unused uint32_t sample_ct, __maybe_unused uint32_t sample_ctav, float* ww) {
#ifdef __LP64__
  const VecF one = VCONST_F(1.0);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    const VecF cur_hdiag = *R_CAST(const VecF*, &(hdiag[sample_offset]));
    const VecF cur_vv = *R_CAST(const VecF*, &(vv[sample_offset]));
    *R_CAST(VecF*, &(ww[sample_offset])) = (one + cur_hdiag) * cur_vv;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    ww[sample_idx] = (S_CAST(float, 1.0) + hdiag[sample_idx]) * vv[sample_idx];
  }
#endif
}

// ---------------------------------------------------------------------------
// FirthRegressionF — penalized IRLS (logistf port)
// ---------------------------------------------------------------------------

BoolErr FirthRegressionF(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* beta, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* ustar, float* delta, float* hdiag, float* ww, float* hh0, float* tmpnxk_buf) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);

  ZeroFArr(predictor_ctav - predictor_ct, &(ustar[predictor_ct]));

  const uint32_t trailing_sample_ct = sample_ctav - sample_ct;
  if (trailing_sample_ct) {
    ZeroFArr(trailing_sample_ct, &(pp[sample_ct]));
    ZeroFArr(trailing_sample_ct, &(vv[sample_ct]));
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      ZeroFArr(trailing_sample_ct, &(tmpnxk_buf[sample_ct + pred_idx * sample_ctav]));
    }
  }

  const uint32_t max_iter = 25;
  const float gconv = S_CAST(float, 0.0001);
  const float xconv = S_CAST(float, 0.0001);
  const double lconv = 0.0001;
  float delta_max = 0.0;
  double loglik_old = 0.0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    ColMajorFmatrixVectorMultiplyStrided(xx, beta, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddFVec(sample_offsets, sample_ctav, pp);
    }
    LogisticSseF(sample_ct, pp);
    double loglik;
    if (ComputeLoglikCheckedF(yy, pp, sample_ct, &loglik)) {
      return 1;
    }
    ComputeVF(pp, sample_ct, vv);

    ComputeHessianF(xx, vv, sample_ct, predictor_ct, hh0);

    if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh0, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    double dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct, predictor_ct);
    loglik += 0.5 * log(dethh);

    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh0, inv_1d_buf, dbl_2d_buf);

    ReflectFmatrix0(predictor_ct, predictor_ctav, hh0);

    FirthComputeHdiagWeightsF(yy, xx, pp, hh0, vv, predictor_ct, predictor_ctav, sample_ct, sample_ctav, hdiag, ww, tmpnxk_buf);

    MultMatrixDxnVectNF(xx, ww, sample_ct, predictor_ct, ustar);
    if (iter_idx) {
      float ustar_max = 0.0;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        const float abs_ustar_cur = fabsf(ustar[pred_uidx]);
        if (abs_ustar_cur > ustar_max) {
          ustar_max = abs_ustar_cur;
        }
      }
      const double loglik_change = loglik - loglik_old;
      if ((delta_max <= xconv) && (ustar_max < gconv) && (loglik_change < lconv)) {
        return 0;
      }
      if (iter_idx > max_iter) {
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    loglik_old = loglik;

    FirthComputeSecondWeightsF(hdiag, vv, sample_ct, sample_ctav, ww);
    ComputeHessianF(xx, ww, sample_ct, predictor_ct, hh);
    if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh, inv_1d_buf, dbl_2d_buf);
    ReflectFmatrix0(predictor_ct, predictor_ctav, hh);

    MultMatrixDxnVectNF(hh, ustar, predictor_ct, predictor_ct, delta);

    delta_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_delta_cur = fabsf(delta[pred_uidx]);
      if (abs_delta_cur > delta_max) {
        delta_max = abs_delta_cur;
      }
    }
    const float maxstep = 5.0;
    if (delta_max > maxstep) {
      const float scaling_factor = maxstep / delta_max;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        delta[pred_uidx] *= scaling_factor;
      }
      delta_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      beta[pred_uidx] += delta[pred_uidx];
    }
  }
}

} // namespace plink2
