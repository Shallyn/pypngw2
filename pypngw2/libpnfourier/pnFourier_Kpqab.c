/**
* Writer: Xiaolin.liu
* shallyn.liu@foxmail.com
**/
#include "pnFourier_Kpqab.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include "core_Kpqab_Approx.h"


#if defined(__STDC_NO_COMPLEX__) || defined(FFTW_NO_COMPLEX)
  #define COMPLEX(r, i)   { (r), (i) }
  #define GETREAL(z)       ((z)[0])
  #define GETIMAG(z)       ((z)[1])
  #define SETCOMP(z,vr,vi)     ((z)[0]=vr,(z)[1]=vi)
#else
  #define COMPLEX(r, i)   ((REAL8 _Complex)((r) + _Complex_I * (i)))
  #define GETREAL(z)       (creal(z))
  #define GETIMAG(z)       (cimag(z))
  #define SETCOMP(z,vr,vi) (z=(vr)+I*(vi))
#endif

/* ------------- Utils ------------- */
static inline REAL8 safe_pow(INT m, REAL8 x)
{
    REAL8 delta = 1.0 - x;                     
    if (delta > 1e-8)       
        return pow(x, m);
    return exp(m * log1p(-delta)); 
}

static INT choose_N_a(INT p_max_abs, INT q, REAL8 e, REAL8 beta,
                    REAL8 tol, INT a)
{
    /* crude but safe – geometric tail + 12 a guard */
    REAL8 z  = fabs((REAL8)q * e);
    INT   k  = (INT)ceil(log(tol) / log(beta));
    REAL8 term = 1.0;
    while (term > tol) { ++k; term *= z / (2.0 * k); }
    INT N = k + 12 * a;
    if (!(N & 1)) ++N;
    if (N < p_max_abs + 12 * a) N = (p_max_abs + 12 * a) | 1;
    return N;
}

static void bessel_fill(INT kmin, INT kmax, REAL8 z, REAL8 *J)
{
    /* Direct libm evaluation – O(N) but <300 calls in practice.  */
    for (INT k = kmin; k <= kmax; ++k) {
        REAL8 v;
        if (k >= 0) {
            v = jn(k, z);
        } else {
            INT kk = -k;
            v = jn(kk, z);
            if (kk & 1) v = -v;      /* J_{-n}(z)=(-1)^n J_n(z) */
        }
        J[k - kmin] = v;
    }
}

/* ------------ Series Method ------------*/
// REAL8 K_pqa0_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     // if (a==0) return jn(p,q*e);  
//     // if(a==0) return jn(p,q*e);
//     REAL8 beta = eval_beta(e);
//     if (q==0 && a==0) return p==0 ? log((1.+sqrt(1.-e*e))/2.) : -pow(beta, abs(p)) / abs(p);
//     REAL8 beta2 = beta * beta;
//     REAL8 z    = q * e;

//     /* choose truncation N so that tail < tol                          */
//     INT p_abs = abs(p);
//     INT N = choose_N_a(p_abs, q, e, beta, atol, a);

//     INT kmin = -N, kmax = N;
//     INT M = kmax - kmin + 1;

//     /* Bessel array J_k(z) (both signs)                                */
//     REAL8 *Jk = malloc(M * sizeof(REAL8));
//     bessel_fill(kmin, kmax, z, Jk);

//     size_t nMax = abs(kmax) + abs(p); // just guess
//     LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     DLaplaceCache *dlc = CreateDLaplaceCache(nMax, a, beta);
//     REAL8 logprefLc = log(1 + beta2);

//     /* series summation                                                */
//     REAL8 sum = 0.0;
//     for (INT k = kmin; k <= kmax; ++k) {
//         INT idx = k - kmin;
//         REAL8 bess = Jk[idx];
//         if (k & 1) bess = -bess;          /* (-1)^k factor */
//         INT n = abs(k + p);             /* |k+p| */
//         REAL8 Lnp = get_Laplace_from_LaplaceCache(n, lc);
//         REAL8 DLnp = get_DLaplace_from_DLaplaceCache(n, dlc);

//         REAL8 term = bess * (DLnp + logprefLc*Lnp);
//         // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
//         sum += term;
//         if (fabs(term) < rtol * fabs(sum)) {   /* early exit */
//             if (k > p_abs) break;            /* tail symmetric */
//         }
//     }
//     free(Jk);
//     DestroyDLaplaceCache(dlc);
//     REAL8 factor = pow(1.0 + beta2, a);
//     return -factor * sum;   /* already includes 1/(2π) via Laplace */
// }

REAL8 K_pqa0_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }

    // if (a==0) return jn(p,q*e);  
    // if(a==0) return jn(p,q*e);
    REAL8 beta = eval_beta(e);
    if (q==0 && a==0) return p==0 ? log((1.+sqrt(1.-e*e))/2.) : -pow(beta, abs(p)) / abs(p);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT p_abs = abs(p);
    INT N = choose_N_a(p_abs, q, e, beta, atol, a);

    INT kmin = -N, kmax = N;
    INT M = kmax - kmin + 1;

    /* Bessel array J_k(z) (both signs)                                */
    // REAL8 *Jk = malloc(M * sizeof(REAL8));
    // bessel_fill(kmin, kmax, z, Jk);
    BesselJCache *jc = CreateBesselJCache(N, z);
    size_t nMax = abs(kmax) + abs(p); // just guess
    LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    DLaplaceCache *dlc = CreateDLaplaceCache(nMax, a, beta);
    REAL8 logprefLc = log(1 + beta2);

    /* series summation         */
    REAL8 sum_p = 0.0;
    for (INT k = 0; ; ++k) {
        // INT idx = k - kmin;
        // REAL8 bess = Jk[idx];
        REAL8 bess = get_BesselJ_from_BesselJCache(k, jc);
        // if (k & 1) bess = -bess;          /* (-1)^k factor */
        INT n = abs(k + p);             /* |k+p| */
        REAL8 Lnp = get_Laplace_from_LaplaceCache(n, lc);
        REAL8 DLnp = get_DLaplace_from_DLaplaceCache(n, dlc);

        REAL8 term = (k & 1 ? -1.0 : 1.0) * bess * (DLnp + logprefLc*Lnp);
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sum_p += term;
        if (fabs(term) <= rtol * fabs(sum_p)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    REAL8 sum_m = 0.0;
    // negative summation
    for (INT k = 1; ; ++k) {
        // INT idx = k - kmin;
        // REAL8 bess = Jk[idx];
        REAL8 bess = get_BesselJ_from_BesselJCache(-k, jc);
        // if (k & 1) bess = -bess;          /* (-1)^k factor */
        INT n = abs(-k + p);             /* |k+p| */
        REAL8 Lnp = get_Laplace_from_LaplaceCache(n, lc);
        REAL8 DLnp = get_DLaplace_from_DLaplaceCache(n, dlc);

        REAL8 term = (k & 1 ? -1.0 : 1.0) * bess * (DLnp + logprefLc*Lnp);
        // print_debug("term_K(-%d) = %.16ee / %.16e\n", k, term, sum_m);
        sum_m += term;
        if (fabs(term) <= rtol * fabs(sum_m)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    STRUCTFREE(jc, BesselJCache);
    STRUCTFREE(lc, LaplaceCache);
    STRUCTFREE(dlc, DLaplaceCache);
    REAL8 factor = pow(1.0 + beta2, a);
    return -factor * (sum_p + sum_m);   /* already includes 1/(2π) via Laplace */
}

// REAL8 K_pqa0_series_cache(INT p, INT q, INT a, REAL8 e, 
//     BesselJCache2D *bc, LaplaceCache2D *lc, DLaplaceCache2D *dlc, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     // if (a==0) return jn(p,q*e);  
//     // if(a==0) return jn(p,q*e);
//     REAL8 beta = eval_beta(e);
//     if (q==0 && a==0) return p==0 ? log((1.+sqrt(1.-e*e))/2.) : -pow(beta, abs(p)) / abs(p);
//     REAL8 beta2 = beta * beta;
//     REAL8 z    = q * e;

//     /* choose truncation N so that tail < tol                          */
//     INT p_abs = abs(p);
//     INT N = choose_N_a(p_abs, q, e, beta, atol, a);

//     INT kmin = -N, kmax = N;
//     INT M = kmax - kmin + 1;

//     /* Bessel array J_k(z) (both signs)                                */
//     // REAL8 *Jk = malloc(M * sizeof(REAL8));
//     // bessel_fill(kmin, kmax, z, Jk);

//     size_t nMax = abs(kmax) + abs(p); // just guess
//     // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     // DLaplaceCache *dlc = CreateDLaplaceCache(nMax, a, beta);
//     REAL8 logprefLc = log(1 + beta2);

//     /* series summation                                                */
//     REAL8 sum = 0.0;
//     for (INT k = kmin; k <= kmax; ++k) {
//         INT idx = k - kmin;
//         // REAL8 bess = Jk[idx];
//         REAL8 bess = get_BesselJ_from_BesselJCache2D(k, q, bc);
//         if (k & 1) bess = -bess;          /* (-1)^k factor */
//         INT n = abs(k + p);             /* |k+p| */
//         // REAL8 Lnp = get_Laplace_from_LaplaceCache(n, lc);
//         // REAL8 DLnp = get_DLaplace_from_DLaplaceCache(n, dlc);
//         REAL8 Lnp = get_Laplace_from_LaplaceCache2D(n, a, lc);
//         REAL8 DLnp = get_DLaplace_from_DLaplaceCache2D(n, a, dlc);

//         REAL8 term = bess * (DLnp + logprefLc*Lnp);
//         // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
//         sum += term;
//         if (fabs(term) < rtol * fabs(sum)) {   /* early exit */
//             if (k > p_abs) break;            /* tail symmetric */
//         }
//     }
//     // free(Jk);
//     // DestroyDLaplaceCache(dlc);
//     REAL8 factor = pow(1.0 + beta2, a);
//     return -factor * sum;   /* already includes 1/(2π) via Laplace */
// }

REAL8 K_pqa0_series_cache(INT p, INT q, INT a, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, DLaplaceCache2D *dlc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }

    // if (a==0) return jn(p,q*e);  
    // if(a==0) return jn(p,q*e);
    REAL8 beta = eval_beta(e);
    if (q==0 && a==0) return p==0 ? log((1.+sqrt(1.-e*e))/2.) : -pow(beta, abs(p)) / abs(p);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT p_abs = abs(p);
    INT N = choose_N_a(p_abs, q, e, beta, atol, a);

    INT kmin = -N, kmax = N;
    INT M = kmax - kmin + 1;

    /* Bessel array J_k(z) (both signs)                                */
    // REAL8 *Jk = malloc(M * sizeof(REAL8));
    // bessel_fill(kmin, kmax, z, Jk);
    // BesselJCache *jc = CreateBesselJCache(N, z);
    size_t nMax = abs(kmax) + abs(p); // just guess
    // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    // DLaplaceCache *dlc = CreateDLaplaceCache(nMax, a, beta);
    REAL8 logprefLc = log(1 + beta2);

    /* series summation         */
    REAL8 sum_p = 0.0;
    for (INT k = 0; ; ++k) {
        // INT idx = k - kmin;
        // REAL8 bess = Jk[idx];
        REAL8 bess = get_BesselJ_from_BesselJCache2D(k, q, bc);
        // if (k & 1) bess = -bess;          /* (-1)^k factor */
        INT n = abs(k + p);             /* |k+p| */
        REAL8 Lnp = get_Laplace_from_LaplaceCache2D(n, a, lc);
        REAL8 DLnp = get_DLaplace_from_DLaplaceCache2D(n, a, dlc);

        REAL8 term = (k & 1 ? -1.0 : 1.0) * bess * (DLnp + logprefLc*Lnp);
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sum_p += term;
        if (fabs(term) <= rtol * fabs(sum_p)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    REAL8 sum_m = 0.0;
    // negative summation
    for (INT k = 1; ; ++k) {
        // INT idx = k - kmin;
        // REAL8 bess = Jk[idx];
        REAL8 bess = get_BesselJ_from_BesselJCache2D(-k, q, bc);
        // if (k & 1) bess = -bess;          /* (-1)^k factor */
        INT n = abs(-k + p);             /* |k+p| */
        REAL8 Lnp = get_Laplace_from_LaplaceCache2D(n, a, lc);
        REAL8 DLnp = get_DLaplace_from_DLaplaceCache2D(n, a, dlc);

        REAL8 term = (k & 1 ? -1.0 : 1.0) * bess * (DLnp + logprefLc*Lnp);
        // print_debug("term_K(-%d) = %.16ee / %.16e\n", k, term, sum_m);
        sum_m += term;
        if (fabs(term) <= rtol * fabs(sum_m)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    // STRUCTFREE(jc, BesselJCache);
    // STRUCTFREE(lc, LaplaceCache);
    // STRUCTFREE(dlc, DLaplaceCache);
    REAL8 factor = pow(1.0 + beta2, a);
    return -factor * (sum_p + sum_m);   /* already includes 1/(2π) via Laplace */
}

// REAL8 Kpa0_Approx(int p, int a, REAL8 e)
// {
//     if(a <= 0 || a>14)
//     {
//         XPrintError("Error - %s: a = %d must be in [1,14]\n", __func__, a);
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     if(e < 0.0 || e >= 1.0)
//     {
//         XPrintError("Error - %s: e = %f must be in [0,1)\n", __func__, e);
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     int absp = abs(p);
//     if (absp < 5 || absp > 500) {
//         XPrintError("Error - %s: |p| = %d must > 5 and < 500\n", __func__, absp);
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     REAL8 e2 = e*e;
//     REAL8 de2 = 1. - e2;
//     REAL8 de = sqrt(de2);
//     REAL8 lnde = log(de);
//     REAL8 lnPref = absp*log(e) - (2.*a-1.)*lnde;
//     size_t indp, inda;
//     indp = absp - 5;
//     inda = a-1;
//     REAL8 mKRightLog0 = COEFF_mkRln0[inda]*lnde;
//     REAL8 lnmKregTaylor = COEFF_lnrL0[indp][inda] + e2*(
//         COEFF_lnrL1[indp][inda] + e2*(
//             COEFF_lnrL2[indp][inda] + e2*(
//                 COEFF_lnrL3[indp][inda] + e2*(
//                     COEFF_lnrL4[indp][inda] + e2*(
//                         COEFF_lnrL5[indp][inda] + e2*COEFF_lnrL6[indp][inda]
//                     )
//                 )
//             )
//         )
//     );
//     REAL8 mKCorr;
//     if (absp > 100 && e < 0.48) mKCorr = 1.;
//     else {
//         REAL8 termU, termV, termln;
//         termU = COEFF_csy0[indp][inda] + e2*(
//             COEFF_csy1[indp][inda] + e2*(
//                 COEFF_csy2[indp][inda] + e2* COEFF_csy3[indp][inda]
//             )
//         );
//         termV = COEFF_csz0[indp][inda] + e2*(
//             COEFF_csz1[indp][inda] + e2*COEFF_csz2[indp][inda]
//         );
//         termln = lnde*de2*(COEFF_cv0[indp][inda] + e2*(
//             COEFF_cv1[indp][inda] + e2*COEFF_cv2[indp][inda]
//         ));
//         REAL8 wind1, wind2, wind3;
//         wind1 = exp(-COEFF_q1[indp][inda]*de2);
//         wind2 = exp(-COEFF_q2[indp][inda]*de2);
//         wind3 = exp(-COEFF_q3[indp][inda]*de2);
//         mKCorr = 1. + e2*e2*e2*e2*e2*e2*e2*(wind1*termU + wind2*de*termV + wind3*termln + de2*de2*de2*de*(COEFF_lbd1[indp][inda] + COEFF_lbd2[indp][inda]*lnde)) / (1. + e2);
//     }
//     // return exp(lnPref + lnmKregTaylor) * mKCorr + mKRightLog0;
//     return mKCorr;
//     // return -(exp(lnmKregTaylor + lnPref) * mKCorr + exp(lnPref) * mKRightLog0);
//     // return -exp(lnmKregTaylor) * mKCorr;
// }

/* ---------- test ------------------------------------------ */
// INT testit()
// {
//     INT    p_min = 0;
//     INT    p_max = 5;
//     INT    q     = 1;
//     INT    a     = 20; 
//     REAL8  e     = 0.99;

//     int P = p_max - p_min + 1;
//     double *vals = malloc(sizeof(double)*P);

//     Jpqa_rangep(p_min, p_max, q, a, e, 1e-16, vals);

//     for (int i=0;i<P;++i)
//         printf("Jnew[%d,%d,%d,%.6g]/(2π) = %.15e\n",
//                p_min+i, q, a, e, vals[i]);

//     free(vals);
//     return 0;
// }

