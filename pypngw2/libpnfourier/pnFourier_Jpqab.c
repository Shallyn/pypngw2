/**
* Writer: Xiaolin.liu
* shallyn.liu@foxmail.com
**/
#include "pnFourier_Jpqab.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define KMIN 5
#define NMIN 5
/* ------------- Utils ------------- */
static inline REAL8 safe_pow(INT m, REAL8 x)
{
    REAL8 delta = 1.0 - x;                     
    if (delta > 1e-8)       
        return pow(x, m);
    return exp(m * log1p(-delta)); 
}

/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                            J^1_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */



/* --------- pre computing Bessel J_k(x) --------- */
static void precompute_bessel(REAL8 *J, INT kMax, REAL8 x)
{
    for (INT k = 0; k <= kMax; ++k) J[k] = jn(k, x);
}

/* -------------- Kahan summation -------------- */
static inline void kahan_add(REAL8 *sum, REAL8 *c, REAL8 term)
{
    REAL8 y = term - *c;
    REAL8 t = *sum + y;
    *c = (t - *sum) - y;
    *sum = t;
}


/* ----------  main  ------------ */
REAL8 J_0q01_series(INT q, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT min_eord = 2;
    INT kMax = (INT)(log(atol) / log(e)) - min_eord;
    kMax = (kMax < KMIN ? KMIN : kMax);
    BesselJCache *jc = CreateBesselJCache(kMax, q*e);
    REAL8 sum = 0.0;
    REAL8 betaPow = beta;
    for (INT k = 0;  ; ++k) {
        REAL8 bess = get_BesselJ_from_BesselJCache(2*k+1, jc);
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = bess * betaPow / (REAL8)(2*k+1);
        betaPow *= beta2;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sum += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sum)) {   /* early exit */
            if (k > min_eord) break;            /* tail symmetric */
        }
    }
    STRUCTFREE(jc, BesselJCache);
    return 2.*sum;
}

REAL8 J_pq01_series(INT p, INT q, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    // if (a==0) return jn(p,q*e);  
    REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT p_abs = abs(p);
    INT min_eord = (p_abs==0 ? 2 : p_abs);
    INT kMax = (INT)(log(atol) / log(e)) - min_eord;
    kMax = (kMax < KMIN ? KMIN : kMax);

    /* Bessel array J_k(z) (both signs)                                */
    // REAL8 *Jk = MYMalloc(M * sizeof(REAL8));
    // bessel_fill(kmin, kmax, z, Jk);
    BesselJCache *jc = CreateBesselJCache(kMax, q*e);

    /* series summation                                                */
    REAL8 sumP = 0.0;
    REAL8 betaPow = beta;
    for (INT k = 1;  ; ++k) {
        REAL8 bessp = get_BesselJ_from_BesselJCache(k+p, jc);
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = bessp * betaPow / (REAL8)k;
        betaPow *= beta;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sumP += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sumP)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    REAL8 sumM = 0.0;
    betaPow = beta;
    for (INT k = 1;  ; ++k) {
        // REAL8 bessp = get_BesselJ_from_BesselJCache(k+p, jc);
        REAL8 bessm = get_BesselJ_from_BesselJCache(p-k, jc);
        // REAL8 diffbess = bessp - bessm;
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = -bessm * betaPow / (REAL8)k;
        betaPow *= beta;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sumM += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sumM)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }
    // MYFree(Jk);
    STRUCTFREE(jc, BesselJCache);
    return sumP + sumM;   /* already includes 1/(2π) via Laplace */
}

// REAL8 J_pqa1_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }

//     // const REAL8 tol = 1e-16;
//     if (p==0 && q==0) return 0.0;
//     if (a==0 && p==0) return J_0q01_series(q, e, atol, rtol);
//     if (a==0) return J_pq01_series(p, q, e, atol, rtol);
//     // REAL8 retsgn = 1.0;
//     // if (p > 0) { // we find p < 0 is more accurate than p > 0
//     //     p = -p;
//     //     q = -q;
//     //     retsgn = -1.0;
//     // }
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     INT kMax = (INT)(fabs(xArg) + 8.0*sqrt(fabs(xArg)+1.0) + 25.0);
//     INT kMin = KMIN;
//     BesselJCache *bc = CreateBesselJCache(kMax, xArg);
//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);

//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                   /* 再多给 10%+4 步保险 */
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;         /* Kahan 累加器 */
//     REAL8 betaPow = 1.0;
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     REAL8 tol_k = atol*0.1;
//     {
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         INT n, k;
//         for ( n = 0;  ; ++n)
//         {
//             REAL8 In = 0.0, cIn = 0.0;    /* 内层和 (d,k) */

//             for (INT dIdx = 0; dIdx < 2; ++dIdx)
//             {
//                 INT d = dIdx ? 1 : -1;
//                 INT m0 = p + d*(n+1);

//                 REAL8 sum_k = 0.0, ck = 0.0;
//                 for ( k = 0;  ; ++k)
//                 {
//                     /* 正 k ------------------ */
//                     // REAL8 termP = ((k & 1) ? -Jtbl[k] : Jtbl[k]) * laplace_na(m0 + k, a, beta);
//                     // REAL8 termP = ((k & 1) ? -Jtbl[k] : Jtbl[k]) * get_Laplace_from_LaplaceCache(m0 + k, lc);
//                     REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * get_Laplace_from_LaplaceCache(m0 + k, lc);
//                     kahan_add(&sum_k, &ck, termP);

//                     /* 负 k (k≠0) ------------ */
//                     if (k) {
//                         // REAL8 termN =  Jtbl[k] * laplace_na(m0 - k, a, beta);  
//                         REAL8 termN =  get_BesselJ_from_BesselJCache(k, bc) * get_Laplace_from_LaplaceCache(m0 - k, lc);
//                         kahan_add(&sum_k, &ck, termN);

//                     }

//                     REAL8 rk = lambda / (k + 1.0);          /* 近似比值 */
//                     if (rk < 0.2 && k>kMin) {
//                         REAL8 tail_est = rk / (1.0 - rk) * fabs(termP);
//                         if (tail_est <= tol_k * fabs(sum_k)) break;
//                     }
//                 } /* k loop */

//                 kahan_add(&In, &cIn, sum_k*d);
//             } /* d loop */

//             REAL8 coeff = betaPow / (REAL8)(n+1.0);
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* 预留给下一个 n */
//             REAL8 r = beta * (n+1.0) / (n+2.0);
//             if (fabs(term) * r <= rtol * (1.0 - r) && n>nMin) {
//                 break;
//             }

//         } /* n loop */

//         kahan_add(&Sn, &corrSn, Sn_loc);
//     } /* end parallel */
//     //free(Jtbl);
//     DestroyBesselJCache(bc);
//     DestroyLaplaceCache(lc);
//     return pref * Sn;
// }

REAL8 J_pqa1_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }

    // const REAL8 tol = 1e-16;
    if (p==0 && q==0) return 0.0;
    if (a==0 && p==0) return J_0q01_series(q, e, atol, rtol);
    if (a==0) return J_pq01_series(p, q, e, atol, rtol);
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    INT kMax = (INT)(fabs(xArg) + 8.0*sqrt(fabs(xArg)+1.0) + 25.0);
    INT kMin = KMIN;
    BesselJCache *bc = CreateBesselJCache(kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);

    /* ------------ init ------------ */
    INT ksum_cum = 0;
    INT ksum_cum_stop = 2;
    REAL8 tol_L = atol;
    REAL8 sum_p, sum_m;
    sum_p = 0.0;
    for (INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta;
        for (INT n = 1; ; n++)
        {
            REAL8 Lterm =  betaPow*(get_Laplace_from_LaplaceCache(m0+n, lc) - get_Laplace_from_LaplaceCache(m0-n, lc))/ (REAL8)n;
            sumL += Lterm;
            betaPow *= beta;
            if (fabs(Lterm) <= tol_L * fabs(sumL)) {
                if (n>2) break;
            }
        } /* n loop p */
        // REAL8 sumL_m = 0.0;
        // betaPow = beta;
        // for (INT n = 1; ; n++)
        // {
        //     REAL8 Lterm =  -betaPow*(get_Laplace_from_LaplaceCache(m0-n, lc))/ (REAL8)n;
        //     sumL_m += Lterm;
        //     betaPow *= beta;
        //     if (fabs(Lterm) <= tol_L * fabs(sumL_m)) {
        //         if (n>2) break;
        //     }
        // } /* n loop p */
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            ksum_cum++;
            // print_debug("termL[%d|%d] = %.16e\n", k, ksum_cum, sumL_p + sumL_m);
            if (k>2 && ksum_cum > ksum_cum_stop) break;
        } else
            ksum_cum = 0;
    } /* k loop */

    sum_m = 0;
    ksum_cum = 0;
    for (INT k = 1;  ; k++)
    {
        INT m0 = p-k;

        REAL8 sumL = 0.0;
        REAL8 betaPow = beta;
        for (INT n = 1; ; n++)
        {
            REAL8 Lterm =  betaPow*(get_Laplace_from_LaplaceCache(m0+n, lc) - get_Laplace_from_LaplaceCache(m0-n, lc))/ (REAL8)n;
            sumL += Lterm;
            betaPow *= beta;
            if (fabs(Lterm) <= tol_L * fabs(sumL)) {
                if (n>2) break;
            }
        } /* n loop */

        // REAL8 sumL_m = 0.0;
        // betaPow = beta;
        // for (INT n = 1; ; n++)
        // {
        //     REAL8 Lterm =  -betaPow*(get_Laplace_from_LaplaceCache(m0-n, lc))/ (REAL8)n;
        //     sumL_m += Lterm;
        //     betaPow *= beta;
        //     if (fabs(Lterm) <= tol_L * fabs(sumL_m)) {
        //         if (n>2) break;
        //     }
        // } /* n loop */

        // print_debug("termL[-%d] = %.16e\n", k, sumL_p + sumL_m);
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(-k, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            ksum_cum++;
            if (k>2 && ksum_cum > ksum_cum_stop) break;
        } else 
            ksum_cum = 0;
    } /* k loop */

    //free(Jtbl);
    STRUCTFREE(bc, BesselJCache);
    STRUCTFREE(lc, LaplaceCache);
    return pref * (sum_p + sum_m);
}

// with cache
REAL8 J_0q01_series_cache(INT q, REAL8 e, BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT min_eord = 2;
    INT kMax = (INT)(log(atol) / log(e)) - min_eord;
    kMax = (kMax < KMIN ? KMIN : kMax);
    // BesselJCache *jc = CreateBesselJCache(kMax, q*e);
    REAL8 sum = 0.0;
    REAL8 betaPow = beta;
    for (INT k = 0;  ; ++k) {
        REAL8 bess = get_BesselJ_from_BesselJCache2D(2*k+1, q, bc);
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = bess * betaPow / (REAL8)(2*k+1);
        betaPow *= beta2;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sum += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sum)) {   /* early exit */
            if (k > min_eord) break;            /* tail symmetric */
        }
    }
    // STRUCTFREE(jc, BesselJCache);
    return 2.*sum;
}

REAL8 J_pq01_series_cache(INT p, INT q, REAL8 e, BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    // if (a==0) return jn(p,q*e);  
    REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta * beta;
    REAL8 z    = q * e;

    /* choose truncation N so that tail < tol                          */
    INT p_abs = abs(p);
    INT min_eord = (p_abs==0 ? 2 : p_abs);
    INT kMax = (INT)(log(atol) / log(e)) - min_eord;
    kMax = (kMax < KMIN ? KMIN : kMax);

    /* Bessel array J_k(z) (both signs)                                */
    // REAL8 *Jk = MYMalloc(M * sizeof(REAL8));
    // bessel_fill(kmin, kmax, z, Jk);
    // BesselJCache *jc = CreateBesselJCache(kMax, q*e);

    /* series summation                                                */
    REAL8 sumP = 0.0;
    REAL8 betaPow = beta;
    for (INT k = 1;  ; ++k) {
        REAL8 bessp = get_BesselJ_from_BesselJCache2D(k+p, q, bc);
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = bessp * betaPow / (REAL8)k;
        betaPow *= beta;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sumP += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sumP)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }

    REAL8 sumM = 0.0;
    betaPow = beta;
    for (INT k = 1;  ; ++k) {
        // REAL8 bessp = get_BesselJ_from_BesselJCache(k+p, jc);
        REAL8 bessm = get_BesselJ_from_BesselJCache2D(p-k, q, bc);
        // REAL8 diffbess = bessp - bessm;
        // REAL8 Lnp = laplace_na(n, a, beta);
        REAL8 term = -bessm * betaPow / (REAL8)k;
        betaPow *= beta;
        // print_debug("Laplace(%d, %d, %e) = %e\n", n, a, beta, Lnp);
        sumM += term;
        // print_debug("[%d]diffbess = %.16e, beta^k/k = %.16e, term = %.16e\n", k, diffbess, betaPow / (REAL8)k, term);
        if (fabs(term) <= rtol * fabs(sumM)) {   /* early exit */
            if (k > p_abs) break;            /* tail symmetric */
        }
    }
    // MYFree(Jk);
    // STRUCTFREE(jc, BesselJCache);
    return sumP + sumM;   /* already includes 1/(2π) via Laplace */
}

// REAL8 J_pqa1_series_cache(INT p, INT q, INT a, REAL8 e, 
//     BesselJCache2D *bc, LaplaceCache2D *lc ,REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (p==0 && q==0) return 0.0;
//     REAL8 retsgn = 1.0;
//     if (p > 0) { // we find p < 0 is more accurate than p > 0
//         p = -p;
//         q = -q;
//         retsgn = -1.0;
//     }
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = retsgn * beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     INT kMax = (INT)(fabs(xArg) + 8.0*sqrt(fabs(xArg)+1.0) + 25.0);
//     INT kMin = KMIN;
//     // BesselJCache *tmpbc = CreateBesselJCache(kMax, xArg);
//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     // LaplaceCache *tmplc = CreateLaplaceCache(nMax, a, beta);

//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                   /* 再多给 10%+4 步保险 */
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;         /* Kahan 累加器 */
//     REAL8 betaPow = 1.0;
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     REAL8 tol_k = atol*0.1;
//     /* ------- 并行 n 层求和 (OpenMP) ------- */
//     {
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         INT n, k;
//         for ( n = 0;  ; ++n)
//         {
//             REAL8 In = 0.0, cIn = 0.0;    /* 内层和 (d,k) */

//             for (INT dIdx = 0; dIdx < 2; ++dIdx)
//             {
//                 INT d = dIdx ? 1 : -1;
//                 INT m0 = p + d*(n+1);

//                 REAL8 sum_k = 0.0, ck = 0.0;
//                 for ( k = 0;  ; ++k)
//                 {
//                     /* 正 k ------------------ */
//                     // REAL8 termP = ((k & 1) ? -Jtbl[k] : Jtbl[k]) * laplace_na(m0 + k, a, beta);
//                     // REAL8 termP = ((k & 1) ? -Jtbl[k] : Jtbl[k]) * get_Laplace_from_LaplaceCache(m0 + k, lc);
//                     REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * get_Laplace_from_LaplaceCache2D(m0 + k, a, lc);
//                     // REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, tmpbc) * get_Laplace_from_LaplaceCache(m0 + k, tmplc);
//                     kahan_add(&sum_k, &ck, termP);

//                     /* 负 k (k≠0) ------------ */
//                     if (k) {
//                         // REAL8 termN =  Jtbl[k] * laplace_na(m0 - k, a, beta);  
//                         REAL8 termN =  get_BesselJ_from_BesselJCache2D(k, q, bc) * get_Laplace_from_LaplaceCache2D(m0 - k, a, lc);
//                         // REAL8 termN =  get_BesselJ_from_BesselJCache(k, tmpbc) * get_Laplace_from_LaplaceCache(m0 - k, tmplc);
//                         kahan_add(&sum_k, &ck, termN);

//                     }

//                     REAL8 rk = lambda / (k + 1.0);          /* 近似比值 */
//                     if (rk < 0.2 && k>kMin) {
//                         REAL8 tail_est = rk / (1.0 - rk) * fabs(termP);
//                         if (tail_est <= tol_k * fabs(sum_k)) break;
//                     }
//                 } /* k loop */

//                 kahan_add(&In, &cIn, sum_k*d);
//             } /* d loop */

//             REAL8 coeff = betaPow / (REAL8)(n+1.0);
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* 预留给下一个 n */
//             REAL8 r = beta * (n+1.0) / (n+2.0);
//             if (fabs(term) * r <= rtol * (1.0 - r) && n>nMin) {
//                 break;
//             }

//         } /* n loop */

//         kahan_add(&Sn, &corrSn, Sn_loc);
//     } /* end parallel */
//     //free(Jtbl);
//     // DestroyBesselJCache(tmpbc);
//     // DestroyLaplaceCache(tmplc);
//     return pref * Sn;
// }

REAL8 J_pqa1_series_cache(INT p, INT q, INT a, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }

    // const REAL8 tol = 1e-16;
    if (p==0 && q==0) return 0.0;
    if (a==0 && p==0) return J_0q01_series_cache(q, e, bc, lc, atol, rtol);
    if (a==0) return J_pq01_series_cache(p, q, e, bc, lc, atol, rtol);
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    INT kMax = (INT)(fabs(xArg) + 8.0*sqrt(fabs(xArg)+1.0) + 25.0);
    INT kMin = KMIN;
    // BesselJCache *bc = CreateBesselJCache(kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);

    /* ------------ init ------------ */
    INT ksum_cum = 0;
    INT ksum_cum_stop = 2;
    REAL8 tol_L = atol;
    REAL8 sum_p, sum_m;
    sum_p = 0.0;
    for (INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        // plus
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta;
        for (INT n = 1; ; n++)
        {
            REAL8 Lterm =  betaPow*(get_Laplace_from_LaplaceCache2D(m0+n, a, lc) - get_Laplace_from_LaplaceCache2D(m0-n, a, lc))/ (REAL8)n;
            sumL += Lterm;
            betaPow *= beta;
            if (fabs(Lterm) <= tol_L * fabs(sumL)) {
                if (n>2) break;
            }
        } /* n loop p */
        // plus
        // REAL8 sumL_m = 0.0;
        // betaPow = beta;
        // for (INT n = 1; ; n++)
        // {
        //     REAL8 Lterm =  -betaPow*(get_Laplace_from_LaplaceCache2D(m0-n, a, lc))/ (REAL8)n;
        //     sumL_m += Lterm;
        //     betaPow *= beta;
        //     if (fabs(Lterm) <= tol_L * fabs(sumL_m)) {
        //         if (n>2) break;
        //     }
        // } /* n loop p */
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            ksum_cum++;
            // print_debug("termL[%d|%d] = %.16e\n", k, ksum_cum, sumL_p + sumL_m);
            if (k>2 && ksum_cum > ksum_cum_stop) break;
        } else
            ksum_cum = 0;
    } /* k loop */

    sum_m = 0;
    ksum_cum = 0;
    for (INT k = 1;  ; k++)
    {
        INT m0 = p-k;

        REAL8 sumL = 0.0;
        REAL8 betaPow = beta;
        for (INT n = 1; ; n++)
        {
            REAL8 Lterm =  betaPow*(get_Laplace_from_LaplaceCache2D(m0+n, a, lc) - get_Laplace_from_LaplaceCache2D(m0-n, a, lc))/ (REAL8)n;
            sumL += Lterm;
            betaPow *= beta;
            if (fabs(Lterm) <= tol_L * fabs(sumL)) {
                if (n>2) break;
            }
        } /* n loop */

        // REAL8 sumL_m = 0.0;
        // betaPow = beta;
        // for (INT n = 1; ; n++)
        // {
        //     REAL8 Lterm =  -betaPow*(get_Laplace_from_LaplaceCache2D(m0-n, a, lc))/ (REAL8)n;
        //     sumL_m += Lterm;
        //     betaPow *= beta;
        //     if (fabs(Lterm) <= tol_L * fabs(sumL_m)) {
        //         if (n>2) break;
        //     }
        // } /* n loop */

        // print_debug("termL[-%d] = %.16e\n", k, sumL_p + sumL_m);
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(-k, q, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            ksum_cum++;
            if (k>2 && ksum_cum > ksum_cum_stop) break;
        } else 
            ksum_cum = 0;
    } /* k loop */

    //free(Jtbl);
    // STRUCTFREE(bc, BesselJCache);
    // STRUCTFREE(lc, LaplaceCache);
    return pref * (sum_p + sum_m);
}

/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                            J^2_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */

/* ----------  main  ------------ */
// REAL8 J_pqa2_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (p==0 && q==0) return 0.0;
//     // const REAL8 tol = 1e-16;
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = beta*beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
//     const INT kMin = KMIN;
//     BesselJCache *bc = CreateBesselJCache(kMax, xArg);

//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;         /* Kahan summation */
//     REAL8 betaPow = 1.0;
//     REAL8 tol_k = atol*1e-1;
//     {
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         for ( int n = 0;  ; n++)
//         {
//             REAL8 In = 0.0, cIn = 0.0;   
//             for (int s=0; s<=n ; s++) {
//                 REAL8 spref = (n+2.)/((s+1)*(n-s+1));
//                 for (INT dIdx = 0; dIdx < 2; ++dIdx)
//                 {
//                     INT d = dIdx ? 1 : -1;
//                     INT shift_s = p + d*(n - 2*s);            /* Δ */
//                     INT shift_n = p + d*(n + 2);              /* Δ */

//                     REAL8 sum_k = 0.0, ck = 0.0;
//                     for ( int k = 0;  ; ++k)
//                     {
//                         /* 正 k ------------------ */
//                         REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * 
//                             (get_Laplace_from_LaplaceCache(k+shift_s, lc) - 
//                              get_Laplace_from_LaplaceCache(k+shift_n, lc));
//                         REAL8 termN = 0.0;

//                         /* 负 k (k≠0) ------------ */
//                         if (k) {
//                             termN =  get_BesselJ_from_BesselJCache(k, bc)
//                                 * (get_Laplace_from_LaplaceCache(-k+shift_s, lc) - 
//                              get_Laplace_from_LaplaceCache(-k+shift_n, lc));   
                            
//                         }
//                         kahan_add(&sum_k, &ck, spref*(termN + termP));
//                         REAL8 rk = lambda * log(n+2) / (k + 1.0);          /* approx */
//                         if (rk < 0.2 && k>kMin) {
//                             REAL8 tail_est = rk / (1.0 - rk) * (fabs(termP) + fabs(termN));
//                             if (tail_est <= tol_k * fabs(sum_k)) break;
//                         }
//                     } /* k loop */
//                     kahan_add(&In, &cIn, sum_k);
//                 } /* d loop */
//             } /* s loop */
//             REAL8 coeff = betaPow / (n+2.);
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* for next n */
//             REAL8 r = beta * (n+1.0) / (n+3.0);
//             if (fabs(term) * r <= rtol * (1.0 - r) && n>nMin) {
//                 break;
//             }

//         } /* n loop */
//         kahan_add(&Sn, &corrSn, Sn_loc);
//     } /* end parallel */
//     DestroyBesselJCache(bc);
//     DestroyLaplaceCache(lc);
//     return -pref * Sn;
// }

REAL8 J_pq02_series(INT p, INT q, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    const REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta*beta;
    const REAL8 xArg = q * e;

    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    REAL8 sum_p = 0.0;
    REAL8 betaPow = beta2;
    INT ncum = 0, ncumMax = 2;
    INT nMin = 2;
    for (INT n=0;; n++) {
        REAL8 term = 0.0;
        for (INT s=0; s<=n; s++) {
            term += (get_BesselJ_from_BesselJCache(p+n-2*s, bc) + get_BesselJ_from_BesselJCache(p-n+2*s, bc)
             - get_BesselJ_from_BesselJCache(p+n+2, bc) - get_BesselJ_from_BesselJCache(p-n-2, bc)) / ((s+1.)*(n-s+1.));
        }
        term *= betaPow;
        sum_p += term;
        betaPow *= beta;
        // print_debug("term[%d] = %.16e, sum = %.16e\n", n, term, sum_p);
        if (fabs(term) <= rtol * fabs(sum_p)) {
            ncum++;
            if (n > nMin && ncum > ncumMax) break;
        } else
            ncum = 0;
    }
    STRUCTFREE(bc, BesselJCache);
    return -(sum_p);
}

REAL8 J_pqa2_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }
    if (p==0 && q==0) return 0.0;
    if (a==0) return J_pq02_series(p, q, e, atol, rtol);
    // const REAL8 tol = 1e-16;
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    const INT kMin = KMIN;
    BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
    // if (!Jtbl) { 
    //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
    //     X_ERROR_REAL8(X_ENOMEM);
    // }
    // precompute_bessel(Jtbl, kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    /* ------------ n upper ------------ */
    INT kcum_max = 2;
    REAL8 beta2 = beta*beta;
    REAL8 tolL = atol*0.1;
    REAL8 sum_p = 0.0, sum_m = 0.0;
    INT kcum = 0;
    for ( INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta2;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s=0; s<=n ; s++) {
                term += (
                    get_Laplace_from_LaplaceCache(m0 + (n-2*s), lc) + 
                    get_Laplace_from_LaplaceCache(m0 - (n-2*s), lc) - 
                    get_Laplace_from_LaplaceCache(m0 + (n+2), lc) - 
                    get_Laplace_from_LaplaceCache(m0 - (n+2), lc) 
                ) / ((s+1.)*(n-s+1.));
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            kcum++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */

    // negative k
    kcum = 0;
    for ( INT k = 1;  ; k++)
    {
        INT m0 = -k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta2;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s=0; s<=n ; s++) {
                term += (
                    get_Laplace_from_LaplaceCache(m0 + (n-2*s), lc) + 
                    get_Laplace_from_LaplaceCache(m0 - (n-2*s), lc) - 
                    get_Laplace_from_LaplaceCache(m0 + (n+2), lc) - 
                    get_Laplace_from_LaplaceCache(m0 - (n+2), lc) 
                ) / ((s+1.)*(n-s+1.));
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(-k, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            kcum ++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */


    STRUCTFREE(bc, BesselJCache);
    STRUCTFREE(lc, LaplaceCache);
    return -pref * (sum_m + sum_p);
}

// with cache
REAL8 J_pq02_series_cache(INT p, INT q, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    const REAL8 beta = eval_beta(e);
    REAL8 beta2 = beta*beta;
    const REAL8 xArg = q * e;

    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    REAL8 sum_p = 0.0;
    REAL8 betaPow = beta2;
    INT ncum = 0, ncumMax = 2;
    INT nMin = 2;
    for (INT n=0;; n++) {
        REAL8 term = 0.0;
        for (INT s=0; s<=n; s++) {
            term += (get_BesselJ_from_BesselJCache2D(p+n-2*s, q, bc) + get_BesselJ_from_BesselJCache2D(p-n+2*s, q, bc)
             - get_BesselJ_from_BesselJCache2D(p+n+2, q, bc) - get_BesselJ_from_BesselJCache2D(p-n-2, q, bc)) / ((s+1.)*(n-s+1.));
        }
        term *= betaPow;
        sum_p += term;
        betaPow *= beta;
        // print_debug("term[%d] = %.16e, sum = %.16e\n", n, term, sum_p);
        if (fabs(term) <= rtol * fabs(sum_p)) {
            ncum++;
            if (n > nMin && ncum > ncumMax) break;
        } else
            ncum = 0;
    }
    // STRUCTFREE(bc, BesselJCache);
    return -(sum_p);
}

REAL8 J_pqa2_series_cache(INT p, INT q, INT a, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }
    if (p==0 && q==0) return 0.0;
    if (a==0) return J_pq02_series_cache(p, q, e, bc, lc, atol, rtol);
    // const REAL8 tol = 1e-16;
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    const INT kMin = KMIN;
    // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
    // if (!Jtbl) { 
    //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
    //     X_ERROR_REAL8(X_ENOMEM);
    // }
    // precompute_bessel(Jtbl, kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    /* ------------ n upper ------------ */
    INT kcum_max = 2;
    REAL8 beta2 = beta*beta;
    REAL8 tolL = atol*0.1;
    REAL8 sum_p = 0.0, sum_m = 0.0;
    INT kcum = 0;
    for ( INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta2;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s=0; s<=n ; s++) {
                term += (
                    get_Laplace_from_LaplaceCache2D(m0 + (n-2*s), a, lc) + 
                    get_Laplace_from_LaplaceCache2D(m0 - (n-2*s), a, lc) - 
                    get_Laplace_from_LaplaceCache2D(m0 + (n+2), a, lc) - 
                    get_Laplace_from_LaplaceCache2D(m0 - (n+2), a, lc) 
                ) / ((s+1.)*(n-s+1.));
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            kcum++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */

    // negative k
    kcum = 0;
    for ( INT k = 1;  ; k++)
    {
        INT m0 = -k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta2;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s=0; s<=n ; s++) {
                term += (
                    get_Laplace_from_LaplaceCache2D(m0 + (n-2*s), a, lc) + 
                    get_Laplace_from_LaplaceCache2D(m0 - (n-2*s), a, lc) - 
                    get_Laplace_from_LaplaceCache2D(m0 + (n+2), a, lc) - 
                    get_Laplace_from_LaplaceCache2D(m0 - (n+2), a, lc) 
                ) / ((s+1.)*(n-s+1.));
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(-k, q, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            kcum ++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */
    // STRUCTFREE(bc, BesselJCache);
    // STRUCTFREE(lc, LaplaceCache);
    return -pref * (sum_m + sum_p);
}

// REAL8 J_pqa2_series_cache(INT p, INT q, INT a, REAL8 e, 
//     BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (p==0 && q==0) return 0.0;
//     // const REAL8 tol = 1e-16;
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = beta*beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/atol)) + 3.0);
//     const INT kMin = KMIN;
//     // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                   /* 再多给 10%+4 步保险 */
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;         /* Kahan 累加器 */
//     REAL8 betaPow = 1.0;
//     REAL8 tol_k = atol*0.1;
//     {
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         for ( int n = 0;  ; n++)
//         {
//             REAL8 In = 0.0, cIn = 0.0;    /* 内层和 (d,k) */
//             for (int s=0; s<=n ; s++) {
//                 REAL8 spref = (n+2.)/((s+1)*(n-s+1));
//                 for (INT dIdx = 0; dIdx < 2; ++dIdx)
//                 {
//                     INT d = dIdx ? 1 : -1;
//                     INT shift_s = p + d*(n - 2*s);            /* Δ */
//                     INT shift_n = p + d*(n + 2);              /* Δ */

//                     REAL8 sum_k = 0.0, ck = 0.0;
//                     for ( int k = 0;  ; ++k)
//                     {
//                         /* 正 k ------------------ */
//                         REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * 
//                             (get_Laplace_from_LaplaceCache2D(k+shift_s, a, lc) - 
//                              get_Laplace_from_LaplaceCache2D(k+shift_n, a, lc));
//                         REAL8 termN = 0.0;

//                         /* 负 k (k≠0) ------------ */
//                         if (k) {
//                             termN =  get_BesselJ_from_BesselJCache2D(k, q, bc)
//                                 * (get_Laplace_from_LaplaceCache2D(-k+shift_s, a, lc) - 
//                              get_Laplace_from_LaplaceCache2D(-k+shift_n, a, lc));     /*   负阶总是正号相加   */
                            
//                         }
//                         kahan_add(&sum_k, &ck, spref*(termN + termP));
//                         REAL8 rk = lambda * log(n+2) / (k + 1.0);          /* 近似比值 */
//                         if (rk < 0.2 && k>kMin) {
//                             REAL8 tail_est = rk / (1.0 - rk) * (fabs(termP) + fabs(termN));
//                             if (tail_est <= tol_k * fabs(sum_k)) break;
//                         }
//                     } /* k loop */
//                     kahan_add(&In, &cIn, sum_k);
//                 } /* d loop */
//             } /* s loop */
//             REAL8 coeff = betaPow / (n+2.);
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* 预留给下一个 n */
//             REAL8 r = beta * (n+1.0) / (n+3.0);
//             if (fabs(term) * r <= rtol * (1.0 - r) && n>nMin) {
//                 break;
//             }

//         } /* n loop */
//         kahan_add(&Sn, &corrSn, Sn_loc);
//     } /* end parallel */
//     // DestroyBesselJCache(bc);
//     // DestroyLaplaceCache(lc);
//     return -pref * Sn;
// }

/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                            J^3_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */

REAL8 J_pq03_series(INT p, INT q, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    const REAL8 beta = eval_beta(e);
    REAL8 beta3 = beta*beta*beta;
    const REAL8 xArg = q * e;

    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    REAL8 sum_p = 0.0;
    REAL8 betaPow = beta3;
    INT ncum = 0, ncumMax = 2;
    INT nMin = 2;
    for (INT n=0;; n++) {
        REAL8 term = 0.0;
        for (INT s1=0; s1<=n; s1++) {
            for (INT s2=0; s2<=s1; s2++){
                term += (
                    get_BesselJ_from_BesselJCache(p+(1+n-2*s2), bc) - 
                    get_BesselJ_from_BesselJCache(p-(1+n-2*s2), bc) +
                    get_BesselJ_from_BesselJCache(p+(1+n-2*s1+2*s2), bc) -
                    get_BesselJ_from_BesselJCache(p-(1+n-2*s1+2*s2), bc) -
                    get_BesselJ_from_BesselJCache(p+(n+3), bc) + 
                    get_BesselJ_from_BesselJCache(p-(n+3), bc) - 
                    get_BesselJ_from_BesselJCache(p+(n-1-2*s1), bc) +
                    get_BesselJ_from_BesselJCache(p-(n-1-2*s1), bc)
                ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
            }
        }
        term *= betaPow;
        sum_p += term;
        betaPow *= beta;
        // print_debug("term[%d] = %.16e, sum = %.16e\n", n, term, sum_p);
        if (fabs(term) <= rtol * fabs(sum_p)) {
            ncum++;
            if (n > nMin && ncum > ncumMax) break;
        } else
            ncum = 0;
    }
    STRUCTFREE(bc, BesselJCache);
    return -(sum_p);
}

REAL8 J_pq03_series_cache(INT p, INT q, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    const REAL8 beta = eval_beta(e);
    REAL8 beta3 = beta*beta*beta;
    const REAL8 xArg = q * e;

    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    REAL8 sum_p = 0.0;
    REAL8 betaPow = beta3;
    INT ncum = 0, ncumMax = 2;
    INT nMin = 2;
    for (INT n=0;; n++) {
        REAL8 term = 0.0;
        for (INT s1=0; s1<=n; s1++) {
            for (INT s2=0; s2<=s1; s2++){
                term += (
                    get_BesselJ_from_BesselJCache2D(p+(1+n-2*s2), q, bc) - 
                    get_BesselJ_from_BesselJCache2D(p-(1+n-2*s2), q, bc) +
                    get_BesselJ_from_BesselJCache2D(p+(1+n-2*s1+2*s2),q , bc) -
                    get_BesselJ_from_BesselJCache2D(p-(1+n-2*s1+2*s2), q, bc) -
                    get_BesselJ_from_BesselJCache2D(p+(n+3), q, bc) + 
                    get_BesselJ_from_BesselJCache2D(p-(n+3), q, bc) - 
                    get_BesselJ_from_BesselJCache2D(p+(n-1-2*s1), q, bc) +
                    get_BesselJ_from_BesselJCache2D(p-(n-1-2*s1), q, bc)
                ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
            }
        }
        term *= betaPow;
        sum_p += term;
        betaPow *= beta;
        // print_debug("term[%d] = %.16e, sum = %.16e\n", n, term, sum_p);
        if (fabs(term) <= rtol * fabs(sum_p)) {
            ncum++;
            if (n > nMin && ncum > ncumMax) break;
        } else
            ncum = 0;
    }
    // STRUCTFREE(bc, BesselJCache);
    return -(sum_p);
}

/* ----------  main  ------------ */
// REAL8 J_pqa3_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (p==0 && q==0) return 0.0;
//     if (a==0) return J_pq03_series(p, q, e, atol, rtol);
//     // const REAL8 tol = 1e-16;
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = beta*beta*beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/atol)) + 3.0);
//     const INT kMin = KMIN;
//     BesselJCache *bc = CreateBesselJCache(kMax, xArg);

//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                   /* 再多给 10%+4 步保险 */
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;  
//     REAL8 betaPow = 1.0;
//     REAL8 tol_k = 0.1*atol;
//     {
//         int n, k;
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         for ( n = 0;  ; n++)
//         {
//             REAL8 In = 0.0, cIn = 0.0;  
//             REAL8 ln2n2 = pow(log(n+3.), 2); 
//             for (int s1=0; s1<=n ; s1++) {
//                 for (int s2=0; s2<=s1; s2++) {
//                     REAL8 spref = 1.0/((s1-s2+1.)*(n-s1+1.)*(s2+1.));
//                     for (INT dIdx = 0; dIdx < 2; ++dIdx)
//                     {
//                         INT d = dIdx ? 1 : -1;
//                         INT shift_p1 = p + d*(1 + n - 2*s2); 
//                         INT shift_p2 = p + d*(1 + n - 2*s1 + 2*s2); 
//                         INT shift_n1 = p + d*(n + 3);              /* Δ */
//                         INT shift_n2 = p + d*(n - 1 - 2*s1);              /* Δ */

//                         REAL8 sum_k = 0.0, ck = 0.0;
//                         for ( k = 0;  ; ++k)
//                         {
//                             /* 正 k ------------------ */
//                             REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * 
//                                 (get_Laplace_from_LaplaceCache(k+shift_p1, lc) + get_Laplace_from_LaplaceCache(k+shift_p2, lc) - 
//                                 get_Laplace_from_LaplaceCache(k+shift_n1, lc) - get_Laplace_from_LaplaceCache(k+shift_n2, lc));
//                             REAL8 termN = 0.0;
//                             /* 负 k (k≠0) ------------ */
//                             if (k) {
//                                 termN =  get_BesselJ_from_BesselJCache(k, bc)
//                                  * (get_Laplace_from_LaplaceCache(-k+shift_p1, lc) + get_Laplace_from_LaplaceCache(-k+shift_p2, lc) - 
//                                 get_Laplace_from_LaplaceCache(-k+shift_n1, lc) - get_Laplace_from_LaplaceCache(-k+shift_n2, lc));
                                
//                             }
//                             kahan_add(&sum_k, &ck, spref*(termN + termP));
//                             // print_debug("term[p=%d, q=%d, a=%d, n=%d, k=%d, s1=%d, s2=%d, d=%d][e=%e] = %.16e, %.16e\n",
//                             //     p, q, a, n, k, s1, s2, d, e, 
//                             //     spref*termN , spref*termP);
//                             REAL8 rk = lambda / (k + 1.0);          /* 近似比值 */
//                             if (rk < 0.2 && k>kMin) {
//                                 REAL8 tail_est = rk * ln2n2 / (1.0 - rk) * (fabs(termP) + fabs(termN)) / (n+3.);
//                                 if (tail_est <= tol_k * fabs(sum_k)) break;
//                             }
//                         } /* k loop */
//                         kahan_add(&In, &cIn, d*sum_k);
//                     } /* d loop */
//                 }
//             } /* s loop */
//             REAL8 coeff = betaPow;
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* 预留给下一个 n */
//             REAL8 r = beta * (n + 1.0) / (n + 3.0);
//             if (ln2n2 * fabs(term) * r  / (n + 3.0) <=  rtol * (1.0 - r) && n>nMin) {
//                 break;
//             }
//         } /* n loop */
//         kahan_add(&Sn, &corrSn, Sn_loc);
//         // print_debug("n = %d, k = %d\n", n, k);
//     } /* end parallel */
//     DestroyBesselJCache(bc);
//     DestroyLaplaceCache(lc);
//     return -pref * Sn;
// }

REAL8 J_pqa3_series(INT p, INT q, INT a, REAL8 e, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }
    if (p==0 && q==0) return 0.0;
    if (a==0) return J_pq03_series(p, q, e, atol, rtol);
    // const REAL8 tol = 1e-16;
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    const INT kMin = KMIN;
    BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
    // if (!Jtbl) { 
    //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
    //     X_ERROR_REAL8(X_ENOMEM);
    // }
    // precompute_bessel(Jtbl, kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    /* ------------ n upper ------------ */
    INT kcum_max = 2;
    REAL8 beta3 = beta*beta*beta;
    REAL8 tolL = atol*0.1;
    REAL8 sum_p = 0.0, sum_m = 0.0;
    INT kcum = 0;
    for ( INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta3;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;

            for (INT s1=0; s1<=n; s1++) {
                for (INT s2=0; s2<=s1; s2++){
                    term += (
                        get_Laplace_from_LaplaceCache(m0+(1+n-2*s2), lc) - 
                        get_Laplace_from_LaplaceCache(m0-(1+n-2*s2), lc) +
                        get_Laplace_from_LaplaceCache(m0+(1+n-2*s1+2*s2), lc) -
                        get_Laplace_from_LaplaceCache(m0-(1+n-2*s1+2*s2), lc) -
                        get_Laplace_from_LaplaceCache(m0+(n+3), lc) + 
                        get_Laplace_from_LaplaceCache(m0-(n+3), lc) - 
                        get_Laplace_from_LaplaceCache(m0+(n-1-2*s1), lc) +
                        get_Laplace_from_LaplaceCache(m0-(n-1-2*s1), lc)
                    ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
                }
            }
                
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(k, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            kcum++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */

    // negative k
    kcum = 0;
    for ( INT k = 1;  ; k++)
    {
        INT m0 = -k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta3;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s1=0; s1<=n; s1++) {
                for (INT s2=0; s2<=s1; s2++){
                    term += (
                        get_Laplace_from_LaplaceCache(m0+(1+n-2*s2), lc) - 
                        get_Laplace_from_LaplaceCache(m0-(1+n-2*s2), lc) +
                        get_Laplace_from_LaplaceCache(m0+(1+n-2*s1+2*s2), lc) -
                        get_Laplace_from_LaplaceCache(m0-(1+n-2*s1+2*s2), lc) -
                        get_Laplace_from_LaplaceCache(m0+(n+3), lc) + 
                        get_Laplace_from_LaplaceCache(m0-(n+3), lc) - 
                        get_Laplace_from_LaplaceCache(m0+(n-1-2*s1), lc) +
                        get_Laplace_from_LaplaceCache(m0-(n-1-2*s1), lc)
                    ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
                }
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache(-k, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            kcum ++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */


    STRUCTFREE(bc, BesselJCache);
    STRUCTFREE(lc, LaplaceCache);
    return -pref * (sum_m + sum_p);
}

// with cache
// REAL8 J_pqa3_series_cache(INT p, INT q, INT a, REAL8 e, 
//     BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
// {
//     if (e <= 0.0 || e >= EMAX) {                      
//         XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (a < 0) {
//         XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
//         X_ERROR_REAL8(X_EINVAL);
//     }
//     if (p==0 && q==0) return 0.0;
//     // const REAL8 tol = 1e-16;
//     const REAL8 beta = eval_beta(e);
//     const REAL8 pref = beta*beta*beta * pow(1.0 + beta*beta, a);
//     const REAL8 xArg = q * e;

//     /* ------------ k upper ------------ */
//     REAL8 lambda = beta * fabs(xArg) / 2.0;
//     const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/atol)) + 3.0);
//     const INT kMin = KMIN;
//     // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

//     // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
//     // if (!Jtbl) { 
//     //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
//     //     X_ERROR_REAL8(X_ENOMEM);
//     // }
//     // precompute_bessel(Jtbl, kMax, xArg);
//     size_t nMax = abs(p) + 16; // just guess
//     size_t nMin = NMIN;
//     // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
//     /* ------------ n upper ------------ */
//     /* nMax = ceil( log( tol·(1-β) ) / log β )，β<1 ⇒ logβ<0 */
//     // INT nMax = (INT)ceil(log(tol*(1.0-beta)) / log(beta));
//     // nMax += nMax/10 + 4;                   /* 再多给 10%+4 步保险 */
//     // if (nMax < 10) nMax = 10;

//     REAL8 Sn = 0.0, corrSn = 0.0;  
//     REAL8 betaPow = 1.0;
//     REAL8 tol_k = 0.1*atol;
//     {
//         int n, k;
//         REAL8 Sn_loc = 0.0, c_loc = 0.0;
//         for ( n = 0;  ; n++)
//         {
//             REAL8 In = 0.0, cIn = 0.0;  
//             REAL8 ln2n2 = pow(log(n+3.), 2); 
//             for (int s1=0; s1<=n ; s1++) {
//                 for (int s2=0; s2<=s1; s2++) {
//                     REAL8 spref = 1.0/((s1-s2+1.)*(n-s1+1.)*(s2+1.));
//                     for (INT dIdx = 0; dIdx < 2; ++dIdx)
//                     {
//                         INT d = dIdx ? 1 : -1;
//                         INT shift_p1 = p + d*(1 + n - 2*s2); 
//                         INT shift_p2 = p + d*(1 + n - 2*s1 + 2*s2); 
//                         INT shift_n1 = p + d*(n + 3);              /* Δ */
//                         INT shift_n2 = p + d*(n - 1 - 2*s1);              /* Δ */

//                         REAL8 sum_k = 0.0, ck = 0.0;
//                         for ( k = 0;  ; ++k)
//                         {
//                             /* 正 k ------------------ */
//                             REAL8 termP = ((k & 1) ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * 
//                                 (get_Laplace_from_LaplaceCache2D(k+shift_p1, a, lc) + get_Laplace_from_LaplaceCache2D(k+shift_p2, a, lc) - 
//                                 get_Laplace_from_LaplaceCache2D(k+shift_n1, a, lc) - get_Laplace_from_LaplaceCache2D(k+shift_n2, a, lc));
//                             REAL8 termN = 0.0;
//                             /* 负 k (k≠0) ------------ */
//                             if (k) {
//                                 termN =  get_BesselJ_from_BesselJCache2D(k, q, bc)
//                                  * (get_Laplace_from_LaplaceCache2D(-k+shift_p1, a, lc) + get_Laplace_from_LaplaceCache2D(-k+shift_p2, a, lc) - 
//                                 get_Laplace_from_LaplaceCache2D(-k+shift_n1, a, lc) - get_Laplace_from_LaplaceCache2D(-k+shift_n2, a, lc));
                                
//                             }
//                             kahan_add(&sum_k, &ck, spref*(termN + termP));
//                             // print_debug("term[p=%d, q=%d, a=%d, n=%d, k=%d, s1=%d, s2=%d, d=%d][e=%e] = %.16e, %.16e\n",
//                             //     p, q, a, n, k, s1, s2, d, e, 
//                             //     spref*termN , spref*termP);
//                             REAL8 rk = lambda / (k + 1.0);          /* 近似比值 */
//                             if (rk < 0.2 && k>kMin) {
//                                 REAL8 tail_est = rk * ln2n2 / (1.0 - rk) * (fabs(termP) + fabs(termN)) / (n+3.);
//                                 if (tail_est <= tol_k * fabs(sum_k)) break;
//                             }
//                         } /* k loop */
//                         kahan_add(&In, &cIn, d*sum_k);
//                     } /* d loop */
//                 }
//             } /* s loop */
//             REAL8 coeff = betaPow;
//             REAL8 term = coeff * In;
//             kahan_add(&Sn_loc, &c_loc, term);
//             betaPow *= beta;   /* 预留给下一个 n */
//             REAL8 r = beta * (n + 1.0) / (n + 3.0);
//             if (ln2n2 * fabs(term) * r  / (n + 3.0) <=  rtol * (1.0 - r) && n> nMin) {
//                 break;
//             }
//         } /* n loop */
//         kahan_add(&Sn, &corrSn, Sn_loc);
//         // print_debug("n = %d, k = %d\n", n, k);
//     } /* end parallel */
//     // DestroyBesselJCache(bc);
//     // DestroyLaplaceCache(lc);
//     return -pref * Sn;
// }

REAL8 J_pqa3_series_cache(INT p, INT q, INT a, REAL8 e, 
    BesselJCache2D *bc, LaplaceCache2D *lc, REAL8 atol, REAL8 rtol)
{
    if (e <= 0.0 || e >= EMAX) {                      
        XPrintError("Error - %s: e = %f must be in (0,%e)\n", __func__, e , EMAX);
        X_ERROR_REAL8(X_EINVAL);
    }
    if (a < 0) {
        XPrintError("Error - %s: power a = %d must be ≥ 0\n", __func__, a );
        X_ERROR_REAL8(X_EINVAL);
    }
    if (p==0 && q==0) return 0.0;
    if (a==0) return J_pq03_series_cache(p, q, e, bc, lc, atol, rtol);
    // const REAL8 tol = 1e-16;
    const REAL8 beta = eval_beta(e);
    const REAL8 pref = pow(1.0 + beta*beta, a);
    const REAL8 xArg = q * e;

    /* ------------ k upper ------------ */
    REAL8 lambda = beta * fabs(xArg) / 2.0;
    const INT kMax = (int)ceil(lambda + sqrt(6.0*lambda*log(1.0/(atol))) + 3.0);
    const INT kMin = KMIN;
    // BesselJCache *bc = CreateBesselJCache(kMax, xArg);

    // REAL8 *Jtbl = (REAL8*)malloc((kMax+1) * sizeof(REAL8));
    // if (!Jtbl) { 
    //     XPrintError("Error - %s: cannot allocate memory for Jtbl\n", __func__);
    //     X_ERROR_REAL8(X_ENOMEM);
    // }
    // precompute_bessel(Jtbl, kMax, xArg);
    size_t nMax = abs(p) + 16; // just guess
    size_t nMin = NMIN;
    // LaplaceCache *lc = CreateLaplaceCache(nMax, a, beta);
    /* ------------ n upper ------------ */
    INT kcum_max = 2;
    REAL8 beta3 = beta*beta*beta;
    REAL8 tolL = atol*0.1;
    REAL8 sum_p = 0.0, sum_m = 0.0;
    INT kcum = 0;
    for ( INT k = 0;  ; k++)
    {
        INT m0 = k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta3;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;

            for (INT s1=0; s1<=n; s1++) {
                for (INT s2=0; s2<=s1; s2++){
                    term += (
                        get_Laplace_from_LaplaceCache2D(m0+(1+n-2*s2), a, lc) - 
                        get_Laplace_from_LaplaceCache2D(m0-(1+n-2*s2), a, lc) +
                        get_Laplace_from_LaplaceCache2D(m0+(1+n-2*s1+2*s2), a, lc) -
                        get_Laplace_from_LaplaceCache2D(m0-(1+n-2*s1+2*s2), a, lc) -
                        get_Laplace_from_LaplaceCache2D(m0+(n+3), a, lc) + 
                        get_Laplace_from_LaplaceCache2D(m0-(n+3), a, lc) - 
                        get_Laplace_from_LaplaceCache2D(m0+(n-1-2*s1), a, lc) +
                        get_Laplace_from_LaplaceCache2D(m0-(n-1-2*s1), a, lc)
                    ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
                }
            }
                
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(k, q, bc) * (sumL);
        sum_p += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_p)) {
            kcum++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */

    // negative k
    kcum = 0;
    for ( INT k = 1;  ; k++)
    {
        INT m0 = -k+p;
        REAL8 sumL = 0.0;
        REAL8 betaPow = beta3;
        for ( INT n = 0; ; n++ ) {
            REAL8 term = 0.0;
            for (INT s1=0; s1<=n; s1++) {
                for (INT s2=0; s2<=s1; s2++){
                    term += (
                        get_Laplace_from_LaplaceCache2D(m0+(1+n-2*s2), a, lc) - 
                        get_Laplace_from_LaplaceCache2D(m0-(1+n-2*s2), a, lc) +
                        get_Laplace_from_LaplaceCache2D(m0+(1+n-2*s1+2*s2), a, lc) -
                        get_Laplace_from_LaplaceCache2D(m0-(1+n-2*s1+2*s2), a, lc) -
                        get_Laplace_from_LaplaceCache2D(m0+(n+3), a, lc) + 
                        get_Laplace_from_LaplaceCache2D(m0-(n+3), a, lc) - 
                        get_Laplace_from_LaplaceCache2D(m0+(n-1-2*s1), a, lc) +
                        get_Laplace_from_LaplaceCache2D(m0-(n-1-2*s1), a, lc)
                    ) / ((s1-s2+1.)*(n-s1+1.)*(s2+1));
                }
            }
            term *= betaPow;
            sumL += term;
            betaPow *= beta;
            if (fabs(term) <= tolL * fabs(sumL)) {
                if (n>nMin) break;
            }
        }
        REAL8 termJ = (k&1 ? -1.0 : 1.0) * get_BesselJ_from_BesselJCache2D(-k, q, bc) * (sumL);
        sum_m += termJ;
        if (fabs(termJ) <= rtol * fabs(sum_m)) {
            kcum ++;
            if (k>kMin && kcum > kcum_max) break;
        } else 
            kcum = 0;
    } /* n loop */


    // STRUCTFREE(bc, BesselJCache);
    // STRUCTFREE(lc, LaplaceCache);
    return -pref * (sum_m + sum_p);
}

/* ----------------  DEMO  ---------------- */
// #ifdef FASTJ_DEMO
// INT testit(void)
// {
//     INT    p = 2, q = 3, a = 2;
//     REAL8 e = 0.3, tol = 1e-16;

//     REAL8 val = fast_J(p,q,a,e,tol);
//     printf("J[p=%d,q=%d,a=%d](e=%g) = %.15e\n", p,q,a,e,val);

//     p = 10, q = 6, a = 8;
//     e = 0.3, tol = 1e-16;

//     val = fast_J(p,q,a,e,tol);
//     printf("J[p=%d,q=%d,a=%d](e=%g) = %.15e\n", p,q,a,e,val);

//     return 0;
// }
// // #endif

