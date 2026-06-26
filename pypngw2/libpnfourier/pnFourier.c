/**
* Writer: Xiaolin.liu
* shallyn.liu@foxmail.com
**/
#include "pnFourier.h"

static INT ExtendPNEllipticCache(size_t need_p, size_t need_q, size_t need_a, PNEllipticCache *cache)
{
    size_t new_abspMax = cache->abspMax ? cache->abspMax : 1;
    size_t new_absqMax = cache->absqMax ? cache->absqMax : 1;
    size_t new_aMax = cache->aMax ? cache->aMax : 1;
    while (new_abspMax < need_p) new_abspMax <<= 1;
    while (new_absqMax < need_q) new_absqMax <<= 1;
    while (new_aMax < need_a) new_aMax <<= 1;

    size_t newSize = (new_abspMax + 1) * (2*new_absqMax + 1) * (new_aMax + 1);
    REAL8 *newBuf = MYMalloc(newSize * sizeof(REAL8));
    if (!newBuf) {
        X_ERROR(X_ENOMEM);
    }
    for (size_t i=0; i<newSize; i++) newBuf[i] = X_REAL8_FAIL_NAN;
    for (size_t ia = 0; ia <= cache->aMax; ia++)
        for (size_t ip = 0; ip <= cache->abspMax; ip++) {
            size_t old_base = ia * cache->stride_a + ip * cache->stride_p;
            size_t new_base = ia * (new_abspMax + 1)*(2*new_absqMax + 1) + ip * (2*new_absqMax + 1);
            memcpy(&newBuf[new_base],
                &(cache->buf[old_base]),
                    (cache->stride_p)*sizeof(REAL8));
        }
    MYFree(cache->buf);
    cache->buf = newBuf;
    cache->aMax = new_aMax;
    cache->abspMax = new_abspMax;
    cache->absqMax = new_absqMax;
    cache->size = newSize;
    cache->stride_a = (cache->abspMax + 1) * (2*cache->absqMax + 1);
    cache->stride_p = 2*cache->absqMax + 1;
    return X_SUCCESS;
}

PNEllipticCache *CreatePNEllipticCache(size_t abspMax, size_t absqMax, size_t aMax, REAL8 e, REAL8 val0)
{
    PNEllipticCache *ret = MYMalloc(sizeof(PNEllipticCache));
    ret->aMax = 0;
    ret->abspMax = 0;
    ret->absqMax = 0;
    ret->stride_a = (ret->abspMax + 1) * (2*ret->absqMax + 1);
    ret->stride_p = 2*ret->absqMax + 1;
    ret->size = (ret->aMax + 1)*(ret->abspMax + 1)*(2*ret->absqMax + 1);
    ret->beta = eval_beta(e);
    ret->e = e;
    ret->buf = MYMalloc( sizeof(REAL8));
    ret->buf[0] = val0; // J[0, 0, 0] = 0

    if (ExtendPNEllipticCache(abspMax, absqMax, aMax, ret)!=X_SUCCESS) {
        X_ERROR_NULL(X_EFAILED);
    }
    return ret;
}

void DestroyPNEllipticCache(PNEllipticCache *cache)
{
    if (!cache)
        return;
    if (cache->buf) {
        MYFree(cache->buf);
        cache->buf = NULL;
    }
    MYFree(cache);
    cache = NULL;
    return;
}

/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                            J^0_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */
REAL8 get_Jpqa0_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) q = -q;
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa0_series(absp, q, a, cache->e, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}

REAL8 get_Jpqa0_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) q = -q;
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa0_series_cache(absp, q, a, cache->e, jc, lc, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
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
REAL8 get_Jpqa1_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    REAL8 pref = 1.0;
    if (p<0) {
        q = -q;
        pref = -1.0;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa1_series(absp, q, a, cache->e, atol, rtol);
        cache->buf[idx] = ret;
    }
    return pref*ret;
}

REAL8 get_Jpqa1_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    REAL8 pref = 1.0;
    if (p<0) {
        q = -q;
        pref = -1.0;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa1_series_cache(absp, q, a, cache->e, jc, lc, atol, rtol);
        cache->buf[idx] = ret;
    }
    return pref*ret;
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
REAL8 get_Jpqa2_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) {
        q = -q;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa2_series(absp, q, a, cache->e, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}

REAL8 get_Jpqa2_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) {
        q = -q;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa2_series_cache(absp, q, a, cache->e, jc, lc, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}


/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                            J^3_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */
REAL8 get_Jpqa3_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) {
        q = -q;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa3_series(absp, q, a, cache->e, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}

REAL8 get_Jpqa3_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    REAL8 pref = 1.0;
    if (p<0) {
        q = -q;
        pref = -1.0;
    }
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = J_pqa3_series_cache(absp, q, a, cache->e, jc, lc, atol, rtol);
        cache->buf[idx] = ret;
    }
    return pref*ret;
}


/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                              K_pqa                               */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */
REAL8 get_Kpqa0_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) q = -q;
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = K_pqa0_series(absp, q, a, cache->e, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}

REAL8 get_Kpqa0_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, DLaplaceCache2D *dlc,
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol)
{
    size_t absp = (p < 0 ? -p : p);
    size_t absq = (q < 0 ? -q : q);
    if (absp > cache->abspMax || absq > cache->absqMax || a > cache->aMax) {
        if (ExtendPNEllipticCache(absp, absq, a, cache) != X_SUCCESS) {
            X_ERROR_REAL8(X_ENOMEM);
        }
    }
    if (p<0) q = -q;
    size_t idx = (size_t)a * cache->stride_a + (size_t)absp * cache->stride_p + cache->absqMax + q;
    REAL8 ret = cache->buf[idx];
    if (X_IS_REAL8_FAIL_NAN(ret)) {
        ret = K_pqa0_series_cache(absp, q, a, cache->e, jc, lc, dlc, atol, rtol);
        cache->buf[idx] = ret;
    }
    return ret;
}

/* ---------------------------------------------------------------- */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/*                           Evaluator                              */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* ---------------------------------------------------------------- */
PNEllipticEvaluator* CreatePNEllipticEvaluator(REAL8 e, REAL8 atol, REAL8 rtol)
{
    PNEllipticEvaluator* ret = (PNEllipticEvaluator*)MYMalloc(sizeof(PNEllipticEvaluator));
    ret->is_zeroe = (e <= 1e-16);
    ret->beta = eval_beta(e);
    ret->e = e;
    ret->atol = atol;
    ret->rtol = rtol;
    ret->cacheL = CreateLaplaceCache2D(2, 2, ret->beta);
    ret->cacheBJ = CreateBesselJCache2D(2, 2, ret->e);
    ret->cacheDL = CreateDLaplaceCache2D(2, 2, ret->beta);

    ret->cacheJ0 = CreatePNEllipticCache(2, 2, 2, ret->e, 1.0);
    ret->cacheJ1 = CreatePNEllipticCache(2, 2, 2, ret->e, 0.0);
    ret->cacheJ2 = CreatePNEllipticCache(2, 2, 2, ret->e, 0.0);
    ret->cacheJ3 = CreatePNEllipticCache(2, 2, 2, ret->e, 0.0);
    ret->cacheK0 = CreatePNEllipticCache(2, 2, 2, ret->e, log((1. + sqrt(1. - e*e)) / 2.));
    return ret;
}

void DestroyPNEllipticEvaluator(PNEllipticEvaluator *cache)
{
    if (!cache) return;
    STRUCTFREE(cache->cacheBJ, BesselJCache2D);
    STRUCTFREE(cache->cacheL,  LaplaceCache2D);
    STRUCTFREE(cache->cacheDL, DLaplaceCache2D);

    STRUCTFREE(cache->cacheJ0, PNEllipticCache);
    STRUCTFREE(cache->cacheJ1, PNEllipticCache);
    STRUCTFREE(cache->cacheJ2, PNEllipticCache);
    STRUCTFREE(cache->cacheJ3, PNEllipticCache);
    STRUCTFREE(cache->cacheK0, PNEllipticCache);

    MYFree(cache);
    cache = NULL;
    return;
}

REAL8 evaluate_Jpqab(INT p, INT q, size_t a, size_t b, PNEllipticEvaluator* cache)
{
    if (cache->is_zeroe) {
        return (p==0 && b==0 ? 1.0 : 0.0);
    }
    switch (b)
    {
        case 0:
            return get_Jpqa0_from_cache(p, q, a, cache->cacheBJ, cache->cacheL, cache->cacheJ0, cache->atol, cache->rtol);
        case 1:
            return get_Jpqa1_from_cache(p, q, a, cache->cacheBJ, cache->cacheL, cache->cacheJ1, cache->atol, cache->rtol);
        case 2:
            return get_Jpqa2_from_cache(p, q, a, cache->cacheBJ, cache->cacheL, cache->cacheJ2, cache->atol, cache->rtol);
        case 3:
            return get_Jpqa3_from_cache(p, q, a, cache->cacheBJ, cache->cacheL, cache->cacheJ3, cache->atol, cache->rtol);
        default:
            X_PRINT_INFO("Unsupported J[%d,%d,%zu,%zu]", p, q, a, b);
            X_ERROR_REAL8(X_EINVAL);
    }
    return X_REAL8_FAIL_NAN;
}

REAL8 evaluate_Kpqa0(INT p, INT q, size_t a, PNEllipticEvaluator* cache)
{
    if (cache->is_zeroe) return 0.0;
    return get_Kpqa0_from_cache(p, q, a, cache->cacheBJ, cache->cacheL, cache->cacheDL, cache->cacheK0, cache->atol, cache->rtol);
}
