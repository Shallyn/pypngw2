/**
* Writer: Xiaolin.liu
* shallyn.liu@foxmail.com
**/

#ifndef __INCLUDE_PNFOURIER__
#define __INCLUDE_PNFOURIER__

#include "pnFourier_Jpqa0.h"
#include "pnFourier_Jpqab.h"
#include "pnFourier_Kpqab.h"

typedef struct {
    REAL8   e;
    REAL8   beta;
    size_t  abspMax;
    size_t  absqMax;
    size_t  aMax;
    size_t  size;   /* len = (2*absqMax + 1)*(abspMax + 1)*(aMax + 1) */
    size_t  stride_a; /* (abspMax + 1)*(2*absqMax + 1) */
    size_t  stride_p; /* 2*absqMax + 1 */
    REAL8   *buf;   /* J[ip, iq, ia] = iq + absqMax + ip*stride_p + ia*stride_a */
}PNEllipticCache;


typedef struct {
    BOOLEAN is_zeroe;
    REAL8 e;
    REAL8 beta;
    REAL8 atol;
    REAL8 rtol;
    LaplaceCache2D *cacheL;
    DLaplaceCache2D *cacheDL;
    BesselJCache2D *cacheBJ;
    PNEllipticCache *cacheJ0;
    PNEllipticCache *cacheJ1;
    PNEllipticCache *cacheJ2;
    PNEllipticCache *cacheJ3;
    PNEllipticCache *cacheK0;
}PNEllipticEvaluator;

PNEllipticCache *CreatePNEllipticCache(size_t abspMax, size_t absqMax, size_t aMax, REAL8 e, REAL8 val0);
void DestroyPNEllipticCache(PNEllipticCache *cache);

REAL8 get_Jpqa0_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol);
REAL8 get_Jpqa0_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol);

REAL8 get_Jpqa1_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol);
REAL8 get_Jpqa1_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol);


REAL8 get_Jpqa2_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol);
REAL8 get_Jpqa2_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol);

REAL8 get_Jpqa3_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol);
REAL8 get_Jpqa3_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, 
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol);

REAL8 get_Kpqa0_from_PNEllipticCache(INT p, INT q, INT a, PNEllipticCache *cache, REAL8 atol, REAL8 rtol);
REAL8 get_Kpqa0_from_cache(INT p, INT q, INT a, 
    BesselJCache2D *jc, LaplaceCache2D *lc, DLaplaceCache2D *dlc,
    PNEllipticCache *cache, REAL8 atol, REAL8 rtol);

REAL8 evaluate_Jpqab(INT p, INT q, size_t a, size_t b, PNEllipticEvaluator* cache);
REAL8 evaluate_Kpqa0(INT p, INT q, size_t a, PNEllipticEvaluator* cache);

PNEllipticEvaluator* CreatePNEllipticEvaluator(REAL8 e, REAL8 atol, REAL8 rtol);
void DestroyPNEllipticEvaluator(PNEllipticEvaluator *cache);

#endif

