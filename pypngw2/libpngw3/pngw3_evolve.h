/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW3_EVOLVE__
#define __INCLUDE_PNGW3_EVOLVE__
#include "pngw3_Utils.h"

void evaluate_DynParamsCache(REAL8 v, REAL8 e, DynParamsCache *cache);
INT evolve_binary_orbit_t(
    PNGW3Core       *pnparams,
    PNGW3Dynamics   **outdyn
);

#endif
