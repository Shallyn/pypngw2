/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW3__
#define __INCLUDE_PNGW3__
#include "pngw3_Utils.h"
#include "pngw3_evolve.h"
#include "src_hlm/pngw3_hlms.h"

#if defined(USE_OPENMP) && defined(_OPENMP)
#include <omp.h>
#else
  #define omp_get_max_threads() 1
  #define omp_get_thread_num()  0
#endif

REAL8 evaluate_chi_from_l(REAL8 l, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
REAL8 evaluate_l_from_chi(REAL8 chi, REAL8 eta, DynParamsCache *dc);
INT evaluate_PACorr_from_tdynamics(REAL8 eta, REAL8 chi0, 
    PNGW3Dynamics *dyn);
INT evaluate_hSphericalModes_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSphericalModes_TD **ret_hlms,
        REAL8 eps_h, REAL8 eps_pe);
INT evaluate_hSphericalModes_pce10_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSphericalModes_TD **ret_hlms);


INT evaluate_h22_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSingleSphericalMode_TD **ret_h22, REAL8 eps_h, REAL8 eps_pe);

INT evaluate_h22_pce10_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSingleSphericalMode_TD **ret_h22, REAL8 eps_h, REAL8 eps_pe);


#endif