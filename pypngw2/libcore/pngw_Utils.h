/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW_UTILS__
#define __INCLUDE_PNGW_UTILS__
#include "pngw_Datatypes.h"

#define INTERNAL_DEBUG_LEVEL get_debug_level()
#define SET_INTERNAL_DEBUG_LEVEL(flag) set_debug_level(flag)

void set_debug_level(INT flag);
INT get_debug_level();
INT get_ilk_from_nm(INT n, INT m);
REAL8 calculate_k(BBHDynVariables *var, BBHDynParams *pms);
REAL8 calculate_k_from_core(BBHCore *core);
BBHDynParams *CreateBBHDynParams(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 kappa1, REAL8 kappa2, 
    REAL8 Theta, REAL8 Phi, REAL8 eini, REAL8 vomini, INT PN_Ord2);
BBHCore *CreateBBHCore(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 kappa1, REAL8 kappa2,
    REAL8 Theta, REAL8 Phi, REAL8 e0, REAL8 vom, INT PN_Ord2);
void SetBBHDynVariables(REAL8 e0, REAL8 vom, BBHDynVariables *var);
void SetBBHCoreDynVariables(REAL8 e0, REAL8 vom, BBHCore *core);

AntennaPatternF *CreateAntennaPatternF(REAL8 psi, REAL8 ra, REAL8 dec,
    REAL8 lbd, REAL8 varphi, REAL8 gamma, REAL8 zeta);
void DestroyAntennaPatternF(AntennaPatternF *apf);

REAL8 calculate_barFplus_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFcross_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFplus_deriv_psi_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFplus_deriv_alpha_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFplus_deriv_delta_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFcross_deriv_psi_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFcross_deriv_alpha_t(REAL8 t_SI, AntennaPatternF *apf);
REAL8 calculate_barFcross_deriv_delta_t(REAL8 t_SI, AntennaPatternF *apf);


#endif
