/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW3_EHFUNC__
#define __INCLUDE_PNGW3_EHFUNC__
#include "pngw3_Utils.h"

REAL8 evaluate_ehfunc_varphi(REAL8 e2);
REAL8 evaluate_ehfunc_beta(REAL8 e2);
REAL8 evaluate_ehfunc_gamma(REAL8 e2);
REAL8 evaluate_ehfunc_varphi01(REAL8 eta, REAL8 e2);
REAL8 evaluate_ehfunc_varphi10(REAL8 e2);
REAL8 evaluate_ehfunc_F(REAL8 e2);
REAL8 evaluate_ehfunc_chi(REAL8 e2);

REAL8 evaluate_ehfunc_varphiTilde(REAL8 e2);
REAL8 evaluate_ehfunc_betaTilde(REAL8 e2);
REAL8 evaluate_ehfunc_gammaTilde(REAL8 e2);
REAL8 evaluate_ehfunc_varphiTilde01(REAL8 eta, REAL8 e2);
REAL8 evaluate_ehfunc_varphiTilde10(REAL8 e2);
REAL8 evaluate_ehfunc_FTilde(REAL8 e2);
REAL8 evaluate_ehfunc_chiTilde(REAL8 e2);

REAL8 evaluate_ehfunc_varphiEPN3(REAL8 e2);
REAL8 evaluate_ehfunc_betaEPN5(REAL8 e2);
REAL8 evaluate_ehfunc_gammaEPN5(REAL8 e2);
REAL8 evaluate_ehfunc_varphiE01PN5(REAL8 eta, REAL8 e2);
REAL8 evaluate_ehfunc_varphiE10PN5(REAL8 eta, REAL8 e2);
REAL8 evaluate_ehfunc_FEPN6(REAL8 e2);
REAL8 evaluate_ehfunc_chiEPN6(REAL8 e2);

#endif