/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW_HLMS__
#define __INCLUDE_PNGW_HLMS__
#include "../pngw3_Utils.h"

COMPLEX16 calculate_h22FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h21FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h20FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);

COMPLEX16 calculate_h22FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h21FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h20FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);

COMPLEX16 calculate_h22_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h21_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h20_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);

INT evaluate_h22(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h21(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h20(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);

COMPLEX16 calculate_h33FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h32FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h31FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h30FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);

COMPLEX16 calculate_h33FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h32FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h31FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h30FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);

COMPLEX16 calculate_h33_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h32_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h31_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h30_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);

INT evaluate_h33(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h32(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h31(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h30(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);

COMPLEX16 calculate_h44FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h43FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h42FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h41FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h40FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);

COMPLEX16 calculate_h44FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h43FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h42FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h41FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h40FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);

COMPLEX16 calculate_h44_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h43_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h42_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h41_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h40_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);

INT evaluate_h44(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h43(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h42(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h41(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h40(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);

COMPLEX16 calculate_h55FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h54FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h53FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h52FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h51FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);
COMPLEX16 calculate_h50FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe);

COMPLEX16 calculate_h55FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h54FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h53FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h52FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h51FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);
COMPLEX16 calculate_h50FourierMode_zeroe(REAL8 eta, DynParamsCache *dc);


COMPLEX16 calculate_h55_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h54_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h53_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h52_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h51_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);
COMPLEX16 calculate_h50_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc);

INT evaluate_h55(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h54(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h53(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h52(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h51(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);
INT evaluate_h50(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,
    REAL8 *hr, REAL8 *hi, REAL8 eps);

#endif