/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pngw3.h"
#include "../libpnfourier/pnFourier.h"
#include <gsl/gsl_sf.h>

void test_func();

void test_func1()
{
    REAL8 atol = 1e-16, rtol = 1e-16;
    REAL8 Jpqa0Ans, Jpqa1Ans, Jpqa2Ans, Jpqa3Ans, Kpqa0Ans;
    INT p, q, a;
    REAL8 e;
    INT is_case1 = 1;
    INT is_case2 = 1;
    INT is_case3 = 1;
    INT is_case4 = 1;
    INT is_case5 = 1;
    INT is_case6 = 1;
    PNEllipticEvaluator *pe = NULL;
    e = 0.9;
    pe = CreatePNEllipticEvaluator(e, atol, rtol);
    // Total Test
    // case 1:
    p = 3, q = -5, a = 7;
    Jpqa0Ans=310931.6876461047;
    Jpqa1Ans=-135182.2663518158;
    Jpqa2Ans=2348.280690146577;
    Jpqa3Ans=44459.83646701042;
    Kpqa0Ans=-718440.0056949718;    
    if (is_case1) {
        print_debug("case1:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }
    
    p = -19, q = -20, a = 6;
    Jpqa0Ans=56970.26096256902;
    Jpqa1Ans=5077.019756795464;
    Jpqa2Ans=-11565.90492753995;
    Jpqa3Ans=-3320.702208368717;
    Kpqa0Ans=-126137.0194186964;
    if (is_case2) {
        print_debug("case2:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }

    p = 0, q = 13, a = 6;
    Jpqa0Ans=13613.89580896292;
    Jpqa1Ans=10185.2646797434;
    Jpqa2Ans=4751.067716858367;
    Jpqa3Ans=-413.2782571845691;
    Kpqa0Ans=-33405.6448563515;   
    if (is_case3) {    
        print_debug("case3:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }   

    p = 9, q = 0, a = 7;
    Jpqa0Ans=248132.5790984056;
    Jpqa1Ans=-126605.68935994;
    Jpqa2Ans=20495.16106928682;
    Jpqa3Ans=31459.00920958112;
    Kpqa0Ans=-581296.8778434412;
    if (is_case4) {    
        print_debug("case4:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }

    p = 0, q = 9, a = 0;
    Jpqa0Ans=0.1475174540443709;
    Jpqa1Ans=0.2730778616346681;
    Jpqa2Ans=-0.1738545032469508;
    Jpqa3Ans=-0.3102584188361371;
    Kpqa0Ans=-0.01321082045424299;  
    if (is_case5) {      
        print_debug("case5:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }
    p = 10, q = -9, a = 0;
    Jpqa0Ans=0.0659429236049341;
    Jpqa1Ans=0.04005097295501511;
    Jpqa2Ans=0.01869454008377534;
    Jpqa3Ans=0.003547388634723625;
    Kpqa0Ans=0.04990993968416432;
    if (is_case6) {      
        print_debug("case6:\nRel dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans)/Jpqa0Ans,
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans)/Jpqa1Ans,
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans)/Jpqa2Ans,
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans)/Jpqa3Ans,
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)/Kpqa0Ans
        );
        print_err("Abs dJpqa0 = %.16e, dJpqa1 = %.16e, dJpqa2 = %.16e, dJpqa3 = %.16e, dKpqa0 = %.16e\n\n",
            (evaluate_Jpqab(p, q, a, 0, pe) - Jpqa0Ans),
            (evaluate_Jpqab(p, q, a, 1, pe) - Jpqa1Ans),
            (evaluate_Jpqab(p, q, a, 2, pe) - Jpqa2Ans),
            (evaluate_Jpqab(p, q, a, 3, pe) - Jpqa3Ans),
            (evaluate_Kpqa0(p, q, a, pe) - Kpqa0Ans)
        );
    }
    STRUCTFREE(pe, PNEllipticEvaluator);
}

void test_hlm_waveform_mode(int p, REAL8 e)
{
    DynParamsCache dc;
    REAL8 eta = 0.2;
    REAL8 v = 1./5.;
    evaluate_DynParamsCache(v, e, &dc);
    PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);
    COMPLEX16 hlm = calculate_h52FourierModep(p, eta, &dc, pe);
    print_debug("hlm_(%d) = (%.16e) + j(%.16e)\n", p, creal(hlm), cimag(hlm));
    STRUCTFREE(pe, PNEllipticEvaluator);
}

void test_PAOscCorrection(REAL8 eta, REAL8 chi, REAL8 e, REAL8 v)
{
    DynParamsCache dc;
    evaluate_DynParamsCache(v, e, &dc);
    evaluate_l_from_chi(chi, eta, &dc);
    // PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);
    // evaluate_KeplerEQCoeff(12, 0.2, &dc, pe);
    // STRUCTFREE(pe, PNEllipticEvaluator);
    // REAL8 eOsc, vOsc, nOsc, clOsc, clbdOsc;
    // evaluate_PACorrections(eta, chi, &dc, &eOsc, &vOsc, &nOsc, &clOsc, &clbdOsc);
    // print_debug("eOsc = %.16e, vOsc = %.16e, nOsc = %.16e\n", eOsc, vOsc, nOsc);
    // print_debug("clOsc = %.16e, clbdOsc = %.16e\n", clOsc, clbdOsc);
    return;
}

void test_evolvet(REAL8 eta, REAL8 e0, REAL8 v0, REAL8 zeta0, 
                  REAL8 phic, REAL8 tc, REAL8 iota)
{
    PNGW3Dynamics *dyn;
    PNGW3HyperParams hpms;
    hpms.EPS_ABS = 1e-16;
    hpms.EPS_REL = 1e-16;
    hpms.PNOrder = 6;
    hpms.vmax = VISCO;
    PNGW3Core *pnparams = CreatePNGW3Core(eta, e0, v0, zeta0, phic, tc, iota, &hpms);
    INT ret = evolve_binary_orbit_t(pnparams, &dyn);
    print_debug("is success ? %d\n", ret);
    print_debug("len = %d\n", dyn->length);

    STRUCTFREE(dyn, PNGW3Dynamics);
    return;
}


void test_hlm_waveform_mode_sum(REAL8 eta, REAL8 l, REAL8 e, REAL8 v)
{
    DynParamsCache dc;
    evaluate_DynParamsCache(v, e, &dc);
    PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);
    REAL8 hr, hi;
    evaluate_h22(eta, l, &dc, pe, &hr, &hi, 1e-10);
    print_debug("hlm = (%.16e) + j(%.16e)\n", hr, hi);
    STRUCTFREE(pe, PNEllipticEvaluator);
}


int test_hlm_td(REAL8 eta, REAL8 e0, REAL8 chi0)
{
    PNGW3HyperParams hpms;
    hpms.EPS_ABS = 1e-16;
    hpms.EPS_REL = 1e-16;
    REAL8 v0 = 0.16683626292147538;
    hpms.vmax = VISCO;
    PNGW3Core *pnparams = CreatePNGW3Core(eta, e0, v0, chi0, 0.0, 0.0, 0.0, &hpms);
    PNGW3Dynamics *dyn = NULL;
    if (evolve_binary_orbit_t(pnparams, &dyn) != X_SUCCESS)
        return X_FAILURE;
    if (evaluate_PACorr_from_tdynamics(pnparams->eta, pnparams->initial_chi, dyn) != X_SUCCESS)
        return X_FAILURE;

    PNGW3HSphericalModes_TD *hlms = NULL;
    if (evaluate_hSphericalModes_from_PADynamics(pnparams->eta, dyn, &hlms,
        1e-5, 1e-8) != X_SUCCESS)
        return X_FAILURE;

    STRUCTFREE(dyn, PNGW3Dynamics);
    STRUCTFREE(hlms, PNGW3HSphericalModes_TD);
    STRUCTFREE(pnparams, PNGW3Core);
    return X_SUCCESS;
}

void test_hlm_zeroe(REAL8 eta, REAL8 v) 
{
    DynParamsCache dc;
    evaluate_DynParamsCache(v, 0.0, &dc);
    COMPLEX16 hlm;
    // L=2
    hlm = calculate_h22FourierMode_zeroe(eta, &dc);
    print_debug("h22 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h21FourierMode_zeroe(eta, &dc);
    print_debug("h21 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    // L=3
    hlm = calculate_h33FourierMode_zeroe(eta, &dc);
    print_debug("h33 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h32FourierMode_zeroe(eta, &dc);
    print_debug("h32 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h31FourierMode_zeroe(eta, &dc);
    print_debug("h31 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    // L=4
    hlm = calculate_h44FourierMode_zeroe(eta, &dc);
    print_debug("h44 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h43FourierMode_zeroe(eta, &dc);
    print_debug("h43 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h42FourierMode_zeroe(eta, &dc);
    print_debug("h42 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h41FourierMode_zeroe(eta, &dc);
    print_debug("h41 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    // L=5
    hlm = calculate_h55FourierMode_zeroe(eta, &dc);
    print_debug("h55 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h54FourierMode_zeroe(eta, &dc);
    print_debug("h54 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h53FourierMode_zeroe(eta, &dc);
    print_debug("h53 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h52FourierMode_zeroe(eta, &dc);
    print_debug("h52 = %.16e , %.16e\n", creal(hlm), cimag(hlm));
    hlm = calculate_h51FourierMode_zeroe(eta, &dc);
    print_debug("h51 = %.16e , %.16e\n", creal(hlm), cimag(hlm));

    return;
}

void test_hlm_cumsum(REAL8 eta, REAL8 v, REAL8 e, REAL8 l)
{
    DynParamsCache dc;
    evaluate_DynParamsCache(v, e, &dc);
    PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);
    int pmax = 10;
    REAL8 hrPos=0.0, hrNeg=0.0, hiPos=0.0, hiNeg=0.0;
    REAL8 cum, eps = 1e-8;
    size_t ici = 0, cons_iter = 3;
    for (int pp = 1; ; pp++) {
        COMPLEX16 coeff = calculate_h54FourierModep(pp+(4), eta, &dc, pe);
        COMPLEX16 term = coeff * cexp(l*I*pp);
        REAL8 abscoeff = cabs(coeff);
        hrPos += creal(term);
        hiPos += cimag(term);
        cum += abscoeff;
        print_debug("[%d] hPos = %.16e, %.16e\n", pp, creal(coeff) , cimag(coeff) );
        if (abscoeff < eps*cum) {
            if (ici >cons_iter)
                break;
            else
                ici++;
        } else 
            ici = 0;
    }
    cum = 0.0;
    for (int pp=1; ; pp++) {
        COMPLEX16 coeff = calculate_h54FourierModep(-pp+(4), eta, &dc, pe);
        COMPLEX16 term = coeff * cexp(-l*I*pp);
        REAL8 abscoeff = cabs(coeff);
        hrNeg += creal(term);
        hiNeg += cimag(term);
        cum += abscoeff;
        print_debug("[-%d] hNeg = %.16e, %.16e\n", pp, creal(coeff) , cimag(coeff) );
        if (abscoeff < eps*cum) {
            if (ici >cons_iter)
                break;
            else
                ici++;
        } else 
            ici = 0;
    }
    print_debug("total h = %.16e, %.16e\n", hrPos + hrNeg, hiPos + hiNeg);

    STRUCTFREE(pe, PNEllipticEvaluator);
}

void test_h22_waveform_mode(int p, REAL8 eta, REAL8 v, REAL8 e)
{
    DynParamsCache dc;
    evaluate_DynParamsCache(v, e, &dc);
    PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);
    COMPLEX16 hlm = calculate_h22FourierModep(p, eta, &dc, pe);
    print_debug("hlm_(%d) = (%.16e) + j(%.16e)\n", p, creal(hlm), cimag(hlm));
    STRUCTFREE(pe, PNEllipticEvaluator);
}

int main()
{
    test_h22_waveform_mode(5, 0.2, 0.25, 0.5);
    // test_hlm_cumsum(0.222222222222222, 0.21222275405251909, 0.4917256553537161, CST_PI/3.);
    // REAL8 eta = 0.2222222222222222;
    // REAL8 l = CST_PI/3.;
    // REAL8 e = 0.2;
    // REAL8 v = 0.16683626292147538;
    // // test_hlm_waveform_mode_sum(eta, l, e, v);
    // test_hlm_td(eta, e, 0);

// test_hlm_zeroe(0.2, 0.15);
    // test_evolvet(0.2, 0.4, 0.2, 0, 0, 0, 0);
    // test_PAOscCorrection(0.2, CST_PI/3., 0.77, 0.2);
    // print_debug("p = -12\n");
    // test_hlm_waveform_mode(-12, 0.9);
    // print_err("\n");
    // print_debug("p = 12\n");
    // test_hlm_waveform_mode(12, 0.9);
    // test_func1();
    // DynParamsCache dc;
    // REAL8 eta = 0.25;
    // REAL8 v = 1./5., e = 0.3;
    // evaluate_DynParamsCache(v, e, &dc);
    // PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-16, 1e-16);

    // print_debug("J0101 = %.16e\n", evaluate_Jpqab(0, 1, 0, 1, pe));
    // print_debug("J = %.16e\n", J_pqa1_series(p-1, p-2, 0, e, 1e-16, 1e-16));
    // COMPLEX16 h22 = calculate_h22InstFourierModep(p, eta, &dc, pe);
    // print_debug("h22 = (%.16e) + j(%.16e)\n", creal(h22), cimag(h22));
    // REAL8 xArg = -1.5;
    // BesselJCache *bc = CreateBesselJCache(2, xArg);
    // int k = 3;
    // // for (int k=-20; k<=20; k++){
    //     REAL8 val1 = get_BesselJ_from_BesselJCache(abs(k), bc);
    //     REAL8 val2 = jn(abs(k), xArg);
    //     print_debug("k=%d, (%.16e, %.16e), ddiff = %.16e\n", k, val1, val2, val1-val2);
    // // }
    // DestroyBesselJCache(bc);
    // REAL8 e = 0.3;
    // REAL8 beta = (1.0 - sqrt(1.0 - e*e)) / e;
    // LaplaceCache2D* lc = CreateLaplaceCache2D(2, 2, beta);
    // DLaplaceCache2D* dlc = CreateDLaplaceCache2D(2, 2, beta);
    // BesselJCache2D* jc = CreateBesselJCache2D(2, 2, e);
    // PNEllipticCache *cache = CreatePNEllipticCache(2, 2, 2, e, 0.0);
    // REAL8 atol = 1e-16, rtol = 1e-16;
    // // REAL8 val1 = J_pqa3_series(-20, -20, 20, e, atol, rtol);
    // // REAL8 val2 = -J_pqa3_series(20, 20, 20, e, atol, rtol);
    // // // REAL8 val2 = J_pqa3_series_cache(-20, -20, 20, cache->e, jc, lc, atol, rtol);
    // // print_debug("%.16e, %.16e, %.16e\n",
    // //     val1, val2, val1-val2
    // // );
    // for (INT ip = -20; ip<=20; ip++)
    //     for (INT iq = -20; iq<=20; iq++) 
    //         for (size_t ia=0; ia<=20; ia++) {
    //         // print_debug("(p, q, a) = (%d, %d, %d)>>\n", ip, iq, ia);
    //         REAL8 val1 = K_pqa0_series(ip, iq, ia, e, atol, rtol);
    //         // REAL8 val2 = get_Jpqa1_from_PNEllipticCache(ip, iq, ia, cache, atol, rtol);
    //         // REAL8 val2 = get_Jpqa2_from_PNEllipticCache(ip, iq, ia, cache, atol, rtol);
    //         REAL8 val2 = get_Kpqa0_from_cache(ip, iq, ia, jc, lc, dlc, cache, atol, rtol);
    //         // REAL8 val2 = 0.0;
    //         if (fabs(val1 - val2)>1e-15*fabs(val1) && fabs(val1 - val2)>1e-15)
    //         print_err("(p, q, a) = (%d, %d, %d)>>(%.16e, %.16e): diff = %.16e\n", 
    //                 ip, iq, ia,
    //                 val1, val2, 
    //                 val1 - val2);
    //     }
    // DestroyLaplaceCache2D(lc);
    // DestroyBesselJCache2D(jc);
    // DestroyPNEllipticCache(cache);

    // REAL8 e = 0.3;
    // REAL8 beta = (1.0 - sqrt(1.0 - e*e)) / e;
    // LaplaceCache2D* lc = CreateLaplaceCache2D(2, 2, beta);
    // BesselJCache2D* jc = CreateBesselJCache2D(2, 2, e);
    // REAL8 atol = 1e-16, rtol = 1e-16;
    // for (INT ip = -10; ip<=10; ip++)
    //     for (INT iq = -10; iq<=10; iq++) 
    //         for (size_t ia=0; ia<=10; ia++){
    //         REAL8 val1 = J_pqa0_series(ip, iq, ia, e, atol, rtol);
    //         REAL8 val2 = J_pqa0_series_cache(ip, iq, ia, e, jc, lc, atol, rtol);
    //         print_debug("(p, q, a) = (%d, %d, %d) = (%.5e, %.5e): diff = %.16e\n", 
    //                 ip, iq, ia,
    //                 val1, val2, val1 - val2);
    //     }
    // DestroyLaplaceCache2D(lc);
    // DestroyBesselJCache2D(jc);

    // REAL8 beta = 0.3;
    // LaplaceCache2D* cache = CreateLaplaceCache2D(2, 2, beta);
    // for (size_t in = 0; in<=20; in++)
    //     for (size_t ia=0; ia<=20; ia++) {
    //         REAL8 val1 = laplace_na(in, ia, beta);
    //         REAL8 val2 = get_Laplace_from_LaplaceCache2D(in, ia, cache);
    //         print_debug("(n, a) = (%d, %d) = (%.5e, %.5e): diff = %.16e\n", in, ia,
    //                 val1, val2, val1 - val2);
    //     }
    // DestroyLaplaceCache2D(cache);

    // REAL8 beta = 0.3;
    // DLaplaceCache2D* cache = CreateDLaplaceCache2D(10, 10, beta);
    // for (size_t in = 0; in<=10; in++)
    //     for (size_t ia=0; ia<=10; ia++) {
    //         REAL8 val1 = Dlaplace_na(in, ia, beta);
    //         REAL8 val2 = get_DLaplace_from_DLaplaceCache2D(in, ia, cache);
    //         print_debug("(n, a) = (%d, %d) = (%.5e, %.5e): diff = %.16e\n", in, ia,
    //                 val1, val2, val1 - val2);
    //     }
    // DestroyDLaplaceCache2D(cache);

    // REAL8 e = 0.3;
    // BesselJCache2D* cache = CreateBesselJCache2D(2, 2, e);
    // for (INT ik = -20; ik<=20; ik++)
    //     for (INT iq=-20; iq<=20; iq++) {
    //         REAL8 val1 = jn(ik, iq*e);
    //         REAL8 val2 = get_BesselJ_from_BesselJCache2D(ik, iq, cache);
    //         print_debug("(k, q) = (%d, %d) = (%.5e, %.5e): diff = %.16e\n", ik, iq,
    //                 val1, val2, val1 - val2);
    //     }
    // DestroyBesselJCache2D(cache);

    // test_func();
    // REAL8Array *arr = CreateREAL8ArrayL(2,2,10);
    // print_debug("n = %.16e\n", test_func(v, e));
    // STRUCTFREE(arr, REAL8Array);
    CheckMemoryLeaks();
    return 1;
}
