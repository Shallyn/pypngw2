/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
**/

#include "pngw3_hlms.h"

// Here we need to add p=0 mode, which is not included in the Fourier expansion. This mode can be calculated by directly evaluating the integrals, which is much more efficient than calculating the Fourier modes and summing them up.
#define DEFFUNC(L, M) \
INT evaluate_h##L##M(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,\
    REAL8 *hr, REAL8 *hi, REAL8 eps)\
{\
    if(pe->is_zeroe) {\
        COMPLEX16 h = cexp(I*l*(-M))*calculate_h##L##M##FourierMode_zeroe(eta, dc);\
        *hr = creal(h);\
        *hi = cimag(h);\
        return X_SUCCESS;\
    }\
    REAL8 hrPos=0.0, hrNeg=0.0, hiPos=0.0, hiNeg=0.0;\
    size_t min_iter = 5;\
    size_t cons_iter = 3, ici = 0;\
    REAL8 cum;\
    cum = 0.0;\
    for (int pp = 1;  ; pp++)\
    {\
        COMPLEX16 coeff = calculate_h##L##M##FourierModep(pp+(M), eta, dc, pe) ;\
        COMPLEX16 term = coeff * cexp(l*I*pp);\
        REAL8 abscoeff = cabs(coeff);\
        hrPos += creal(term);\
        hiPos += cimag(term);\
        cum += abscoeff;\
        if (pp > min_iter && abscoeff < eps*cum) {\
            if (ici >cons_iter)\
                break;\
            else\
                ici++;\
        } else \
            ici = 0;\
    }\
    cum = 0.0;\
    for (int pp = 1;  ; pp++)\
    {\
        COMPLEX16 coeff = calculate_h##L##M##FourierModep(-pp+(M), eta, dc, pe) ;\
        COMPLEX16 term = coeff * cexp(-l*I*pp);\
        REAL8 abscoeff = cabs(coeff);\
        hrNeg += creal(term);\
        hiNeg += cimag(term);\
        cum += abscoeff;\
        if (pp > min_iter && abscoeff < eps*cum) {\
            if (ici >cons_iter)\
                break;\
            else\
                ici++;\
        } else \
            ici = 0;\
    }\
    *hr = hrPos + hrNeg;\
    *hi = hiPos + hiNeg;\
    return X_SUCCESS;\
}

// #define DEFFUNC(L, M) \
// INT evaluate_h##L##M(REAL8 eta, REAL8 l, DynParamsCache *dc, PNEllipticEvaluator *pe,\
//     REAL8 *hr, REAL8 *hi, REAL8 eps)\
// {\
//     if(pe->is_zeroe) {\
//         COMPLEX16 h = cexp(I*l*(-M))*calculate_h##L##M##FourierMode_zeroe(eta, dc);\
//         *hr = creal(h);\
//         *hi = cimag(h);\
//         return X_SUCCESS;\
//     }\
//     REAL8 hrPos=0.0, hrNeg=0.0, hiPos=0.0, hiNeg=0.0;\
//     size_t min_iter = 5;\
//     size_t cons_iter = 3, ici = 0;\
//     for (int pp = 1; pp<10  ; pp++)\
//     {\
//         COMPLEX16 coeff = calculate_h##L##M##FourierModep(pp+(M), eta, dc, pe) ;\
//         COMPLEX16 term = coeff * cexp(l*I*pp);\
//         REAL8 abscoeff = cabs(coeff);\
//         hrPos += creal(term);\
//         hiPos += cimag(term);\
//     }\
//     for (int pp = 1; pp<10  ; pp++)\
//     {\
//         COMPLEX16 coeff = calculate_h##L##M##FourierModep(-pp+(M), eta, dc, pe) ;\
//         COMPLEX16 term = coeff * cexp(-l*I*pp);\
//         REAL8 abscoeff = cabs(coeff);\
//         hrNeg += creal(term);\
//         hiNeg += cimag(term);\
//     }\
//     *hr = hrPos + hrNeg;\
//     *hi = hiPos + hiNeg;\
//     return X_SUCCESS;\
// }


// h22
COMPLEX16 calculate_h22FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe)
{
    COMPLEX16 ret;
    INT pm = p - 2;
    if (pm==0) return 0.0; // FIXME: this is not zero
    INT pm2 = pm*pm;
    INT pm3 = pm2*pm;
    INT pm4 = pm2*pm2;
    INT pm5 = pm4*pm;
    REAL8 *ePow = dc->ePow;
    REAL8 *sq1m2Pow = dc->sq1mesqPow;
    REAL8 m1p3eta = -1. + 3.*eta;
    REAL8 y = sq1m2Pow[1];
    REAL8 yp = y+1.0, myp = -y+1.0;
    REAL8 yp2 = yp*yp, myp2 = myp*myp;
    REAL8 yp3 = yp2*yp, myp3 = myp2*myp;
    REAL8 yp4 = yp3*yp, myp4 = myp3*myp;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
#define EVALEJ(pp, qq, aa, bb) evaluate_Jpqab((pp), (qq), (aa), (bb), pe)
#define EVALEK(pp, qq, aa, bb) evaluate_Kpqa0((pp), (qq), (aa), pe)

// #define EVALEJ(p, q, a, b) 1.0
#define SQ1M2(n) sq1m2Pow[n]
#define EPOW(n) ePow[n]
    REAL8 PN0, PN2, PN4, ImPN5, PN6;
    PN0 = 2*EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*SQ1M2(1) - 2*EVALEJ(pm,pm,1,0)*SQ1M2(2) + EVALEJ(pm,pm,0,0)*(-1 - 2*SQ1M2(1) - SQ1M2(2) - 2*pm*SQ1M2(3)) + 2*EVALEJ(pm,pm,2,0)*SQ1M2(4);
    PN0 /= ePow[2];

    PN2 = -3*myp2*pm*EPOW(1)*EVALEJ(-1 + pm,pm,0,1) - 
        3*pm*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,1) + 
        3*pm*yp2*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*SQ1M2(1) - 
        3*pm*yp2*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*SQ1M2(1) - 
        12*pm*EVALEJ(pm,pm,1,1)*SQ1M2(2) - 
        12*pm2*EVALEJ(pm,pm,0,1)*SQ1M2(3) - 
        (pm*EVALEJ(pm,pm,1,0)*SQ1M2(2)*
        (-90 + 81*eta + 252*SQ1M2(1) - (23 + 8*eta)*SQ1M2(2) + 
        2*(62 + 3*eta)*pm*SQ1M2(3)))/21. + 
        (EVALEJ(pm,pm,0,0)*(-252*yp2 + 
        2*pm2*(343 - 84*eta - 252*SQ1M2(1) + 11*(2 + eta)*SQ1M2(2))*
        SQ1M2(3) + pm*(1 + SQ1M2(1))*
        (-168 - 84*(-8 + eta)*SQ1M2(1) + (279 - 67*eta)*SQ1M2(2) + 
        (-57 + 17*eta)*SQ1M2(3))))/42. + 
        (EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*
        (252*SQ1M2(1) + pm*(126 + 42*(-6 + eta)*SQ1M2(1) + 126*SQ1M2(2) + 
        (-111 + 25*eta)*SQ1M2(3))))/21. + 
        12*pm*EVALEJ(pm,pm,2,1)*SQ1M2(4) + 
        (pm*EVALEJ(pm,pm,2,0)*(9*(-94 + 9*eta) + 252*SQ1M2(1) - 
        (19 + 20*eta)*SQ1M2(2) + 2*(-1 + 3*eta)*pm*SQ1M2(3))*SQ1M2(4))/
        21. + (2*(121 + 15*eta)*pm*EVALEJ(pm,pm,3,0)*SQ1M2(6))/21. - 
        (2*(-1 + 3*eta)*pm*EVALEJ(pm,pm,4,0)*SQ1M2(8))/7.;    
    PN2 /= ePow[2]*pm*sq1m2Pow[2];

    PN4 = -9*myp2*pm2*EPOW(1)*EVALEJ(-1 + pm,pm,0,2) - 
        9*pm2*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,2) + 
        18*pm2*yp2*EPOW(1)*EVALEJ(-1 + pm,pm,1,1)*SQ1M2(1) - 
        18*pm2*yp2*EPOW(1)*EVALEJ(1 + pm,pm,1,1)*SQ1M2(1) - 
        9*myp2*pm2*EVALEJ(-2 + pm,pm,0,1)*(-1 + SQ1M2(2)) + 
        9*pm2*yp2*EVALEJ(2 + pm,pm,0,1)*(-1 + SQ1M2(2)) - 
        36*pm2*EVALEJ(pm,pm,1,2)*SQ1M2(2) - 
        36*pm3*EVALEJ(pm,pm,0,2)*SQ1M2(3) + 
        (pm2*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,1)*
        (189*(-3 + 2*eta) - 168*(-4 + eta)*SQ1M2(1) + 
        (243 - 148*eta)*SQ1M2(2) + 21*(-5 + 2*eta)*pm*SQ1M2(3)))/28. + 
        (myp2*pm2*EPOW(1)*EVALEJ(-1 + pm,pm,0,1)*
        (189*(-3 + 2*eta) + 168*(-4 + eta)*SQ1M2(1) + 
        (243 - 148*eta)*SQ1M2(2) + 21*(-5 + 2*eta)*pm*SQ1M2(3)))/28. + 
        (pm2*EVALEJ(pm,pm,1,1)*SQ1M2(2)*
        (-555 + 216*eta - 504*SQ1M2(1) + (403 - 166*eta)*SQ1M2(2) + 
        (-353 + 30*eta)*pm*SQ1M2(3)))/7. + 
        36*pm2*EVALEJ(pm,pm,2,2)*SQ1M2(4) - 
        (pm2*EVALEJ(pm,pm,2,1)*(957 + 216*eta - 504*SQ1M2(1) + 
        (395 - 142*eta)*SQ1M2(2) + (-101 + 30*eta)*pm*SQ1M2(3))*SQ1M2(4))
        /7. - (pm2*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*SQ1M2(1)*
        (756*(-3 + 2*eta) + 24*(-185 + 72*eta)*SQ1M2(1) - 
        8*(-6 + 53*eta)*SQ1M2(2) + 
        (3224 + 7*eta2*pm - eta*(1328 + 105*pm))*SQ1M2(3) + 
        (972 - 592*eta)*SQ1M2(4) + 7*(-15 + eta)*eta*pm*SQ1M2(5)))/112.\
        + (pm2*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*SQ1M2(1)*
        (756*(-3 + 2*eta) + 24*(-185 + 72*eta)*SQ1M2(1) - 
        8*(-6 + 53*eta)*SQ1M2(2) + 
        (3224 + 7*eta2*pm - eta*(1328 + 105*pm))*SQ1M2(3) + 
        (972 - 592*eta)*SQ1M2(4) + 7*(-15 + eta)*eta*pm*SQ1M2(5)))/112.\
        + (pm2*EVALEJ(pm,pm,2,0)*SQ1M2(4)*
        (-27*(7692 - 4043*eta + 21*eta2) - 864*(431 + 9*eta)*SQ1M2(1) + 
        16*(14419 - 2570*eta + 386*eta2)*SQ1M2(2) - 
        4*(22896 - 7452*eta + (13435 + 10*eta + 128*eta2)*pm)*SQ1M2(3) + 
        3*(574 + 1271*eta - 133*eta2 + 288*(-1 + 3*eta)*pm)*SQ1M2(4) + 
        6*(-137 + 452*eta + 44*eta2)*pm*SQ1M2(5)))/1512. + 
        (pm2*EVALEJ(pm,pm,0,1)*SQ1M2(1)*
        (252 + 504*SQ1M2(1) + 7*(36 + (-7 + 30*eta)*pm)*SQ1M2(2) - 
        504*pm*SQ1M2(3) + (401 - 160*eta)*pm*SQ1M2(4) + 
        21*(-5 + 2*eta)*pm2*SQ1M2(5)))/7. + 
        (EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*
        (12096*SQ1M2(1) + 24*pm*
        (378 + (231 - 294*eta)*SQ1M2(1) + 378*SQ1M2(2) + 
        (-579 + 232*eta)*SQ1M2(3)) + 
        21*(-15 + eta)*eta*pm3*(1 + SQ1M2(2))*SQ1M2(4) + 
        2*pm2*(1134*(3 - 2*eta) - 42*(101 - 50*eta + 2*eta2)*SQ1M2(1) + 
        12*(-1119 + 95*eta)*SQ1M2(2) + 
        (3510 - 2252*eta + 296*eta2)*SQ1M2(3) + 
        6*(-313 + 176*eta)*SQ1M2(4) + 
        5*(-17 + 9*eta + 7*eta2)*SQ1M2(5))))/168. + 
        (4*(121 + 15*eta)*pm2*EVALEJ(pm,pm,3,1)*SQ1M2(6))/7. - 
        (pm2*EVALEJ(pm,pm,3,0)*(138733 + 59068*eta + 5132*eta2 - 
        148176*SQ1M2(1) + 3*(13839 - 9904*eta + 24*eta2)*SQ1M2(2) + 
        2*(-2707 + 1412*eta + 1060*eta2)*pm*SQ1M2(3))*SQ1M2(6))/1512. + 
        (pm2*EVALEJ(pm,pm,1,0)*SQ1M2(2)*
        (216*(279 + 227*eta) + 1512*(-13 + 18*eta)*SQ1M2(1) - 
        6*(-14459 + 6829*eta + 725*eta2)*SQ1M2(2) + 
        8*(21195 - 5184*eta + (53 + 3272*eta + 64*eta2)*pm)*SQ1M2(3) + 
        12*(873 - 1178*eta + 111*eta2 - 36*(121 + 15*eta)*pm)*SQ1M2(4) + 
        2*(7547 - 2692*eta + 292*eta2)*pm*SQ1M2(5) + 
        189*(-15 + eta)*eta*pm2*SQ1M2(6)))/1512. + 
        (EVALEJ(pm,pm,0,0)*(-18144*yp2 - 
        36*pm*(1 + SQ1M2(1))*(189*(7 - 2*eta) - 
        21*(-7 + 10*eta)*SQ1M2(1) + (-411 + 316*eta)*SQ1M2(2) + 
        (-243 + 148*eta)*SQ1M2(3)) - 
        pm3*SQ1M2(3)*(12*(-4415 + 148*eta + 42*eta2) + 
        9*(-6144 - 825*eta + 7*eta2)*SQ1M2(1) + 
        4*(6237 + 709*eta + 114*eta2)*SQ1M2(2) + 
        9*(-3440 + 1143*eta + 7*eta2)*SQ1M2(3) + 
        2*(295 - 951*eta + 149*eta2)*SQ1M2(4)) - 
        126*(-15 + eta)*eta*pm4*SQ1M2(6) + 
        3*pm2*(336*(-14 + 9*eta) + 84*(245 - 50*eta + 2*eta2)*SQ1M2(1) - 
        24*(-1927 + 419*eta)*SQ1M2(2) + 
        (4908 + 4840*eta - 592*eta2)*SQ1M2(3) + 
        (-12203 + 5643*eta - 503*eta2)*SQ1M2(4) - 
        2*(-1801 + 741*eta + 35*eta2)*SQ1M2(5) + 
        (-155 + 259*eta + 9*eta2)*SQ1M2(6))))/504. - 
        (12*(-1 + 3*eta)*pm2*EVALEJ(pm,pm,4,1)*SQ1M2(8))/7. + 
        (pm2*EVALEJ(pm,pm,4,0)*(112007 - 31888*eta + 16*eta2 - 
        2592*(-1 + 3*eta)*SQ1M2(1) - 
        3*(-1477 + 4772*eta + 1324*eta2)*SQ1M2(2) + 
        6*(-131 + 412*eta + 212*eta2)*pm*SQ1M2(3))*SQ1M2(8))/1512. + 
        ((-12793 + 23168*eta + 6664*eta2)*pm2*EVALEJ(pm,pm,5,0)*SQ1M2(10))/
        756. - (5*(-131 + 412*eta + 212*eta2)*pm2*EVALEJ(pm,pm,6,0)*
        SQ1M2(12))/252.;   
    PN4 /= ePow[2]*pm2*sq1m2Pow[4];

    ImPN5 = (-122*eta*pm*EVALEJ(pm,pm,0,0)*(1 + SQ1M2(2)))/35. + 
        (eta*EVALEJ(pm,pm,1,0)*(-732*SQ1M2(1) + 5*pm*(73 + 109*SQ1M2(2))))/
        105. + (eta*EVALEJ(pm,pm,2,0)*SQ1M2(1)*
        (7116 + 540*SQ1M2(2) + pm*SQ1M2(1)*(-1345 + 2547*SQ1M2(2))))/315.\
        - (2*eta*EVALEJ(pm,pm,3,0)*SQ1M2(1)*
        (1610 - 1301*SQ1M2(2) + 1205*pm*SQ1M2(3)))/105. + 
        (3398*eta*EVALEJ(pm,pm,5,0)*SQ1M2(5))/21. + 
        EVALEJ(pm,pm,4,0)*(-44*eta*SQ1M2(3) - (2074*eta*SQ1M2(5))/35. + 
        14*eta*pm*SQ1M2(6)) - 70*eta*EVALEJ(pm,pm,6,0)*SQ1M2(7);
    ImPN5 /= ePow[2];

    // PN6 = (-27*myp2*pm3*EPOW(3)*EVALEJ(-3 + pm,pm,0,1))/2. - 
    //     18*myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,3) - 
    //     18*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,3) - 
    //     (27*pm3*yp2*EPOW(3)*EVALEJ(3 + pm,pm,0,1))/2. + 
    //     54*pm3*yp2*EPOW(1)*EVALEJ(-1 + pm,pm,1,2)*SQ1M2(1) - 
    //     54*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,1,2)*SQ1M2(1) - 
    //     27*myp2*pm3*EVALEJ(-2 + pm,pm,0,2)*(-1 + SQ1M2(2)) + 
    //     27*pm3*yp2*EVALEJ(2 + pm,pm,0,2)*(-1 + SQ1M2(2)) - 
    //     72*pm3*EVALEJ(pm,pm,1,3)*SQ1M2(2) - 
    //     72*pm4*EVALEJ(pm,pm,0,3)*SQ1M2(3) - 
    //     (9*pm3*yp3*EVALEJ(2 + pm,pm,0,1)*(-1 + SQ1M2(1))*
    //     (-434 + 252*eta - 56*(-4 + eta)*SQ1M2(1) + 
    //     (200 - 110*eta)*SQ1M2(2) + 7*(-5 + 2*eta)*pm*SQ1M2(3)))/28. - 
    //     (9*myp3*pm3*EVALEJ(-2 + pm,pm,0,1)*(1 + SQ1M2(1))*
    //     (-434 + 252*eta + 56*(-4 + eta)*SQ1M2(1) + 
    //     (200 - 110*eta)*SQ1M2(2) + 7*(-5 + 2*eta)*pm*SQ1M2(3)))/28. + 
    //     (9*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,2)*
    //     (7*(-31 + 18*eta) - 28*(-4 + eta)*SQ1M2(1) + 
    //     (100 - 55*eta)*SQ1M2(2) + 7*(-5 + 2*eta)*pm*SQ1M2(3)))/14. + 
    //     (9*myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,2)*
    //     (7*(-31 + 18*eta) + 28*(-4 + eta)*SQ1M2(1) + 
    //     (100 - 55*eta)*SQ1M2(2) + 7*(-5 + 2*eta)*pm*SQ1M2(3)))/14. + 
    //     (6*pm3*EVALEJ(pm,pm,1,2)*SQ1M2(2)*
    //     (-645 + 297*eta - 252*SQ1M2(1) + (380 - 174*eta)*SQ1M2(2) + 
    //     (-229 + 36*eta)*pm*SQ1M2(3)))/7. + 
    //     72*pm3*EVALEJ(pm,pm,2,3)*SQ1M2(4) - 
    //     (6*pm3*EVALEJ(pm,pm,2,2)*(3*(37 + 99*eta) - 252*SQ1M2(1) + 
    //     (376 - 162*eta)*SQ1M2(2) + (-103 + 36*eta)*pm*SQ1M2(3))*SQ1M2(4))
    //     /7. - (3*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,1,1)*SQ1M2(1)*
    //     (168*(-31 + 18*eta) + 48*(-215 + 99*eta)*SQ1M2(1) + 
    //     24*(-61 + 15*eta)*SQ1M2(2) + 
    //     (6080 - 2784*eta + 7*(-60 + 9*eta + eta2)*pm)*SQ1M2(3) + 
    //     24*(100 - 55*eta + 7*(-5 + 2*eta)*pm)*SQ1M2(4) + 
    //     7*(-60 + 9*eta + eta2)*pm*SQ1M2(5)))/56. + 
    //     (3*pm3*EPOW(1)*EVALEJ(1 + pm,pm,1,1)*SQ1M2(1)*
    //     (168*(-31 + 18*eta) + 48*(-215 + 99*eta)*SQ1M2(1) + 
    //     24*(-61 + 15*eta)*SQ1M2(2) + 
    //     (6080 - 2784*eta + 7*(-60 + 9*eta + eta2)*pm)*SQ1M2(3) + 
    //     24*(100 - 55*eta + 7*(-5 + 2*eta)*pm)*SQ1M2(4) + 
    //     7*(-60 + 9*eta + eta2)*pm*SQ1M2(5)))/56. + 
    //     (pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,1)*
    //     (-21*(14304 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) + 
    //     4032*(269 - 109*eta + 20*eta2)*SQ1M2(1) + 
    //     3*(42480 + 40128*eta2 + eta*(-112816 + 861*CST_PISQ))*
    //     SQ1M2(2) - (576*(261 - 241*eta + 68*eta2) + 
    //     7*(2880 - 12512*eta + 2016*eta2 + 123*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) + 16*(-834 + 738*eta - 896*eta2 + 
    //     63*(80 - 37*eta + 7*eta2)*pm)*SQ1M2(4) + 
    //     144*(45 - 69*eta + 40*eta2)*pm*SQ1M2(5)))/2688. + 
    //     (myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,1)*
    //     (-21*(14304 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
    //     4032*(185 - 133*eta + 20*eta2)*SQ1M2(1) + 
    //     3*(-38160 + 40128*eta2 + eta*(-80560 + 861*CST_PISQ))*
    //     SQ1M2(2) + (576*(611 - 381*eta + 68*eta2) - 
    //     7*(2880 - 12512*eta + 2016*eta2 + 123*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) - 16*(834 - 738*eta + 896*eta2 + 
    //     63*(80 - 37*eta + 7*eta2)*pm)*SQ1M2(4) + 
    //     144*(45 - 69*eta + 40*eta2)*pm*SQ1M2(5)))/2688. + 
    //     (6*pm3*EVALEJ(pm,pm,0,2)*SQ1M2(1)*
    //     (126 + 252*SQ1M2(1) + 14*(9 + 7*(-4 + 3*eta)*pm)*SQ1M2(2) - 
    //     252*pm*SQ1M2(3) + (379 - 171*eta)*pm*SQ1M2(4) + 
    //     21*(-5 + 2*eta)*pm2*SQ1M2(5)))/7. + 
    //     (12*(121 + 15*eta)*pm3*EVALEJ(pm,pm,3,2)*SQ1M2(6))/7. - 
    //     (pm3*EVALEJ(pm,pm,3,1)*(-13727 + 118576*eta + 14852*eta2 - 
    //     148176*SQ1M2(1) - 9*(-12841 + 6476*eta + 512*eta2)*SQ1M2(2) + 
    //     2*(-13597 + 4418*eta + 1600*eta2)*pm*SQ1M2(3))*SQ1M2(6))/252. - 
    //     (pm3*EVALEJ(pm,pm,2,1)*SQ1M2(4)*
    //     (-9*(-465504 + 25560*eta2 + eta*(73784 + 4305*CST_PISQ)) + 
    //     1728*(779 + 498*eta)*SQ1M2(1) + 
    //     (-2351216 + 251264*eta2 + eta*(-26576 + 7749*CST_PISQ))*
    //     SQ1M2(2) + (-65664*(-27 + 11*eta) + 
    //     (928160 - 50128*eta + 640*eta2 - 2583*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) + 24*(-9064 + 2297*eta - 1557*eta2 + 
    //     72*(-101 + 30*eta)*pm)*SQ1M2(4) + 
    //     96*(790 - 563*eta + 71*eta2)*pm*SQ1M2(5) + 
    //     288*(5 - 17*eta + 6*eta2)*pm2*SQ1M2(6)))/2016. + 
    //     (pm3*EVALEJ(pm,pm,1,1)*SQ1M2(2)*
    //     (-9*(88896 + 26064*eta2 + 5*eta*(-56912 + 861*CST_PISQ)) + 
    //     48384*(-37 + 21*eta)*SQ1M2(1) + 
    //     3*(686432 + 98992*eta2 + eta*(-565888 + 2583*CST_PISQ))*
    //     SQ1M2(2) + (1728*(1387 - 472*eta) + 
    //     (-677728 + 703280*eta + 640*eta2 - 2583*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) + 48*(-2703 + 144*eta - 935*eta2 + 
    //     36*(-347 + 12*eta)*pm)*SQ1M2(4) + 
    //     32*(15667 - 7577*eta + 533*eta2)*pm*SQ1M2(5) + 
    //     72*(-1240 + 121*eta + 45*eta2)*pm2*SQ1M2(6)))/2016. - 
    //     (pm3*EVALEJ(pm,pm,3,0)*SQ1M2(6)*
    //     (5*(1133324608 - 91111840*eta2 + 6995200*eta3 - 
    //     3*eta*(37396952 + 5540535*CST_PISQ)) + 
    //     3300*(889120 + 47600*eta2 + eta*(354736 - 2583*CST_PISQ))*
    //     SQ1M2(1) + (-1806118452 + 360538600*eta2 + 6335800*eta3 + 
    //     eta*(169338860 - 57536325*CST_PISQ))*SQ1M2(2) + 
    //     40*(-660*(-80921 + 46651*eta + 5333*eta2) + 
    //     (-4621101 - 12940195*eta + 313145*eta2 + 502300*eta3 + 
    //         852390*eta*CST_PISQ)*pm)*SQ1M2(3) + 
    //     20*(-1867077 - 7052475*eta + 212170*eta2 + 823470*eta3 + 
    //     330*(-12073 + 9380*eta + 8284*eta2)*pm)*SQ1M2(4) - 
    //     40*(-224082 + 130225*eta + 597880*eta2 + 140180*eta3)*pm*
    //     SQ1M2(5) + 4950*eta*(251 - 787*eta + 109*eta2)*pm2*SQ1M2(6)))/
    //     3.3264e6 - (pm3*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*SQ1M2(1)*
    //     (-105*(16032 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
    //     30*(109056 + 26064*eta2 + eta*(-292624 + 4305*CST_PISQ))*
    //     SQ1M2(1) + 180*(5500 + 3708*eta2 + 
    //     eta*(5644 - 287*CST_PISQ))*SQ1M2(2) + 
    //     (10*(767072 - 590584*eta + 98488*eta2 + 2583*eta*CST_PISQ) + 
    //     3*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) + 5*(863472 + 27712*eta2 + 
    //     3*eta*(-143056 + 861*CST_PISQ))*SQ1M2(4) + 
    //     (-80*(2886 - 225*eta + 1933*eta2) + 
    //     3*(-89600 + 228164*eta - 38100*eta2 + 1520*eta3 + 
    //         1435*eta*CST_PISQ)*pm)*SQ1M2(5) - 
    //     160*(417 - 369*eta + 448*eta2)*SQ1M2(6) + 
    //     60*eta*(2053 - 701*eta + 20*eta2)*pm*SQ1M2(7)))/13440. + 
    //     (pm3*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*SQ1M2(1)*
    //     (-105*(16032 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
    //     30*(109056 + 26064*eta2 + eta*(-292624 + 4305*CST_PISQ))*
    //     SQ1M2(1) + 180*(5500 + 3708*eta2 + 
    //     eta*(5644 - 287*CST_PISQ))*SQ1M2(2) + 
    //     (10*(767072 - 590584*eta + 98488*eta2 + 2583*eta*CST_PISQ) + 
    //     3*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ)*pm)*
    //     SQ1M2(3) + 5*(863472 + 27712*eta2 + 
    //     3*eta*(-143056 + 861*CST_PISQ))*SQ1M2(4) + 
    //     (-80*(2886 - 225*eta + 1933*eta2) + 
    //     3*(-89600 + 228164*eta - 38100*eta2 + 1520*eta3 + 
    //         1435*eta*CST_PISQ)*pm)*SQ1M2(5) - 
    //     160*(417 - 369*eta + 448*eta2)*SQ1M2(6) + 
    //     60*eta*(2053 - 701*eta + 20*eta2)*pm*SQ1M2(7)))/13440. + 
    //     (pm3*EVALEJ(pm,pm,2,0)*SQ1M2(4)*
    //     (-2970*(-1744720 + 930120*eta2 + 
    //     3*eta*(-3140224 + 61705*CST_PISQ)) + 
    //     29700*(-580424 + 26024*eta2 + 21*eta*(16168 + 205*CST_PISQ))*
    //     SQ1M2(1) + 3*(4551646832 + 1083711100*eta2 + 7218800*eta3 - 
    //     15*eta*(456878756 + 142065*CST_PISQ))*SQ1M2(2) - 
    //     660*(15*(-2222808 + 364064*eta + 71232*eta2 + 
    //         4305*eta*CST_PISQ) + 
    //     2*(1720902 - 2202184*eta + 34640*eta2 + 5800*eta3 - 
    //         38745*eta*CST_PISQ)*pm)*SQ1M2(3) + 
    //     3*(-2157399888 + 204968960*eta + 105237400*eta2 + 8236600*eta3 - 
    //     6392925*eta*CST_PISQ + 
    //     4400*(-215471 + 35047*eta + 875*eta2)*pm)*SQ1M2(4) + 
    //     12*(19800*(5123 - 4022*eta + 822*eta2) + 
    //     (55960648 + 4346000*eta2 + 319600*eta3 + 
    //         75*eta*(-1064176 + 28413*CST_PISQ))*pm)*SQ1M2(5) + 
    //     100*(202734 + 458991*eta - 369951*eta2 + 41596*eta3 + 
    //     396*(-925 + 2828*eta + 2188*eta2)*pm + 
    //     198*eta*(-10539 + 1115*eta + 135*eta2)*pm2)*SQ1M2(6) - 
    //     600*(22107 - 93633*eta + 51152*eta2 + 4146*eta3)*pm*SQ1M2(7)))/
    //     9.9792e6 + (pm3*EVALEJ(pm,pm,0,1)*SQ1M2(1)*
    //     (-12096*(-35 + 16*eta) - 24192*(-29 + 16*eta)*SQ1M2(1) - 
    //     3*(960*(-17 + 9*eta) + 
    //     (-140496 - 161696*eta + 26544*eta2 + 4305*eta*CST_PISQ)*pm)
    //     *SQ1M2(2) + 864*(-456 + 276*eta + (-83 + 354*eta)*pm)*
    //     SQ1M2(3) + 3*(960*(-58 + 27*eta) + 
    //     (35952 - 135008*eta + 34704*eta2 + 861*eta*CST_PISQ)*pm)*
    //     SQ1M2(4) + pm*(288*(2169 - 914*eta) - 
    //     7*(-6000 - 8240*eta + 1008*eta2 + 123*eta*CST_PISQ)*pm)*
    //     SQ1M2(5) + 16*pm*(-4693 + 3272*eta - 1228*eta2 + 
    //     756*(-5 + 2*eta)*pm)*SQ1M2(6) + 
    //     48*(530 - 395*eta + 132*eta2)*pm2*SQ1M2(7)))/672. + 
    //     (EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*
    //     (862202880*SQ1M2(1) + 10264320*pm*
    //     (77 - 7*(-23 + 16*eta)*SQ1M2(1) + 77*SQ1M2(2) + 
    //     3*(-52 + 23*eta)*SQ1M2(3)) + 
    //     297*pm4*SQ1M2(4)*(3*(-89600 + 69504*eta + 560*eta2 + 
    //         1435*eta*CST_PISQ) + 20160*(-15 + eta)*eta*SQ1M2(1) + 
    //     (-268800 + 391892*eta - 67540*eta2 + 4840*eta3 + 
    //         4305*eta*CST_PISQ)*SQ1M2(2) + 
    //     20*eta*(-911 - 269*eta + 74*eta2)*SQ1M2(4)) + 
    //     2970*pm2*(-18144*(-31 + 18*eta) + 
    //     21*(6912 - 26080*eta + 5136*eta2 + 615*eta*CST_PISQ)*
    //     SQ1M2(1) - 864*(419 + 73*eta)*SQ1M2(2) - 
    //     3*(84624 - 148144*eta + 47040*eta2 + 861*eta*CST_PISQ)*
    //     SQ1M2(3) + 1440*(-250 + 127*eta)*SQ1M2(4) + 
    //     32*(4341 - 3168*eta + 1060*eta2)*SQ1M2(5)) + 
    //     pm3*(31185*(16032 + 7056*eta2 + 
    //         eta*(-37696 + 615*CST_PISQ)) - 
    //     297*(2877280 + 569520*eta2 - 6720*eta3 + 
    //         eta*(-4362048 + 73185*CST_PISQ))*SQ1M2(1) + 
    //     2970*(-1383192 - 84888*eta2 + 
    //         eta*(766184 + 4305*CST_PISQ))*SQ1M2(2) - 
    //     297*(-1354800 - 604960*eta2 + 19680*eta3 + 
    //         3*eta*(261584 + 7175*CST_PISQ))*SQ1M2(3) - 
    //     495*(-5131728 + 12288*eta2 + 
    //         eta*(1966544 + 9471*CST_PISQ))*SQ1M2(4) + 
    //     7920*(-37491 + 12897*eta - 6682*eta2 + 959*eta3)*SQ1M2(5) + 
    //     95040*(141 - 214*eta + 276*eta2)*SQ1M2(6) + 
    //     40*(119313 - 36855*eta + 21438*eta2 + 3185*eta3)*SQ1M2(7))))/
    //     1.99584e6 - (36*(-1 + 3*eta)*pm3*EVALEJ(pm,pm,4,2)*SQ1M2(8))/7. + 
    //     (pm3*EVALEJ(pm,pm,4,1)*(115787 - 45172*eta + 5848*eta2 - 
    //     2592*(-1 + 3*eta)*SQ1M2(1) - 
    //     3*(-865 + 2624*eta + 2260*eta2)*SQ1M2(2) + 
    //     6*(-41 + 106*eta + 320*eta2)*pm*SQ1M2(3))*SQ1M2(8))/252. - 
    //     (pm3*EVALEJ(pm,pm,4,0)*(5*(-142054480 + 149287700*eta2 + 
    //     3601020*eta3 + eta*(-626853244 + 25145505*CST_PISQ)) + 
    //     26400*(-98047 + 41927*eta + 4843*eta2)*SQ1M2(1) + 
    //     (-81297604 - 45616800*eta2 + 14787200*eta3 + 
    //     2025*eta*(-1117084 + 28413*CST_PISQ))*SQ1M2(2) - 
    //     2*(-79200*(-395 + 1222*eta + 812*eta2) + 
    //     (-73879604 - 183572560*eta - 24213760*eta2 + 4353200*eta3 + 
    //         6392925*eta*CST_PISQ)*pm)*SQ1M2(3) - 
    //     300*(84987 - 307600*eta + 39867*eta2 + 60625*eta3 + 
    //     264*(-131 + 412*eta + 212*eta2)*pm)*SQ1M2(4) + 
    //     600*(8410 - 35055*eta + 16980*eta2 + 8956*eta3)*pm*SQ1M2(5))*
    //     SQ1M2(8))/3.3264e6 + (pm3*EVALEJ(pm,pm,1,0)*SQ1M2(2)*
    //     (-14850*(-410368 + 74544*eta2 + 
    //     3*eta*(42256 + 1435*CST_PISQ)) - 
    //     14850*(-127152 + 51408*eta2 + eta*(-511520 + 11193*CST_PISQ))*
    //     SQ1M2(1) + 7425*(-879552 + 266176*eta2 + 
    //     91*eta*(-6592 + 123*CST_PISQ))*SQ1M2(2) + 
    //     165*(30*(2064656 + 230560*eta2 + 
    //         3*eta*(-716848 + 3731*CST_PISQ)) + 
    //     (27573920 - 3328840*eta2 + 46400*eta3 - 
    //         9*eta*(-2083432 + 73185*CST_PISQ))*pm)*SQ1M2(3) + 
    //     45*(52126488 + 50788392*eta - 19553480*eta2 + 217520*eta3 + 
    //     426195*eta*CST_PISQ + 
    //     440*(-23015 + 107182*eta + 13178*eta2)*pm)*SQ1M2(4) + 
    //     33*(-2400*(40447 - 19796*eta + 5084*eta2) + 
    //     (4203952 + 19461000*eta2 + 614400*eta3 - 
    //         25*eta*(3315272 + 7749*CST_PISQ))*pm)*SQ1M2(5) + 
    //     5*(8*(-902016 + 4528350*eta - 2318595*eta2 + 134015*eta3) - 
    //     11880*(-43483 + 17392*eta + 1884*eta2)*pm + 
    //     891*(-89600 + 168174*eta - 24870*eta2 + 750*eta3 + 
    //         1435*eta*CST_PISQ)*pm2)*SQ1M2(6) + 
    //     120*pm*(-449037 - 20990*eta3 + eta*(150665 - 935550*pm) + 
    //     90*eta2*(-1054 + 693*pm))*SQ1M2(7) + 
    //     14850*eta*(619 - 1111*eta + 233*eta2)*pm2*SQ1M2(8)))/9.9792e6 + 
    //     (EVALEJ(pm,pm,0,0)*(-431101440*yp2 - 
    //     5132160*pm*(1 + SQ1M2(1))*
    //     (371 - 126*eta + (175 - 98*eta)*SQ1M2(1) + 
    //     (-128 + 83*eta)*SQ1M2(2) + 5*(-20 + 11*eta)*SQ1M2(3)) - 
    //     1485*pm2*(1 + SQ1M2(1))*
    //     (21*(74208 + 7056*eta2 + eta*(-68800 + 615*CST_PISQ)) + 
    //     21*(-8544 + 3216*eta2 + eta*(-22528 + 615*CST_PISQ))*
    //     SQ1M2(1) - 3*(532464 + 67008*eta2 + 
    //         eta*(-273712 + 861*CST_PISQ))*SQ1M2(2) - 
    //     3*(-87120 + 27072*eta2 + eta*(-74032 + 861*CST_PISQ))*
    //     SQ1M2(3) + 32*(-9843 + 1989*eta + 1672*eta2)*SQ1M2(4) + 
    //     32*(417 - 369*eta + 448*eta2)*SQ1M2(5)) - 
    //     594*pm5*(-268800 + 4305*eta*CST_PISQ + 
    //     10*eta3*(176 + 271*SQ1M2(2) + 91*SQ1M2(4)) - 
    //     10*eta2*(440 - 1008*SQ1M2(1) - 355*SQ1M2(2) + 511*SQ1M2(4)) + 
    //     2*eta*(69806 - 75600*SQ1M2(1) - 31655*SQ1M2(2) + 805*SQ1M2(4))
    //     )*SQ1M2(6) - pm4*SQ1M2(3)*
    //     (-297*(4503680 + 359120*eta2 + 
    //         21*eta*(-405344 + 4305*CST_PISQ)) + 
    //     891*(-3670560 + 145280*eta2 + 
    //         eta*(528464 + 35875*CST_PISQ))*SQ1M2(1) - 
    //     495*(-2192688 - 91248*eta2 + 2592*eta3 + 
    //         41*eta*(48512 + 315*CST_PISQ))*SQ1M2(2) + 
    //     99*(11799680 - 1898060*eta2 + 14520*eta3 + 
    //         eta*(4536556 - 73185*CST_PISQ))*SQ1M2(3) + 
    //     264*(-3226761 + 1103085*eta + 46655*eta2 + 6600*eta3)*
    //     SQ1M2(4) + 5940*(38736 - 17215*eta + 7619*eta2 + 74*eta3)*
    //     SQ1M2(5) + 40*(100224 - 421200*eta + 169233*eta2 + 
    //         2465*eta3)*SQ1M2(6)) - 
    //     pm3*(20790*(17392 + 7056*eta2 + 
    //         5*eta*(-7712 + 123*CST_PISQ)) - 
    //     297*(10578400 + 569520*eta2 - 6720*eta3 + 
    //         eta*(-8797248 + 73185*CST_PISQ))*SQ1M2(1) - 
    //     594*(9713320 + 973920*eta2 + 
    //         eta*(-9706328 + 81795*CST_PISQ))*SQ1M2(2) - 
    //     99*(-15706800 - 1956000*eta2 + 59040*eta3 + 
    //         eta*(-3234544 + 107625*CST_PISQ))*SQ1M2(3) + 
    //     990*(4746096 + 469984*eta2 + 
    //         7*eta*(-431296 + 615*CST_PISQ))*SQ1M2(4) + 
    //     990*(832776 - 103808*eta2 + 7672*eta3 + 
    //         7*eta*(-95000 + 123*CST_PISQ))*SQ1M2(5) + 
    //     20*(-1587843 + 7855029*eta - 4460688*eta2 + 191483*eta3)*
    //     SQ1M2(6) + 40*(5036445 - 3667383*eta + 978966*eta2 + 
    //         3185*eta3)*SQ1M2(7) + 
    //     20*(-54135 - 262575*eta + 16488*eta2 + 2195*eta3)*SQ1M2(8))))/
    //     1.99584e6 + ((-12793 + 23168*eta + 6664*eta2)*pm3*EVALEJ(pm,pm,5,1)*
    //     SQ1M2(10))/126. + (pm3*EVALEJ(pm,pm,5,0)*
    //     (-98168468 + 41434420*eta2 + 15160500*eta3 + 
    //     eta*(-777794432 + 18326385*CST_PISQ) + 
    //     2640*(-13285 + 24440*eta + 10504*eta2)*SQ1M2(1) - 
    //     20*(-363402 + 72587*eta + 1400513*eta2 + 301453*eta3)*SQ1M2(2) + 
    //     8*(-253672 + 364525*eta + 535510*eta2 + 104650*eta3)*pm*SQ1M2(3))
    //     *SQ1M2(10))/332640. - (5*(-131 + 412*eta + 212*eta2)*pm3*
    //     EVALEJ(pm,pm,6,1)*SQ1M2(12))/42. - 
    //     (pm3*EVALEJ(pm,pm,6,0)*(-59360228 - 53380240*eta2 + 4296320*eta3 + 
    //     205*eta*(-1504444 + 31185*CST_PISQ) + 
    //     39600*(-131 + 412*eta + 212*eta2)*SQ1M2(1) - 
    //     600*(4401 - 19477*eta + 13313*eta2 + 6571*eta3)*SQ1M2(2) + 
    //     600*(56 - 557*eta + 1378*eta2 + 598*eta3)*pm*SQ1M2(3))*SQ1M2(12))
    //     /332640. - ((-165848 + 169110*eta + 567245*eta2 + 86665*eta3)*pm3*
    //     EVALEJ(pm,pm,7,0)*SQ1M2(14))/4620. + 
    //     (5*(56 - 557*eta + 1378*eta2 + 598*eta3)*pm3*EVALEJ(pm,pm,8,0)*
    //     SQ1M2(16))/396.;
    REAL8 PN6_lnv2x0_eta0 = (-214*pm3*SQ1M2(6)*(-5*EVALEJ(pm,pm,3,0) + 2*pm*EVALEJ(pm,pm,1,0)*SQ1M2(1) + 
        7*EVALEJ(pm,pm,3,0)*SQ1M2(2) - 25*EVALEJ(pm,pm,4,0)*SQ1M2(2) - 10*pm*EVALEJ(pm,pm,3,0)*SQ1M2(3) + 
        2*EVALEJ(pm,pm,2,0)*(1 + SQ1M2(2) + pm*SQ1M2(3)) - 21*EVALEJ(pm,pm,4,0)*SQ1M2(4) + 
        70*EVALEJ(pm,pm,5,0)*SQ1M2(4) + 6*pm*EVALEJ(pm,pm,4,0)*SQ1M2(5) - 30*EVALEJ(pm,pm,6,0)*SQ1M2(6)))/
        105.;
    PN6_lnv2x0_eta0 *= 1./(ePow[2]*pm3*sq1m2Pow[6]);

    PN6 = (-27*myp2*pm3*EPOW(3)*EVALEJ(-3 + pm,pm,0,1))/2. - 18*myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,3) - 
        18*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,3) - (27*pm3*yp2*EPOW(3)*EVALEJ(3 + pm,pm,0,1))/2. + 
        54*pm3*yp2*EPOW(1)*EVALEJ(-1 + pm,pm,1,2)*SQ1M2(1) - 
        54*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,1,2)*SQ1M2(1) - 
        27*myp2*pm3*EVALEJ(-2 + pm,pm,0,2)*(-1 + SQ1M2(2)) + 
        27*pm3*yp2*EVALEJ(2 + pm,pm,0,2)*(-1 + SQ1M2(2)) - 72*pm3*EVALEJ(pm,pm,1,3)*SQ1M2(2) - 
        72*pm4*EVALEJ(pm,pm,0,3)*SQ1M2(3) - (9*pm3*yp3*EVALEJ(2 + pm,pm,0,1)*(-1 + SQ1M2(1))*
        (-434 + 252*eta - 56*(-4 + eta)*SQ1M2(1) + (200 - 110*eta)*SQ1M2(2) + 
        7*(-5 + 2*eta)*pm*SQ1M2(3)))/28. - 
        (9*myp3*pm3*EVALEJ(-2 + pm,pm,0,1)*(1 + SQ1M2(1))*
        (-434 + 252*eta + 56*(-4 + eta)*SQ1M2(1) + (200 - 110*eta)*SQ1M2(2) + 
        7*(-5 + 2*eta)*pm*SQ1M2(3)))/28. + 
        (9*pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,2)*
        (7*(-31 + 18*eta) - 28*(-4 + eta)*SQ1M2(1) + (100 - 55*eta)*SQ1M2(2) + 
        7*(-5 + 2*eta)*pm*SQ1M2(3)))/14. + 
        (9*myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,2)*
        (7*(-31 + 18*eta) + 28*(-4 + eta)*SQ1M2(1) + (100 - 55*eta)*SQ1M2(2) + 
        7*(-5 + 2*eta)*pm*SQ1M2(3)))/14. + 
        (6*pm3*EVALEJ(pm,pm,1,2)*SQ1M2(2)*(-645 + 297*eta - 252*SQ1M2(1) + (380 - 174*eta)*SQ1M2(2) + 
        (-229 + 36*eta)*pm*SQ1M2(3)))/7. + 72*pm3*EVALEJ(pm,pm,2,3)*SQ1M2(4) + 
        (6*pm3*EVALEJ(pm,pm,2,2)*(-3*(37 + 99*eta) + 252*SQ1M2(1) + 2*(-188 + 81*eta)*SQ1M2(2) + 
        (103 - 36*eta)*pm*SQ1M2(3))*SQ1M2(4))/7. - 
        (3*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,1,1)*SQ1M2(1)*
        (168*(-31 + 18*eta) + 48*(-215 + 99*eta)*SQ1M2(1) + 24*(-61 + 15*eta)*SQ1M2(2) + 
        (6080 - 2784*eta + 7*(-60 + 9*eta + eta2)*pm)*SQ1M2(3) + 
        24*(100 - 55*eta + 7*(-5 + 2*eta)*pm)*SQ1M2(4) + 7*(-60 + 9*eta + eta2)*pm*SQ1M2(5)))/
        56. + (3*pm3*EPOW(1)*EVALEJ(1 + pm,pm,1,1)*SQ1M2(1)*
        (168*(-31 + 18*eta) + 48*(-215 + 99*eta)*SQ1M2(1) + 24*(-61 + 15*eta)*SQ1M2(2) + 
        (6080 - 2784*eta + 7*(-60 + 9*eta + eta2)*pm)*SQ1M2(3) + 
        24*(100 - 55*eta + 7*(-5 + 2*eta)*pm)*SQ1M2(4) + 7*(-60 + 9*eta + eta2)*pm*SQ1M2(5)))/
        56. + (pm3*yp2*EPOW(1)*EVALEJ(1 + pm,pm,0,1)*
        (-21*(14304 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) + 
        4032*(269 - 109*eta + 20*eta2)*SQ1M2(1) + 
        3*(42480 + 40128*eta2 + eta*(-112816 + 861*CST_PISQ))*SQ1M2(2) - 
        (576*(261 - 241*eta + 68*eta2) + 
        7*(2880 - 12512*eta + 2016*eta2 + 123*eta*CST_PISQ)*pm)*SQ1M2(3) + 
        16*(-834 + 738*eta - 896*eta2 + 63*(80 - 37*eta + 7*eta2)*pm)*SQ1M2(4) + 
        144*(45 - 69*eta + 40*eta2)*pm*SQ1M2(5)))/2688. + 
        (myp2*pm3*EPOW(1)*EVALEJ(-1 + pm,pm,0,1)*
        (-21*(14304 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
        4032*(185 - 133*eta + 20*eta2)*SQ1M2(1) + 
        3*(-38160 + 40128*eta2 + eta*(-80560 + 861*CST_PISQ))*SQ1M2(2) + 
        (576*(611 - 381*eta + 68*eta2) - 
        7*(2880 - 12512*eta + 2016*eta2 + 123*eta*CST_PISQ)*pm)*SQ1M2(3) - 
        16*(834 - 738*eta + 896*eta2 + 63*(80 - 37*eta + 7*eta2)*pm)*SQ1M2(4) + 
        144*(45 - 69*eta + 40*eta2)*pm*SQ1M2(5)))/2688. + 
        (6*pm3*EVALEJ(pm,pm,0,2)*SQ1M2(1)*(126 + 252*SQ1M2(1) + 14*(9 + 7*(-4 + 3*eta)*pm)*SQ1M2(2) - 
        252*pm*SQ1M2(3) + (379 - 171*eta)*pm*SQ1M2(4) + 21*(-5 + 2*eta)*pm2*SQ1M2(5)))/7. + 
        (12*(121 + 15*eta)*pm3*EVALEJ(pm,pm,3,2)*SQ1M2(6))/7. + 
        (428*pm3*EVALEK(pm,pm,2,0)*(1 + SQ1M2(2) + pm*SQ1M2(3))*SQ1M2(6))/105. - 
        (214*pm3*EVALEK(pm,pm,3,0)*(5 - 7*SQ1M2(2) + 10*pm*SQ1M2(3))*SQ1M2(6))/105. + 
        (pm3*EVALEJ(pm,pm,3,1)*(13727 - 118576*eta - 14852*eta2 + 148176*SQ1M2(1) + 
        9*(-12841 + 6476*eta + 512*eta2)*SQ1M2(2) - 2*(-13597 + 4418*eta + 1600*eta2)*pm*SQ1M2(3)
        )*SQ1M2(6))/252. + (pm3*EVALEJ(pm,pm,2,1)*SQ1M2(4)*
        (9*(-465504 + 25560*eta2 + eta*(73784 + 4305*CST_PISQ)) - 
        1728*(779 + 498*eta)*SQ1M2(1) + 
        (2351216 - 251264*eta2 + eta*(26576 - 7749*CST_PISQ))*SQ1M2(2) + 
        (65664*(-27 + 11*eta) + (-928160 + 50128*eta - 640*eta2 + 2583*eta*CST_PISQ)*pm)*
        SQ1M2(3) + 24*(9064 - 2297*eta + 1557*eta2 - 72*(-101 + 30*eta)*pm)*SQ1M2(4) - 
        96*(790 - 563*eta + 71*eta2)*pm*SQ1M2(5) - 288*(5 - 17*eta + 6*eta2)*pm2*SQ1M2(6)))/2016.
        + (pm3*EVALEJ(pm,pm,1,1)*SQ1M2(2)*(-9*
        (88896 + 26064*eta2 + 5*eta*(-56912 + 861*CST_PISQ)) + 
        48384*(-37 + 21*eta)*SQ1M2(1) + 
        3*(686432 + 98992*eta2 + eta*(-565888 + 2583*CST_PISQ))*SQ1M2(2) + 
        (1728*(1387 - 472*eta) + (-677728 + 703280*eta + 640*eta2 - 2583*eta*CST_PISQ)*pm)*
        SQ1M2(3) + 48*(-2703 + 144*eta - 935*eta2 + 36*(-347 + 12*eta)*pm)*SQ1M2(4) + 
        32*(15667 - 7577*eta + 533*eta2)*pm*SQ1M2(5) + 
        72*(-1240 + 121*eta + 45*eta2)*pm2*SQ1M2(6)))/2016. - 
        (pm3*EVALEJ(pm,pm,3,0)*SQ1M2(6)*(5*(1133324608 - 91111840*eta2 + 6995200*eta3 - 
        3*eta*(37396952 + 5540535*CST_PISQ)) + 
        3300*(889120 + 47600*eta2 + eta*(354736 - 2583*CST_PISQ))*SQ1M2(1) + 
        (-1806118452 + 360538600*eta2 + 6335800*eta3 + eta*(169338860 - 57536325*CST_PISQ))*
        SQ1M2(2) + 40*(-660*(-80921 + 46651*eta + 5333*eta2) + 
        (-4197381 - 12940195*eta + 313145*eta2 + 502300*eta3 + 852390*eta*CST_PISQ)*pm)*
        SQ1M2(3) + 20*(-1867077 - 7052475*eta + 212170*eta2 + 823470*eta3 + 
        330*(-12073 + 9380*eta + 8284*eta2)*pm)*SQ1M2(4) - 
        40*(-224082 + 130225*eta + 597880*eta2 + 140180*eta3)*pm*SQ1M2(5) + 
        4950*eta*(251 - 787*eta + 109*eta2)*pm2*SQ1M2(6)))/3.3264e6 + 
        (428*pm4*EVALEK(pm,pm,1,0)*SQ1M2(7))/105. - 
        (pm3*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*SQ1M2(1)*
        (-105*(16032 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
        30*(109056 + 26064*eta2 + eta*(-292624 + 4305*CST_PISQ))*SQ1M2(1) + 
        180*(5500 + 3708*eta2 + eta*(5644 - 287*CST_PISQ))*SQ1M2(2) + 
        (10*(767072 - 590584*eta + 98488*eta2 + 2583*eta*CST_PISQ) + 
        3*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ)*pm)*SQ1M2(3) + 
        5*(863472 + 27712*eta2 + 3*eta*(-143056 + 861*CST_PISQ))*SQ1M2(4) + 
        (-80*(2886 - 225*eta + 1933*eta2) + 
        3*(-89600 + 228164*eta - 38100*eta2 + 1520*eta3 + 1435*eta*CST_PISQ)*pm)*SQ1M2(5)\
        - 160*(417 - 369*eta + 448*eta2)*SQ1M2(6) + 
        60*eta*(2053 - 701*eta + 20*eta2)*pm*SQ1M2(7)))/13440. + 
        (pm3*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*SQ1M2(1)*
        (-105*(16032 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
        30*(109056 + 26064*eta2 + eta*(-292624 + 4305*CST_PISQ))*SQ1M2(1) + 
        180*(5500 + 3708*eta2 + eta*(5644 - 287*CST_PISQ))*SQ1M2(2) + 
        (10*(767072 - 590584*eta + 98488*eta2 + 2583*eta*CST_PISQ) + 
        3*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ)*pm)*SQ1M2(3) + 
        5*(863472 + 27712*eta2 + 3*eta*(-143056 + 861*CST_PISQ))*SQ1M2(4) + 
        (-80*(2886 - 225*eta + 1933*eta2) + 
        3*(-89600 + 228164*eta - 38100*eta2 + 1520*eta3 + 1435*eta*CST_PISQ)*pm)*SQ1M2(5)\
        - 160*(417 - 369*eta + 448*eta2)*SQ1M2(6) + 
        60*eta*(2053 - 701*eta + 20*eta2)*pm*SQ1M2(7)))/13440. + 
        (pm3*EVALEJ(pm,pm,2,0)*SQ1M2(4)*(-2970*
        (-1744720 + 930120*eta2 + 3*eta*(-3140224 + 61705*CST_PISQ)) + 
        29700*(-580424 + 26024*eta2 + 21*eta*(16168 + 205*CST_PISQ))*SQ1M2(1) + 
        3*(4551646832 + 1083711100*eta2 + 7218800*eta3 - 15*eta*(456878756 + 142065*CST_PISQ))*
        SQ1M2(2) - 660*(15*(-2222808 + 364064*eta + 71232*eta2 + 4305*eta*CST_PISQ) + 
        2*(1720902 - 2202184*eta + 34640*eta2 + 5800*eta3 - 38745*eta*CST_PISQ)*pm)*
        SQ1M2(3) + 3*(-2157399888 + 204968960*eta + 105237400*eta2 + 8236600*eta3 - 
        6392925*eta*CST_PISQ + 4400*(-215471 + 35047*eta + 875*eta2)*pm)*SQ1M2(4) + 
        12*(19800*(5123 - 4022*eta + 822*eta2) + 
        (57090568 + 4346000*eta2 + 319600*eta3 + 75*eta*(-1064176 + 28413*CST_PISQ))*pm)*
        SQ1M2(5) + 100*(202734 + 458991*eta - 369951*eta2 + 41596*eta3 + 
        396*(-925 + 2828*eta + 2188*eta2)*pm + 198*eta*(-10539 + 1115*eta + 135*eta2)*pm2)*
        SQ1M2(6) - 600*(22107 - 93633*eta + 51152*eta2 + 4146*eta3)*pm*SQ1M2(7)))/9.9792e6 + 
        (pm3*EVALEJ(pm,pm,0,1)*SQ1M2(1)*(-12096*(-35 + 16*eta) - 24192*(-29 + 16*eta)*SQ1M2(1) - 
        3*(960*(-17 + 9*eta) + (-140496 - 161696*eta + 26544*eta2 + 4305*eta*CST_PISQ)*pm)*
        SQ1M2(2) + 864*(-456 + 276*eta + (-83 + 354*eta)*pm)*SQ1M2(3) + 
        3*(960*(-58 + 27*eta) + (35952 - 135008*eta + 34704*eta2 + 861*eta*CST_PISQ)*pm)*
        SQ1M2(4) + pm*(288*(2169 - 914*eta) - 
        7*(-6000 - 8240*eta + 1008*eta2 + 123*eta*CST_PISQ)*pm)*SQ1M2(5) + 
        16*pm*(-4693 + 3272*eta - 1228*eta2 + 756*(-5 + 2*eta)*pm)*SQ1M2(6) + 
        48*(530 - 395*eta + 132*eta2)*pm2*SQ1M2(7)))/672. + 
        (EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*(862202880*SQ1M2(1) + 
        10264320*pm*(77 - 7*(-23 + 16*eta)*SQ1M2(1) + 77*SQ1M2(2) + 3*(-52 + 23*eta)*SQ1M2(3)) + 
        297*pm4*SQ1M2(4)*(3*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ) + 
        20160*(-15 + eta)*eta*SQ1M2(1) + 
        (-268800 + 391892*eta - 67540*eta2 + 4840*eta3 + 4305*eta*CST_PISQ)*SQ1M2(2) + 
        20*eta*(-911 - 269*eta + 74*eta2)*SQ1M2(4)) + 
        2970*pm2*(-18144*(-31 + 18*eta) + 
        21*(6912 - 26080*eta + 5136*eta2 + 615*eta*CST_PISQ)*SQ1M2(1) - 
        864*(419 + 73*eta)*SQ1M2(2) - 
        3*(84624 - 148144*eta + 47040*eta2 + 861*eta*CST_PISQ)*SQ1M2(3) + 
        1440*(-250 + 127*eta)*SQ1M2(4) + 32*(4341 - 3168*eta + 1060*eta2)*SQ1M2(5)) + 
        pm3*(31185*(16032 + 7056*eta2 + eta*(-37696 + 615*CST_PISQ)) - 
        297*(2877280 + 569520*eta2 - 6720*eta3 + eta*(-4362048 + 73185*CST_PISQ))*
        SQ1M2(1) + 2970*(-1383192 - 84888*eta2 + eta*(766184 + 4305*CST_PISQ))*SQ1M2(2) - 
        297*(-1354800 - 604960*eta2 + 19680*eta3 + 3*eta*(261584 + 7175*CST_PISQ))*
        SQ1M2(3) - 495*(-5131728 + 12288*eta2 + eta*(1966544 + 9471*CST_PISQ))*SQ1M2(4) + 
        7920*(-37491 + 12897*eta - 6682*eta2 + 959*eta3)*SQ1M2(5) + 
        95040*(141 - 214*eta + 276*eta2)*SQ1M2(6) + 
        40*(119313 - 36855*eta + 21438*eta2 + 3185*eta3)*SQ1M2(7))))/1.99584e6 - 
        (36*(-1 + 3*eta)*pm3*EVALEJ(pm,pm,4,2)*SQ1M2(8))/7. + 
        (214*pm3*EVALEK(pm,pm,4,0)*(-25 - 21*SQ1M2(2) + 6*pm*SQ1M2(3))*SQ1M2(8))/105. + 
        (pm3*EVALEJ(pm,pm,4,1)*(115787 - 45172*eta + 5848*eta2 - 2592*(-1 + 3*eta)*SQ1M2(1) - 
        3*(-865 + 2624*eta + 2260*eta2)*SQ1M2(2) + 6*(-41 + 106*eta + 320*eta2)*pm*SQ1M2(3))*
        SQ1M2(8))/252. - (pm3*EVALEJ(pm,pm,4,0)*
        (5*(-142054480 + 149287700*eta2 + 3601020*eta3 + eta*(-626853244 + 25145505*CST_PISQ)) + 
        26400*(-98047 + 41927*eta + 4843*eta2)*SQ1M2(1) + 
        (-81297604 - 45616800*eta2 + 14787200*eta3 + 2025*eta*(-1117084 + 28413*CST_PISQ))*
        SQ1M2(2) - 2*(-79200*(-395 + 1222*eta + 812*eta2) + 
        (-69811892 - 183572560*eta - 24213760*eta2 + 4353200*eta3 + 6392925*eta*CST_PISQ)*
        pm)*SQ1M2(3) - 300*(84987 - 307600*eta + 39867*eta2 + 60625*eta3 + 
        264*(-131 + 412*eta + 212*eta2)*pm)*SQ1M2(4) + 
        600*(8410 - 35055*eta + 16980*eta2 + 8956*eta3)*pm*SQ1M2(5))*SQ1M2(8))/3.3264e6 + 
        (pm3*EVALEJ(pm,pm,1,0)*SQ1M2(2)*(-14850*
        (-410368 + 74544*eta2 + 3*eta*(42256 + 1435*CST_PISQ)) - 
        14850*(-127152 + 51408*eta2 + eta*(-511520 + 11193*CST_PISQ))*SQ1M2(1) + 
        7425*(-879552 + 266176*eta2 + 91*eta*(-6592 + 123*CST_PISQ))*SQ1M2(2) + 
        165*(30*(2064656 + 230560*eta2 + 3*eta*(-716848 + 3731*CST_PISQ)) + 
        (27573920 - 3328840*eta2 + 46400*eta3 - 9*eta*(-2083432 + 73185*CST_PISQ))*pm)*
        SQ1M2(3) + 45*(52126488 + 50788392*eta - 19553480*eta2 + 217520*eta3 + 
        426195*eta*CST_PISQ + 440*(-23015 + 107182*eta + 13178*eta2)*pm)*SQ1M2(4) + 
        33*(-2400*(40447 - 19796*eta + 5084*eta2) + 
        (4820272 + 19461000*eta2 + 614400*eta3 - 25*eta*(3315272 + 7749*CST_PISQ))*pm)*
        SQ1M2(5) + 5*(8*(-902016 + 4528350*eta - 2318595*eta2 + 134015*eta3) - 
        11880*(-43483 + 17392*eta + 1884*eta2)*pm + 
        891*(-89600 + 168174*eta - 24870*eta2 + 750*eta3 + 1435*eta*CST_PISQ)*pm2)*SQ1M2(6)
        + 120*pm*(-449037 - 20990*eta3 + eta*(150665 - 935550*pm) + 90*eta2*(-1054 + 693*pm))*
        SQ1M2(7) + 14850*eta*(619 - 1111*eta + 233*eta2)*pm2*SQ1M2(8)))/9.9792e6 + 
        (EVALEJ(pm,pm,0,0)*(-431101440*yp2 - 
        5132160*pm*(1 + SQ1M2(1))*(371 - 126*eta + (175 - 98*eta)*SQ1M2(1) + 
        (-128 + 83*eta)*SQ1M2(2) + 5*(-20 + 11*eta)*SQ1M2(3)) - 
        1485*pm2*(1 + SQ1M2(1))*(21*(74208 + 7056*eta2 + eta*(-68800 + 615*CST_PISQ)) + 
        21*(-8544 + 3216*eta2 + eta*(-22528 + 615*CST_PISQ))*SQ1M2(1) - 
        3*(532464 + 67008*eta2 + eta*(-273712 + 861*CST_PISQ))*SQ1M2(2) - 
        3*(-87120 + 27072*eta2 + eta*(-74032 + 861*CST_PISQ))*SQ1M2(3) + 
        32*(-9843 + 1989*eta + 1672*eta2)*SQ1M2(4) + 32*(417 - 369*eta + 448*eta2)*SQ1M2(5))\
        - 594*pm5*(-268800 + 4305*eta*CST_PISQ + 10*eta3*(176 + 271*SQ1M2(2) + 91*SQ1M2(4)) - 
        10*eta2*(440 - 1008*SQ1M2(1) - 355*SQ1M2(2) + 511*SQ1M2(4)) + 
        2*eta*(69806 - 75600*SQ1M2(1) - 31655*SQ1M2(2) + 805*SQ1M2(4)))*SQ1M2(6) - 
        pm4*SQ1M2(3)*(-297*(4503680 + 359120*eta2 + 21*eta*(-405344 + 4305*CST_PISQ)) + 
        891*(-3670560 + 145280*eta2 + eta*(528464 + 35875*CST_PISQ))*SQ1M2(1) - 
        495*(-2192688 - 91248*eta2 + 2592*eta3 + 41*eta*(48512 + 315*CST_PISQ))*SQ1M2(2) + 
        99*(11799680 - 1898060*eta2 + 14520*eta3 + eta*(4536556 - 73185*CST_PISQ))*
        SQ1M2(3) + 264*(-3226761 + 1103085*eta + 46655*eta2 + 6600*eta3)*SQ1M2(4) + 
        5940*(38736 - 17215*eta + 7619*eta2 + 74*eta3)*SQ1M2(5) + 
        40*(100224 - 421200*eta + 169233*eta2 + 2465*eta3)*SQ1M2(6)) - 
        pm3*(20790*(17392 + 7056*eta2 + 5*eta*(-7712 + 123*CST_PISQ)) - 
        297*(10578400 + 569520*eta2 - 6720*eta3 + eta*(-8797248 + 73185*CST_PISQ))*
        SQ1M2(1) - 594*(9713320 + 973920*eta2 + eta*(-9706328 + 81795*CST_PISQ))*
        SQ1M2(2) - 99*(-15706800 - 1956000*eta2 + 59040*eta3 + 
            eta*(-3234544 + 107625*CST_PISQ))*SQ1M2(3) + 
        990*(4746096 + 469984*eta2 + 7*eta*(-431296 + 615*CST_PISQ))*SQ1M2(4) + 
        990*(832776 - 103808*eta2 + 7672*eta3 + 7*eta*(-95000 + 123*CST_PISQ))*SQ1M2(5) + 
        20*(-1587843 + 7855029*eta - 4460688*eta2 + 191483*eta3)*SQ1M2(6) + 
        40*(5036445 - 3667383*eta + 978966*eta2 + 3185*eta3)*SQ1M2(7) + 
        20*(-54135 - 262575*eta + 16488*eta2 + 2195*eta3)*SQ1M2(8))))/1.99584e6 + 
        ((-12793 + 23168*eta + 6664*eta2)*pm3*EVALEJ(pm,pm,5,1)*SQ1M2(10))/126. + 
        (428*pm3*EVALEK(pm,pm,5,0)*SQ1M2(10))/3. + 
        (pm3*EVALEJ(pm,pm,5,0)*(-98168468 + 41434420*eta2 + 15160500*eta3 + 
        eta*(-777794432 + 18326385*CST_PISQ) + 
        2640*(-13285 + 24440*eta + 10504*eta2)*SQ1M2(1) - 
        20*(-363402 + 72587*eta + 1400513*eta2 + 301453*eta3)*SQ1M2(2) + 
        8*(-253672 + 364525*eta + 535510*eta2 + 104650*eta3)*pm*SQ1M2(3))*SQ1M2(10))/332640. - 
        (5*(-131 + 412*eta + 212*eta2)*pm3*EVALEJ(pm,pm,6,1)*SQ1M2(12))/42. - 
        (428*pm3*EVALEK(pm,pm,6,0)*SQ1M2(12))/7. + 
        (pm3*EVALEJ(pm,pm,6,0)*(59360228 + 53380240*eta2 - 4296320*eta3 - 
        205*eta*(-1504444 + 31185*CST_PISQ) - 39600*(-131 + 412*eta + 212*eta2)*SQ1M2(1) + 
        600*(4401 - 19477*eta + 13313*eta2 + 6571*eta3)*SQ1M2(2) - 
        600*(56 - 557*eta + 1378*eta2 + 598*eta3)*pm*SQ1M2(3))*SQ1M2(12))/332640. - 
        ((-165848 + 169110*eta + 567245*eta2 + 86665*eta3)*pm3*EVALEJ(pm,pm,7,0)*SQ1M2(14))/4620. + 
        (5*(56 - 557*eta + 1378*eta2 + 598*eta3)*pm3*EVALEJ(pm,pm,8,0)*SQ1M2(16))/396.;
    PN6 /= ePow[2]*pm3*SQ1M2(6);
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];
    PN6 = PN6 + PN6_lnv2x0_eta0*log(v2/CONST_X0);

    // tail part
    COMPLEX16 PN3_tail, PN5_tail, PN6_tail;
    COMPLEX16 lnIpvO2 =  pm>0 ? (log(0.5*pm) + 0.5*I*CST_PI) : (log(-0.5*pm) - 0.5*I*CST_PI);
    REAL8 lnxOx0p = log(v2/CONST_X0P);
    COMPLEX16 pref_tail, pref_tail_PN6;
    pref_tail = 3.*lnxOx0p + 2*lnIpvO2;
    PN3_tail = pm*EVALEJ(pm,pm,1,0)*SQ1M2(1) + 
        EVALEJ(pm,pm,3,0)*(-2.5 + (7*SQ1M2(2))/2. - 5*pm*SQ1M2(3)) + 
        EVALEJ(pm,pm,2,0)*(1 + SQ1M2(2) + pm*SQ1M2(3)) + 
        (EVALEJ(pm,pm,4,0)*SQ1M2(2)*(-25 - 21*SQ1M2(2) + 6*pm*SQ1M2(3)))/2. + 
        35*EVALEJ(pm,pm,5,0)*SQ1M2(4) - 15*EVALEJ(pm,pm,6,0)*SQ1M2(6);
    PN3_tail *= pref_tail * I / (ePow[2]*pm);

    PN5_tail = pref_tail * (
        6*pm2*EVALEJ(pm,pm,1,1)*SQ1M2(1) - 
        6*pm2*EVALEJ(pm,pm,0,0)*(1 + SQ1M2(2)) + 
        6*pm*EVALEJ(pm,pm,2,1)*(1 + SQ1M2(2) + pm*SQ1M2(3)) + 
        3*pm*EVALEJ(pm,pm,4,1)*SQ1M2(2)*(-25 - 21*SQ1M2(2) + 6*pm*SQ1M2(3)) - 
        3*pm*EVALEJ(pm,pm,3,1)*(5 - 7*SQ1M2(2) + 10*pm*SQ1M2(3)) + 
        EVALEJ(pm,pm,1,0)*(-6*pm*SQ1M2(1) + 
        (pm2*(315 + 42*(-9 + eta)*SQ1M2(1) - 441*SQ1M2(2) + 
        (-111 + 4*eta)*SQ1M2(3)))/42.) + 
        210*pm*EVALEJ(pm,pm,5,1)*SQ1M2(4) + 
        (EVALEJ(pm,pm,5,0)*(1260 + 
        pm*(135*(-20 + 3*eta) + 1692*SQ1M2(1) - 
        (791 + 277*eta)*SQ1M2(2)) + 2*(67 + 15*eta)*pm2*SQ1M2(3))*
        SQ1M2(4))/6. + (EVALEJ(pm,pm,2,0)*
        (252*(1 + SQ1M2(2)) + pm2*SQ1M2(1)*
        (210 + 1050*SQ1M2(1) + 14*(-25 + 9*eta)*SQ1M2(2) + 
        882*SQ1M2(3) + (67 - 68*eta)*SQ1M2(4)) + 
        pm*(-210 + 756*SQ1M2(1) + 5*(-165 + 26*eta)*SQ1M2(2) - 
        504*SQ1M2(3) + (57 - 38*eta)*SQ1M2(4))))/42. + 
        EVALEJ(pm,pm,3,0)*(3*(-5 + 7*SQ1M2(2)) + 
        (pm2*(245 - 60*eta - 315*SQ1M2(1) + (86 + 28*eta)*SQ1M2(2))*
        SQ1M2(3))/6. + (pm*(1050 + 7*(219 + 46*eta)*SQ1M2(2) + 
        9072*SQ1M2(3) - (561 + 158*eta)*SQ1M2(4)))/84.) - 
        90*pm*EVALEJ(pm,pm,6,1)*SQ1M2(6) - 
        (EVALEJ(pm,pm,6,0)*(3780 - 
        7*pm*(2884 + 15*eta - 540*SQ1M2(1) + 15*(1 + 11*eta)*SQ1M2(2)) + 
        30*(-1 + 3*eta)*pm2*SQ1M2(3))*SQ1M2(6))/42. + 
        EVALEJ(pm,pm,4,0)*((-3*(50 + (-43 + 18*eta)*pm)*SQ1M2(2))/2. - 
        210*pm*SQ1M2(3) + (-63 + (159.91666666666666 - (53*eta)/6.)*pm)*
        SQ1M2(4) + pm*(-90 + (-65.23333333333333 + 4*eta)*pm)*SQ1M2(5) + 
        (pm*(209 + 290*eta + 504*pm)*SQ1M2(6))/28. - 
        (1 + 4*eta)*pm2*SQ1M2(7)) - 
        3*(43 + 15*eta)*pm*EVALEJ(pm,pm,7,0)*SQ1M2(8) + 
        5*(-1 + 3*eta)*pm*EVALEJ(pm,pm,8,0)*SQ1M2(10)
    ) + (
        -6*pm*(2 + pm)*EVALEJ(pm,pm,1,0)*SQ1M2(1) - 
        6*(2 + pm)*EVALEJ(pm,pm,2,0)*(1 + SQ1M2(2) + pm*SQ1M2(3)) + 
        EVALEJ(pm,pm,3,0)*(15*(2 + pm) - 21*(2 + pm)*SQ1M2(2) + 
        (pm*(1680 + pm*(840 - 1804*CST_LN2 + 713*(2*CST_LN2) + 63*(6*CST_LN2)))*
        SQ1M2(3))/28.) + (EVALEJ(pm,pm,5,0)*
        (-35280 + pm*(29562*CST_LN2 - 14886*(2*CST_LN2) + 35*(-504 + (6*CST_LN2))))*
        SQ1M2(4))/84. + EVALEJ(pm,pm,4,0)*
        (75*(2 + pm)*SQ1M2(2) + 63*(2 + pm)*SQ1M2(4) + 
        (pm*(-1260 + pm*(-630 + 1006*CST_LN2 - 524*(2*CST_LN2) + 7*(6*CST_LN2)))*
        SQ1M2(5))/35.) + (EVALEJ(pm,pm,6,0)*
        (7560 + pm*(3780 - 6498*CST_LN2 + 3144*(2*CST_LN2) + 35*(6*CST_LN2)))*
        SQ1M2(6))/42.
    );
    PN5_tail *= I / (ePow[2]*pm2*sq1m2Pow[2]);

    pref_tail_PN6 = 515063 + 396900*lnxOx0p*lnxOx0p + 179760*lnIpvO2 + 
        176400*lnIpvO2*lnIpvO2 + 2520*lnxOx0p*(107 + 210*lnIpvO2) + 
        29400*CST_PISQ;
    PN6_tail = (pm*EVALEJ(pm,pm,2,0)*(1 + SQ1M2(2)))/88200. + 
        (EVALEJ(pm,pm,3,0)*(12*SQ1M2(1) + pm*(-5 + 7*SQ1M2(2))))/176400. - 
        (EVALEJ(pm,pm,4,0)*SQ1M2(1)*
        (28 - 24*SQ1M2(2) + pm*SQ1M2(1)*(25 + 21*SQ1M2(2))))/176400. + 
        (EVALEJ(pm,pm,5,0)*(-17 + 5*pm*SQ1M2(1))*SQ1M2(3))/12600. - 
        (EVALEJ(pm,pm,7,0)*SQ1M2(5))/280. + 
        (EVALEJ(pm,pm,6,0)*(44*SQ1M2(3) + 21*SQ1M2(5) - 3*pm*SQ1M2(6)))/
        17640. + (EVALEJ(pm,pm,8,0)*SQ1M2(7))/840.;
    PN6_tail *= pref_tail_PN6 / (ePow[2]*pm);
    // print_debug("pm = %d\n", pm);
    // print_debug("pref_tail_PN6 = %.16e + i %.16e\n", creal(pref_tail_PN6), cimag(pref_tail_PN6));
    // print_debug("PN0 = %.16e\n", PN0);
    // print_debug("PN2 = %.16e\n", dc->vPow[2]*PN2);
    // print_debug("PN3_tail = %.16e, %.16e\n", dc->vPow[3]*creal(PN3_tail), dc->vPow[3]*cimag(PN3_tail) );
    // print_debug("PN4 = %.16e\n", dc->vPow[4]*PN4);
    // print_debug("ImPN5 = %.16e\n", dc->vPow[5]*ImPN5);
    // print_debug("PN5_tail = %.16e, %.16e\n", dc->vPow[5]*creal(PN5_tail), dc->vPow[5]*cimag(PN5_tail) );
    // print_debug("PN6 = %.16e\n", dc->vPow[6]*PN6);
    // print_debug("PN6_tail = %.16e, %.16e\n", dc->vPow[6]*creal(PN6_tail), dc->vPow[6]*cimag(PN6_tail) );
    // return v2*PN0;
    return (PN0 + 
        v2 * (PN2 + 
            v*(PN3_tail + 
                v*(PN4 + 
                    v*( I*ImPN5 + PN5_tail +
                        v*( PN6 + PN6_tail) )))));
    // return (PN0 + 
    //     v2 * (PN2 + 
    //         v*(PN3_tail + 
    //             v*(PN4 + 
    //                 v*( I*ImPN5 + PN5_tail
    //                 )))));
#undef EVALEJ
#undef EVALEK
#undef SQ1M2
#undef EPOW
}

COMPLEX16 calculate_h22FourierMode_zeroe(REAL8 eta, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];

    REAL8 PN0 = 2.0;
    REAL8 PN2 = (-107 + 55*eta)/21.;
    REAL8 PN4 = (-2173 - 7483*eta + 2047*eta2)/756.;
    REAL8 ImPN5 = ((-112.)/5.)*eta;
    REAL8 lnxOx0 = log(v2/CONST_X0);
    REAL8 PN6 = (287761194. - 41727750.*eta - 36469800.*eta2 + 5731750.*eta3 + 2130975.*eta*CST_PISQ)/2494800.;
    PN6 = PN6 + 1712.*lnxOx0/105. ;
    REAL8 lnxOx0p = log(v2/CONST_X0P);
    COMPLEX16 PN3_tail = 2.*(3.*I*lnxOx0p + CST_PI);
    COMPLEX16 PN5_tail = (-107. + 34.*eta)*(3.*I*lnxOx0p + CST_PI)/21.;
    COMPLEX16 PN6_tail = -515063./22050. - 18*lnxOx0p*lnxOx0p + lnxOx0p*(-428./35. + (12.*I)*CST_PI) + ((428.*I)/105.)*CST_PI + (2.*CST_PISQ)/3.;
    return PN0 + 
        v2*(PN2 + 
            v*(PN3_tail + 
                v*(PN4 + 
                    v*(I*ImPN5 + PN5_tail + 
                        v*(PN6 + PN6_tail)))));
    // return (PN0 + 
    //     v2 * (PN2 + 
    //         v*(PN3_tail + 
    //             v*(PN4 + 
    //                 v*( I*ImPN5 + PN5_tail
    //                 )))));

}

COMPLEX16 calculate_h22_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];
    REAL8 e = dc->ePow[1], e2 = dc->ePow[2];
    COMPLEX16 eilvec[50], eilnvec[50];
    eilvec[0] = 1.0;
    eilnvec[0] = 1.0;
    for (size_t p=1; p<50; p++) {
        eilvec[p] = cexp(I*((REAL8)p)*l);
        eilnvec[p] = cexp(-I*((REAL8)p)*l);
    }
#define EIL(x) eilvec[x]
#define EILN(x) eilnvec[x]
    COMPLEX16 PN0 = (2*EILN(2) + e*((-3*EILN(1))/2. + (9*EILN(3))/2. + e*(-5*EILN(2) + 8*EILN(4) + e*((-7*EIL(1))/48. + (13*EILN(1))/16. - (171*EILN(3))/16. + (625*EILN(5))/48. + e*(-0.125*EIL(2) + (23*EILN(2))/8. - 20*EILN(4) + (81*EILN(6))/4. + e*((-47*EIL(1))/768. - (153*EIL(3))/1280. + (5*EILN(1))/384. + (963*EILN(3))/128. - (26875*EILN(5))/768. + (117649*EILN(7))/3840. + e*((-11*EIL(2))/240. - (11*EIL(4))/90. - (65*EILN(2))/144. + (101*EILN(4))/6. - (2349*EILN(6))/40. + (2048*EILN(8))/45. + e*((-1091*EIL(1))/30720. - (63*EIL(3))/2048. - (1875*EIL(5))/14336. + (227*EILN(1))/6144. - (4311*EILN(3))/2048. + (210625*EILN(5))/6144. - (588245*EILN(7))/6144. + (4782969*EILN(9))/71680. + e*((-173*EIL(2))/5760. - (17*EIL(4))/1260. - (81*EIL(6))/560. + (85*EILN(2))/1152. - (1177*EILN(4))/180. + (20979*EILN(6))/320. - (48128*EILN(8))/315. + (390625*EILN(10))/4032. + e*((-10315*EIL(1))/442368. - (15507*EIL(3))/573440. + (33125*EIL(5))/4.128768e6 - (4353013*EIL(7))/2.654208e7 + (34349*EILN(1))/1.47456e6 + (58689*EILN(3))/163840. - (7439375*EILN(5))/442368. + (88354399*EILN(7))/737280. - (110008287*EILN(9))/458752. + (25937424601*EILN(11))/1.8579456e8 + e*((-177*EIL(2))/8960. - (53*EIL(4))/2016. + (81*EIL(6))/2240. - (128*EIL(8))/675. + (1007*EILN(2))/57600. + (2159*EILN(4))/1440. - (4941*EILN(6))/128. + (66688*EILN(8))/315. - (8984375*EILN(10))/24192. + (34992*EILN(12))/175.)))))))))));
    COMPLEX16 PN2 = eta*((55*EILN(2))/21. + e*((-9*EILN(1))/28. + (21*EILN(3))/4. + e*((-23*EILN(2))/42. + (172*EILN(4))/21. + e*((965*EIL(1))/2016. - (575*EILN(1))/672. - (165*EILN(3))/224. + (23125*EILN(5))/2016. + e*((193*EIL(2))/336. - (727*EILN(2))/336. - (2*EILN(4))/7. + (837*EILN(6))/56. + e*((8809*EIL(1))/32256. + (2553*EIL(3))/3584. - (3691*EILN(1))/16128. - (9255*EILN(3))/1792. + (63125*EILN(5))/32256. + (84035*EILN(7))/4608. + e*((313*EIL(2))/1440. + (3391*EIL(4))/3780. + (235*EILN(2))/864. - (2929*EILN(4))/252. + (4563*EILN(6))/560. + (19456*EILN(8))/945. + e*((288577*EIL(1))/1.29024e6 + (18507*EIL(3))/143360. + (2061875*EIL(5))/1.806336e6 - (58169*EILN(1))/258048. + (60159*EILN(3))/28672. - (6441875*EILN(5))/258048. + (453789*EILN(7))/20480. + (20726199*EILN(9))/1.00352e6 + e*((9665*EIL(2))/48384. - (37*EIL(4))/1512. + (2295*EIL(6))/1568. - (12461*EILN(2))/48384. + (10981*EILN(4))/1512. - (46197*EILN(6))/896. + (67072*EILN(8))/1323. + (390625*EILN(10))/24192. + e*((17788873*EIL(1))/9.289728e7 + (1605111*EIL(3))/8.02816e6 - (48446875*EIL(5))/1.73408256e8 + (300694037*EIL(7))/1.5925248e8 - (11857051*EILN(1))/6.193152e7 - (1525269*EILN(3))/2.29376e6 + (53059375*EILN(5))/2.654208e6 - (455284823*EILN(7))/4.42368e6 + (3375181791*EILN(9))/3.211264e7 + (25937424601*EILN(11))/7.80337152e9 + e*((111577*EIL(2))/677376. + (96653*EIL(4))/423360. - (21519*EIL(6))/31360. + (16192*EIL(8))/6615. - (11129*EILN(2))/69120. - (146123*EILN(4))/60480. + (3060369*EILN(6))/62720. - (263744*EILN(8))/1323. + (23046875*EILN(10))/112896. - (5832*EILN(12))/245.)))))))))))
        +((-107*EILN(2))/21. + e*((-347*EILN(1))/28. + (27*EILN(3))/4. + e*((-615*EILN(2))/14. + (722*EILN(4))/21. + e*((-6085*EIL(1))/2016. + (2987*EILN(1))/672. - (21771*EILN(3))/224. + (169975*EILN(5))/2016. + e*((-257*EIL(2))/84. + (4895*EILN(2))/336. - (1362*EILN(4))/7. + (18783*EILN(6))/112. + e*((-26227*EIL(1))/10752. - (63657*EIL(3))/17920. + (23047*EILN(1))/16128. + (96255*EILN(3))/1792. - (3955675*EILN(5))/10752. + (6946289*EILN(7))/23040. + e*((-1241*EIL(2))/720. - (6509*EIL(4))/1512. - (4283*EILN(2))/864. + (73691*EILN(4))/504. - (149427*EILN(6))/224. + (965207*EILN(8))/1890. + e*((-579133*EIL(1))/258048. - (180567*EIL(3))/143360. - (3212225*EIL(5))/602112. + (403517*EILN(1))/258048. - (506511*EILN(3))/28672. + (87800575*EILN(5))/258048. - (215713141*EILN(7))/184320. + (166973481*EILN(9))/200704. + e*((-46997*EIL(2))/30240. - (1511*EIL(4))/2520. - (841329*EIL(6))/125440. - (2977*EILN(2))/1344. - (426299*EILN(4))/7560. + (6460983*EILN(6))/8960. - (4407986*EILN(8))/2205. + (127586975*EILN(10))/96768. + e*((-192930401*EIL(1))/9.289728e7 - (11485683*EIL(3))/8.02816e6 + (77812075*EIL(5))/1.73408256e8 - (1356615841*EIL(7))/1.5925248e8 + (19754003*EILN(1))/1.2386304e7 - (3614283*EILN(3))/2.29376e6 - (422579675*EILN(5))/2.654208e6 + (6347896079*EILN(7))/4.42368e6 - (107271707427*EILN(9))/3.211264e7 + (3195309060779*EILN(11))/1.560674304e9 + e*((-2193913*EIL(2))/1.69344e6 - (252085*EIL(4))/169344. + (527847*EIL(6))/250880. - (69242507*EIL(8))/6.3504e6 - (389417*EILN(2))/172800. + (29763*EILN(4))/4480. - (10098747*EILN(6))/25088. + (1154748377*EILN(8))/423360. - (22273710925*EILN(10))/4.064256e6 + (245372523*EILN(12))/78400.)))))))))));

    COMPLEX16 PN4 = eta2*((2047*EILN(2))/756. + e*((-167*EILN(1))/432. + (63*EILN(3))/16. + e*((-11*EILN(2))/1512. + (25*EILN(4))/27. + e*((-45181*EIL(1))/72576. + (7631*EILN(1))/24192. + (6357*EILN(3))/896. - (873125*EILN(5))/72576. + e*((-37117*EIL(2))/36288. - (29207*EILN(2))/12096. + (38495*EILN(4))/1134. - (10179*EILN(6))/224. + e*((-92627*EIL(1))/165888. - (117753*EIL(3))/71680. + (54497*EILN(1))/82944. - (116877*EILN(3))/7168. + (17830625*EILN(5))/165888. - (97699091*EILN(7))/829440. + e*((-162853*EIL(2))/362880. - (347813*EIL(4))/136080. + (466033*EILN(2))/217728. - (576533*EILN(4))/9072. + (627399*EILN(6))/2240. - (2216192*EILN(8))/8505. + e*((-6187901*EIL(1))/9.289728e6 - (98379*EIL(3))/573440. - (757275625*EIL(5))/1.95084288e8 + (18374435*EILN(1))/2.7869184e7 + (181935*EILN(3))/16384. - (1812143125*EILN(5))/9.289728e6 + (12921608161*EILN(7))/1.990656e7 - (422495595*EILN(9))/802816. + e*((-5295443*EIL(2))/8.70912e6 + (129943*EIL(4))/272160. - (25947*EIL(6))/4480. + (552619*EILN(2))/1.741824e6 + (12676649*EILN(4))/272160. - (9297801*EILN(6))/17920. + (11818112*EILN(8))/8505. - (871484375*EILN(10))/870912. + e*((-504800089*EIL(1))/6.68860416e8 - (3027093*EIL(3))/4.58752e6 + (11137384375*EIL(5))/6.242697216e9 - (48889932721*EIL(7))/5.73308928e9 + (1683111103*EILN(1))/2.22953472e9 - (20999247*EILN(3))/9.17504e6 + (15210128125*EILN(5))/9.5551488e7 - (199326633373*EILN(7))/1.5925248e8 + (72131955489*EILN(9))/2.5690112e7 - (73117599950219*EILN(11))/4.013162496e10 + e*((-235113733*EIL(2))/3.6578304e8 - (13402727*EIL(4))/1.524096e7 + (525501*EIL(6))/125440. - (33298352*EIL(8))/2.679075e6 + (58157569*EILN(2))/8.70912e7 - (104909341*EILN(4))/6.53184e6 + (117888021*EILN(6))/250880. - (167347088*EILN(8))/59535. + (595745703125*EILN(10))/1.09734912e8 - (3932226*EILN(12))/1225.)))))))))))
        +eta*((-1069*EILN(2))/108. + e*((68801*EILN(1))/3024. - (4401*EILN(3))/112. + e*((119087*EILN(2))/1512. - (17704*EILN(4))/189. + e*((772885*EIL(1))/72576. + (820657*EILN(1))/24192. + (118227*EILN(3))/896. - (13685275*EILN(5))/72576. + e*((61693*EIL(2))/5184. + (1328891*EILN(2))/12096. + (266311*EILN(4))/1134. - (39261*EILN(6))/112. + e*((18029933*EIL(1))/1.161216e6 + (235653*EIL(3))/14336. + (2444575*EILN(1))/82944. + (130659*EILN(3))/1024. + (524032825*EILN(5))/1.161216e6 - (103354769*EILN(7))/165888. + e*((3509791*EIL(2))/362880. + (229381*EIL(4))/9720. + (27078971*EILN(2))/217728. + (353665*EILN(4))/4536. + (142887*EILN(6))/160. - (36603169*EILN(8))/34020. + e*((202057829*EIL(1))/9.289728e6 + (3672603*EIL(3))/573440. + (6648483505*EIL(5))/1.95084288e8 + (680994445*EILN(1))/2.7869184e7 + (20274927*EILN(3))/114688. - (1197309755*EILN(5))/9.289728e6 + (35008940183*EILN(7))/1.990656e7 - (1461536325*EILN(9))/802816. + e*((22278913*EIL(2))/1.741824e6 + (2110013*EIL(4))/1.90512e6 + (12321021*EIL(6))/250880. + (229407323*EILN(2))/1.741824e6 + (71620813*EILN(4))/272160. - (444201*EILN(6))/640. + (161681837*EILN(8))/47628. - (36974620775*EILN(10))/1.2192768e7 + e*((18439985905*EIL(1))/6.68860416e8 + (72350151*EIL(3))/6.422528e6 - (55217414815*EIL(5))/6.242697216e9 + (404077391809*EIL(7))/5.73308928e9 + (42623026193*EILN(1))/2.22953472e9 + (1524563703*EILN(3))/9.17504e6 + (323688312965*EILN(5))/6.68860416e8 - (64413584431*EILN(7))/3.1850496e7 + (164287746639*EILN(9))/2.5690112e7 - (1401019653133243*EILN(11))/2.8092137472e11 + e*((5227772287*EIL(2))/3.6578304e8 + (46500989*EIL(4))/3.81024e6 - (2695747*EIL(6))/100352. + (34536950011*EIL(8))/3.429216e8 + (12438635141*EILN(2))/8.70912e7 + (105559135*EILN(4))/653184. + (16690287*EILN(6))/15680. - (37268977379*EILN(8))/7.62048e6 + (2581176072725*EILN(10))/2.19469824e8 - (1272131627*EILN(12))/156800.)))))))))))
        +((-2173*EILN(2))/756. + e*((-12172*EILN(1))/189. + 36*EILN(3) + e*((-493351*EILN(2))/1512. + (34997*EILN(4))/189. + e*((128267*EIL(1))/9072. + (138311*EILN(1))/3024. - (83637*EILN(3))/112. + (653215*EILN(5))/1296. + e*((788509*EIL(2))/36288. - (1075003*EILN(2))/12096. - (810100*EILN(4))/567. + (121383*EILN(6))/112. + e*((1145587*EIL(1))/36288. + (115911*EIL(3))/4480. + (9706883*EILN(1))/72576. - (125829*EILN(3))/896. - (93734485*EILN(5))/36288. + (106906457*EILN(7))/51840. + e*((15823573*EIL(2))/362880. + (300283*EIL(4))/9720. - (4288603*EILN(2))/217728. + (117037*EILN(4))/4536. - (2543769*EILN(6))/560. + (61980769*EILN(8))/17010. + e*((313766687*EIL(1))/5.80608e6 + (3230271*EIL(3))/71680. + (912511655*EIL(5))/2.4385536e7 + (740835881*EILN(1))/3.483648e6 - (329265*EILN(3))/2048. + (721209905*EILN(5))/1.161216e6 - (19485763199*EILN(7))/2.48832e6 + (3074814891*EILN(9))/501760. + e*((627362447*EIL(2))/8.70912e6 + (6246337*EIL(4))/136080. + (1439127*EIL(6))/31360. + (75022415*EILN(2))/1.741824e6 - (40906771*EILN(4))/136080. + (18360261*EILN(6))/8960. - (6319072693*EILN(8))/476280. + (1238053355*EILN(10))/124416. + e*((33421263553*EIL(1))/4.1803776e8 + (57610881*EIL(3))/802816. + (4526581835*EIL(5))/9.7542144e7 + (10180675769*EIL(7))/1.7915904e8 + (40328664119*EILN(1))/1.3934592e8 - (62398179*EILN(3))/573440. - (45803229275*EILN(5))/8.3607552e7 + (20050668187*EILN(7))/3.981312e6 - (44397572337*EILN(9))/2.00704e6 + (138199671373093*EILN(11))/8.77879296e9 + e*((7675116689*EIL(2))/7.3156608e7 + (309712789*EIL(4))/4.35456e6 + (11546379*EIL(6))/250880. + (24311017411*EIL(8))/3.429216e8 + (9205006841*EILN(2))/8.70912e7 - (1424149817*EILN(4))/6.53184e6 - (142764651*EILN(6))/125440. + (16525678577*EILN(8))/1.524096e6 - (142410234185*EILN(10))/3.919104e6 + (15309258417*EILN(12))/627200.)))))))))));
    COMPLEX16 ImPN5 = eta*((-112*EILN(2))/5. + e*((1591*EILN(1))/210. - (22303*EILN(3))/210. + e*((7934*EILN(2))/105. - (2276*EILN(4))/7. + e*((185*EIL(1))/336. - (6523*EILN(1))/1680. + (547871*EILN(3))/1680. - (272749*EILN(5))/336. + e*((569*EIL(2))/630. - (799*EILN(2))/14. + (321137*EILN(4))/315. - (188618*EILN(6))/105. + e*((100049*EIL(1))/80640. + (30979*EIL(3))/19200. + (31033*EILN(1))/13440. - (4000889*EILN(3))/13440. + (43100411*EILN(5))/16128. - (70441019*EILN(7))/19200. + e*((3649*EIL(2))/3150. + (11587*EIL(4))/4200. + (8489*EILN(2))/504. - (2721431*EILN(4))/2520. + (26297843*EILN(6))/4200. - (22286939*EILN(8))/3150. + e*((436249*EIL(1))/215040. + (1039957*EIL(3))/1.0752e6 + (1742809*EIL(5))/387072. + (4562717*EILN(1))/1.93536e6 + (24771743*EILN(3))/215040. - (46017315*EILN(5))/14336. + (18753180241*EILN(7))/1.3824e6 - (187350839*EILN(9))/14336. + e*((11471*EIL(2))/5600. + (116749*EIL(4))/264600. + (166847*EIL(6))/23520. + (28769*EILN(2))/30240. + (1345231*EILN(4))/2520. - (20137571*EILN(6))/2400. + (7339008049*EILN(8))/264600. - (164692189*EILN(10))/7056. + e*((71465701*EIL(1))/2.58048e7 + (8514333*EIL(3))/4.01408e6 - (66957463*EIL(5))/8.6704128e7 + (1445484221*EIL(7))/1.327104e8 + (472103167*EILN(1))/1.548288e8 - (15581443*EILN(3))/819200. + (5961007273*EILN(5))/3.096576e6 - (44313428831*EILN(7))/2.21184e6 + (4353520139413*EILN(9))/8.02816e7 - (263872722844709*EILN(11))/6.5028096e9 + e*((4510817*EIL(2))/1.5876e6 + (1647077*EIL(4))/705600. - (379801*EIL(6))/117600. + (3125792389*EIL(8))/1.90512e8 + (579563*EILN(2))/168000. - (16103491*EILN(4))/113400. + (923356709*EILN(6))/156800. - (63096721069*EILN(8))/1.4112e6 + (155929986337*EILN(10))/1.524096e6 - (81148493863*EILN(12))/1.176e6)))))))))));

    COMPLEX16 PN6 = eta3*((114635*EILN(2))/49896. + e*((-1685*EILN(1))/7392. + (44837*EILN(3))/22176. + e*((4853*EILN(2))/99792. - (48721*EILN(4))/12474. + e*((1873057*EIL(1))/4.790016e6 - (6421*EILN(1))/228096. + (321877*EILN(3))/177408. - (90319375*EILN(5))/4.790016e6 + e*((658349*EIL(2))/798336. + (2206069*EILN(2))/798336. + (21155*EILN(4))/8316. - (221005*EILN(6))/4928. + e*((4807399*EIL(1))/6.967296e6 + (4775839*EIL(3))/2.838528e6 - (33236759*EILN(1))/3.8320128e7 + (23898079*EILN(3))/1.419264e6 - (885284375*EILN(5))/7.6640256e7 - (862787345*EILN(7))/1.0948608e7 + e*((13093763*EIL(2))/2.395008e7 + (5796055*EIL(4))/1.796256e6 - (52359799*EILN(2))/1.4370048e7 + (41497411*EILN(4))/598752. - (74509*EILN(6))/896. - (29439104*EILN(8))/280665. + e*((3796518941*EIL(1))/3.06561024e9 + (8386901*EIL(3))/1.1354112e8 + (25185274375*EIL(5))/4.291854336e9 - (68148095*EILN(1))/5.5738368e7 - (469812895*EILN(3))/2.2708224e7 + (143666545625*EILN(5))/6.13122048e8 - (46756586597*EILN(7))/1.4598144e8 - (7052281119*EILN(9))/8.830976e7 + e*((625277249*EIL(2))/5.7480192e8 - (158269747*EIL(4))/1.2573792e8 + (7062661*EIL(6))/689920. - (6185071*EILN(2))/1.6422912e7 - (1643681483*EILN(4))/1.796256e7 + (274305809*EILN(6))/394240. - (1900991296*EILN(8))/1.964655e6 + (37710546875*EILN(10))/4.02361344e8 + e*((7595600363*EIL(1))/4.013162496e9 + (7466017441*EIL(3))/6.35830272e9 - (163345253125*EIL(5))/3.7456183296e10 + (6541995386521*EIL(7))/3.7838389248e11 - (723606359*EILN(1))/3.82205952e8 + (11583484141*EILN(3))/1.81665792e9 - (14613552371875*EILN(5))/4.4144787456e10 + (19734762643421*EILN(7))/1.051066368e10 - (1435834243971*EILN(9))/5.65182464e8 + (207237664614299*EILN(11))/3.37105649664e11 + e*((12340976569*EIL(2))/8.04722688e9 + (155755547*EIL(4))/9.144576e7 - (4263151*EIL(6))/394240. + (55922632*EIL(8))/1.964655e6 - (1849354723*EILN(2))/1.14960384e9 + (6286453937*EILN(4))/1.4370048e8 - (5741163097*EILN(6))/5.51936e6 + (1319992328*EILN(8))/280665. - (1628212890625*EILN(10))/2.68240896e8 + (5026401*EILN(12))/2695.)))))))))))
        +eta2*((-20261*EILN(2))/1386. + e*((-319607*EILN(1))/16632. - (3291*EILN(3))/154. + e*((-55691*EILN(2))/924. - (380645*EILN(4))/8316. + e*((-356087*EIL(1))/22176. - (2356135*EILN(1))/14784. + (229821*EILN(3))/2464. - (905225*EILN(5))/6336. + e*((-10821985*EIL(2))/798336. - (513127*EILN(2))/924. + (3372377*EILN(4))/6237. - (598393*EILN(6))/1408. + e*((-2383333*EIL(1))/48384. - (1050689*EIL(3))/56320. - (915565001*EILN(1))/3.193344e6 - (17266255*EILN(3))/19712. + (161636675*EILN(5))/101376. - (315283591*EILN(7))/285120. + e*((-6155309*EIL(2))/295680. - (60643133*EIL(4))/1.99584e6 - (49179419*EILN(2))/49896. - (67614643*EILN(4))/44352. + (34838477*EILN(6))/8960. - (1283130287*EILN(8))/498960. + e*((-1532216617*EIL(1))/1.419264e7 - (17470897*EIL(3))/3.15392e6 - (565203305*EIL(5))/1.0948608e7 - (5507823173*EILN(1))/1.3934592e7 - (398084415*EILN(3))/315392. - (7509447235*EILN(5))/2.433024e6 + (472177361369*EILN(7))/5.474304e7 - (17333576739*EILN(9))/3.15392e6 + e*((-8937400771*EIL(2))/1.9160064e8 + (53521067*EIL(4))/3.49272e6 - (138025231*EIL(6))/1.57696e6 - (2719883069*EILN(2))/1.824768e6 - (1930207337*EILN(4))/1.49688e6 - (5328510043*EILN(6))/788480. + (1497492869*EILN(8))/83160. - (5912078855275*EILN(10))/5.36481792e8 + e*((-52558675997*EIL(1))/2.7869184e8 - (1043919235*EIL(3))/3.5323904e7 + (4577062535*EIL(5))/8.6704128e7 - (1536139571029*EIL(7))/1.051066368e10 - (2217653603*EILN(1))/4.58752e6 - (3097479907*EILN(3))/1.57696e6 - (893767789295*EILN(5))/1.226244096e9 - (145528016459*EILN(7))/9.732096e6 + (25398831144339*EILN(9))/7.0647808e8 - (246236411130193*EILN(11))/1.170505728e10 + e*((-605553479791*EIL(2))/8.04722688e9 - (34777181*EIL(4))/1.354752e6 + (1157660171*EIL(6))/9.46176e6 - (3617311562453*EIL(8))/1.50885504e10 - (1299142242323*EILN(2))/6.386688e8 - (672620874011*EILN(4))/2.8740096e8 + (3121709335*EILN(6))/2.207744e6 - (516980278471*EILN(8))/1.596672e7 + (668903191944725*EILN(10))/9.656672256e9 - (799500197047*EILN(12))/2.06976e7)))))))))))
        +eta*((-16.72588985088985 + (41*CST_PISQ)/48.)*EILN(2) + e*((761.7193722943723 - (287*CST_PISQ)/16.)*EILN(1) + (-905.6077922077922 + (615*CST_PISQ)/32.)*EILN(3) + e*((2955.6335076960077 - (4879*CST_PISQ)/96.)*EILN(2) + (-2684.846103896104 + (533*CST_PISQ)/12.)*EILN(4) + e*((-68.3341019320186 - (26527*CST_PISQ)/9216.)*EIL(1) + (1450.5800850769601 - (118367*CST_PISQ)/3072.)*EILN(1) + (4064.7379058441556 - (41205*CST_PISQ)/1024.)*EILN(3) + (-5470.049665804874 + (607825*CST_PISQ)/9216.)*EILN(5) + e*((-196.02890061327562 - (8405*CST_PISQ)/3072.)*EIL(2) + (6958.09359743266 - (24149*CST_PISQ)/192.)*EILN(2) + (4229.110822510823 + (3895*CST_PISQ)/192.)*EILN(4) + (-9356.110409902598 + (64575*CST_PISQ)/1024.)*EILN(6) + e*((-219.19612130731923 - (318119*CST_PISQ)/36864.)*EIL(1) + (-308.8556742086039 - (62853*CST_PISQ)/20480.)*EIL(3) + (1073.9139304127586 - (2032985*CST_PISQ)/36864.)*EILN(1) + (11850.321022727272 - (773547*CST_PISQ)/4096.)*EILN(3) + (2901.2547403286335 + (3430675*CST_PISQ)/18432.)*EILN(5) + (-14298.237323758418 - (443989*CST_PISQ)/92160.)*EILN(7) + e*((-498.3178004248838 - (689497*CST_PISQ)/92160.)*EIL(2) + (-458.9947122347817 - (7831*CST_PISQ)/2160.)*EIL(4) + (8964.919418933048 - (1247753*CST_PISQ)/6912.)*EILN(2) + (18045.81233114879 - (179375*CST_PISQ)/576.)*EILN(4) + (-1392.25734375 + (5904369*CST_PISQ)/10240.)*EILN(6) + (-19934.701756587034 - (1837087*CST_PISQ)/8640.)*EILN(8) + e*((-506.10852300754317 - (20201971*CST_PISQ)/1.179648e6)*EIL(1) + (-626.6239077085024 - (5093061*CST_PISQ)/655360.)*EIL(3) + (-668.4674264266549 - (12102175*CST_PISQ)/2.752512e6)*EIL(5) + (-38.21274736999835 - (83518681*CST_PISQ)/1.179648e6)*EILN(1) + (14391.475189605317 - (29071665*CST_PISQ)/131072.)*EILN(3) + (28010.20895139625 - (716715875*CST_PISQ)/1.179648e6)*EILN(5) + (-11685.130675406408 + (8398445699*CST_PISQ)/5.89824e6)*EILN(7) + (-25296.45318822291 - (635815413*CST_PISQ)/917504.)*EILN(9) + e*((-1015.4603047902136 - (6184235*CST_PISQ)/442368.)*EIL(2) + (-740.2109268898356 - (8176589*CST_PISQ)/967680.)*EIL(4) + (-962.5111929462836 - (1243161*CST_PISQ)/229376.)*EIL(6) + (10736.749242090214 - (110326777*CST_PISQ)/442368.)*EILN(2) + (18262.182968374218 - (5790389*CST_PISQ)/27648.)*EILN(4) + (46020.782143618104 - (107464239*CST_PISQ)/81920.)*EILN(6) + (-33760.01003013888 + (152879939*CST_PISQ)/48384.)*EILN(8) + (-28332.744098992647 - (10362281575*CST_PISQ)/6.193152e6)*EILN(10) + e*((-967.5726015434586 - (5992167589*CST_PISQ)/2.1233664e8)*EIL(1) + (-1173.8553677996633 - (256666191*CST_PISQ)/1.835008e7)*EIL(3) + (-837.5218497427122 - (1851511825*CST_PISQ)/1.98180864e8)*EIL(5) + (-1374.95403405499 - (4303634741*CST_PISQ)/6.3700992e8)*EIL(7) + (-1868.325468210874 - (5993201773*CST_PISQ)/7.077888e7)*EILN(1) + (18012.22671196479 - (422136369*CST_PISQ)/1.31072e6)*EILN(3) + (19845.067934911116 - (3057311575*CST_PISQ)/4.2467328e7)*EILN(5) + (80084.4501431625 - (205508659489*CST_PISQ)/7.077888e7)*EILN(7) + (-78182.40967386022 + (60205143699*CST_PISQ)/9.17504e6)*EILN(9) + (-25148.48001789796 - (31510983563807*CST_PISQ)/8.91813888e9)*EILN(11) + e*((-1793.7705548433598 - (75600023*CST_PISQ)/3.44064e6)*EIL(2) + (-1312.996479204722 - (28432639*CST_PISQ)/1.93536e6)*EIL(4) + (-892.791979423339 - (3402549*CST_PISQ)/327680.)*EIL(6) + (-1952.2434823176254 - (246494419*CST_PISQ)/2.90304e7)*EIL(8) + (12078.931039469388 - (7102524587*CST_PISQ)/2.21184e7)*EILN(2) + (23786.60290893183 - (52629773*CST_PISQ)/138240.)*EILN(4) + (16152.068562242725 + (456804657*CST_PISQ)/1.14688e6)*EILN(6) + (145164.52742040646 - (1749769669*CST_PISQ)/276480.)*EILN(8) + (-163680.48478499433 + (482957643725*CST_PISQ)/3.7158912e7)*EILN(10) + (-8796.993381913846 - (2471923989*CST_PISQ)/358400.)*EILN(12))))))))))))
        +((761273*EILN(2))/6600. + e*((-152382067*EILN(1))/415800. + (40398301*EILN(3))/46200. + e*((-147962753*EILN(2))/55440. + (677312767*EILN(4))/207900. + e*((159694639*EIL(1))/554400. + (1187699041*EILN(1))/1.6632e6 - (446587649*EILN(3))/61600. + (409798571*EILN(5))/47520. + e*((605752727*EIL(2))/1.4256e6 - (1639986547*EILN(2))/739200. - (1916431603*EILN(4))/124740. + (14114214613*EILN(6))/739200. + e*((5134349299*EIL(1))/4.8384e6 + (5353429111*EIL(3))/9.856e6 + (8133496487*EILN(1))/2.28096e6 - (6385205531*EILN(3))/985600. - (64567910591*EILN(5))/2.128896e6 + (4336294770377*EILN(7))/1.14048e8 + e*((67756091159*EIL(2))/4.9896e7 + (11824992307*EIL(4))/1.6632e7 + (4155388231*EILN(2))/7.98336e6 - (99612938233*EILN(4))/9.9792e6 - (13059282769*EILN(6))/224000. + (1762744976863*EILN(8))/2.4948e7 + e*((422365569391*EIL(1))/1.77408e8 + (5903295829*EIL(3))/3.9424e6 + (181713347579*EIL(5))/1.9160064e8 + (671224938311*EILN(1))/8.70912e7 - (52412498951*EILN(3))/1.18272e7 - (209370926179*EILN(5))/2.128896e7 - (3013625336251*EILN(7))/2.737152e7 + (2466676751741*EILN(9))/1.9712e7 + e*((6882216100099*EIL(2))/2.395008e9 + (582294694181*EIL(4))/3.49272e8 + (50483941351*EIL(6))/3.9424e7 + (17959715311*EILN(2))/3.99168e6 - (1473852135317*EILN(4))/1.49688e8 + (925457277*EILN(6))/2.816e6 - (5097022799249*EILN(8))/2.4948e7 + (573926881083763*EILN(10))/2.68240896e9 + e*((2394746934463*EIL(1))/5.5738368e8 + (13196697566809*EIL(3))/4.415488e9 + (5794647883247*EIL(5))/3.121348608e9 + (1068369967161823*EIL(7))/6.13122048e11 + (364385232709751*EILN(1))/2.7869184e10 - (1528957700793*EILN(3))/1.261568e9 - (20498589752735*EILN(5))/1.226244096e9 + (5192997978681389*EILN(7))/1.53280512e11 - (1314755932093689*EILN(9))/3.5323904e9 + (138974438431062131*EILN(11))/3.90168576e11 + e*((36169149491593*EIL(2))/7.185024e9 + (15310024523*EIL(4))/4.8384e6 + (1123977163709*EIL(6))/5.51936e8 + (1797970898569537*EIL(8))/7.5442752e11 + (77919078658453*EILN(2))/7.98336e9 - (49331634352141*EILN(4))/7.185024e9 - (4753304951261*EILN(6))/1.655808e8 + (655707869783921*EILN(8))/5.588352e9 - (918292495952987*EILN(10))/1.379524608e9 + (42926349119957*EILN(12))/7.392e7)))))))))));

    REAL8 lgxOx0p = log(v2/CONST_X0P);
    COMPLEX16 PN3_tail = (CST_2PI*(EILN(2) + e*((-3*EILN(1))/8. + (27*EILN(3))/8. + e*((-5*EILN(2))/2. + 8*EILN(4) + e*((-7*EIL(1))/192. + (13*EILN(1))/64. - (513*EILN(3))/64. + (3125*EILN(5))/192. + e*(-0.0625*EIL(2) + (23*EILN(2))/16. - 20*EILN(4) + (243*EILN(6))/8. + e*((-47*EIL(1))/3072. - (459*EIL(3))/5120. + (5*EILN(1))/1536. + (2889*EILN(3))/512. - (134375*EILN(5))/3072. + (823543*EILN(7))/15360. + e*((-11*EIL(2))/480. - (11*EIL(4))/90. - (65*EILN(2))/288. + (101*EILN(4))/6. - (7047*EILN(6))/80. + (4096*EILN(8))/45. + e*((-1091*EIL(1))/122880. - (189*EIL(3))/8192. - (9375*EIL(5))/57344. + (227*EILN(1))/24576. - (12933*EILN(3))/8192. + (1053125*EILN(5))/24576. - (4117715*EILN(7))/24576. + (43046721*EILN(9))/286720. + e*((-173*EIL(2))/11520. - (17*EIL(4))/1260. - (243*EIL(6))/1120. + (85*EILN(2))/2304. - (1177*EILN(4))/180. + (62937*EILN(6))/640. - (96256*EILN(8))/315. + (1953125*EILN(10))/8064. + e*((-10315*EIL(1))/1.769472e6 - (46521*EIL(3))/2.29376e6 + (165625*EIL(5))/1.6515072e7 - (30471091*EIL(7))/1.0616832e8 + (34349*EILN(1))/5.89824e6 + (176067*EILN(3))/655360. - (37196875*EILN(5))/1.769472e6 + (618480793*EILN(7))/2.94912e6 - (990074583*EILN(9))/1.835008e6 + (285311670611*EILN(11))/7.4317824e8 + e*((-177*EIL(2))/17920. - (53*EIL(4))/2016. + (243*EIL(6))/4480. - (256*EIL(8))/675. + (1007*EILN(2))/115200. + (2159*EILN(4))/1440. - (14823*EILN(6))/256. + (133376*EILN(8))/315. - (44921875*EILN(10))/48384. + (104976*EILN(12))/175.)))))))))))
        +I*(6*lgxOx0p*EILN(2) + e*(lgxOx0p*((-9*EILN(1))/4. + (81*EILN(3))/4.) + (3*EILN(1)*CST_LN2)/2. + EILN(3)*((-27*CST_LN2)/2. + (27*CST_LN3)/2.) + e*(lgxOx0p*(-15*EILN(2) + 48*EILN(4)) + 32*EILN(4)*CST_LN2 + e*(lgxOx0p*((7*EIL(1))/32. + (39*EILN(1))/32. - (1539*EILN(3))/32. + (3125*EILN(5))/32.) - (7*EIL(1)*CST_LN2)/48. - (13*EILN(1)*CST_LN2)/16. + EILN(3)*((513*CST_LN2)/16. - (513*CST_LN3)/16.) + EILN(5)*((-3125*CST_LN2)/48. + (3125*CST_LN5)/48.) + e*(lgxOx0p*((3*EIL(2))/8. + (69*EILN(2))/8. - 120*EILN(4) + (729*EILN(6))/4.) - 80*EILN(4)*CST_LN2 + (243*EILN(6)*CST_LN3)/2. + e*(lgxOx0p*((47*EIL(1))/512. + (1377*EIL(3))/2560. + (5*EILN(1))/256. + (8667*EILN(3))/256. - (134375*EILN(5))/512. + (823543*EILN(7))/2560.) - (47*EIL(1)*CST_LN2)/768. - (5*EILN(1)*CST_LN2)/384. + EIL(3)*((-459*CST_LN2)/1280. + (459*CST_LN3)/1280.) + EILN(3)*((-2889*CST_LN2)/128. + (2889*CST_LN3)/128.) + EILN(5)*((134375*CST_LN2)/768. - (134375*CST_LN5)/768.) + EILN(7)*((-823543*CST_LN2)/3840. + (823543*CST_LN7)/3840.) + e*(lgxOx0p*((11*EIL(2))/80. + (11*EIL(4))/15. - (65*EILN(2))/48. + 101*EILN(4) - (21141*EILN(6))/40. + (8192*EILN(8))/15.) + (22*EIL(4)*CST_LN2)/45. + (202*EILN(4)*CST_LN2)/3. + (32768*EILN(8)*CST_LN2)/45. - (7047*EILN(6)*CST_LN3)/20. + e*(lgxOx0p*((1091*EIL(1))/20480. + (567*EIL(3))/4096. + (28125*EIL(5))/28672. + (227*EILN(1))/4096. - (38799*EILN(3))/4096. + (1053125*EILN(5))/4096. - (4117715*EILN(7))/4096. + (129140163*EILN(9))/143360.) - (1091*EIL(1)*CST_LN2)/30720. - (227*EILN(1)*CST_LN2)/6144. + EILN(3)*((12933*CST_LN2)/2048. - (12933*CST_LN3)/2048.) + EIL(3)*((-189*CST_LN2)/2048. + (189*CST_LN3)/2048.) + EILN(9)*((-43046721*CST_LN2)/71680. + (43046721*CST_LN3)/35840.) + EIL(5)*((-9375*CST_LN2)/14336. + (9375*CST_LN5)/14336.) + EILN(5)*((-1053125*CST_LN2)/6144. + (1053125*CST_LN5)/6144.) + EILN(7)*((4117715*CST_LN2)/6144. - (4117715*CST_LN7)/6144.) + e*(lgxOx0p*((173*EIL(2))/1920. + (17*EIL(4))/210. + (729*EIL(6))/560. + (85*EILN(2))/384. - (1177*EILN(4))/30. + (188811*EILN(6))/320. - (192512*EILN(8))/105. + (1953125*EILN(10))/1344.) + (17*EIL(4)*CST_LN2)/315. - (1177*EILN(4)*CST_LN2)/45. - (770048*EILN(8)*CST_LN2)/315. + (243*EIL(6)*CST_LN3)/280. + (62937*EILN(6)*CST_LN3)/160. + (1953125*EILN(10)*CST_LN5)/2016. + e*(lgxOx0p*((10315*EIL(1))/294912. + (139563*EIL(3))/1.14688e6 - (165625*EIL(5))/2.752512e6 + (30471091*EIL(7))/1.769472e7 + (34349*EILN(1))/983040. + (528201*EILN(3))/327680. - (37196875*EILN(5))/294912. + (618480793*EILN(7))/491520. - (2970223749*EILN(9))/917504. + (285311670611*EILN(11))/1.2386304e8) - (10315*EIL(1)*CST_LN2)/442368. - (34349*EILN(1)*CST_LN2)/1.47456e6 + EILN(9)*((990074583*CST_LN2)/458752. - (990074583*CST_LN3)/229376.) + EIL(3)*((-46521*CST_LN2)/573440. + (46521*CST_LN3)/573440.) + EILN(3)*((-176067*CST_LN2)/163840. + (176067*CST_LN3)/163840.) + EILN(5)*((37196875*CST_LN2)/442368. - (37196875*CST_LN5)/442368.) + EIL(5)*((165625*CST_LN2)/4.128768e6 - (165625*CST_LN5)/4.128768e6) + e*(lgxOx0p*((531*EIL(2))/8960. + (53*EIL(4))/336. - (729*EIL(6))/2240. + (512*EIL(8))/225. + (1007*EILN(2))/19200. + (2159*EILN(4))/240. - (44469*EILN(6))/128. + (266752*EILN(8))/105. - (44921875*EILN(10))/8064. + (629856*EILN(12))/175.) + (53*EIL(4)*CST_LN2)/504. + (2048*EIL(8)*CST_LN2)/675. + (2159*EILN(4)*CST_LN2)/360. + (1067008*EILN(8)*CST_LN2)/315. - (243*EIL(6)*CST_LN3)/1120. - (14823*EILN(6)*CST_LN3)/64. + EILN(12)*((419904*CST_LN2)/175. + (419904*CST_LN3)/175.) - (44921875*EILN(10)*CST_LN5)/12096.) + EIL(7)*((-30471091*CST_LN2)/2.654208e7 + (30471091*CST_LN7)/2.654208e7) + EILN(7)*((-618480793*CST_LN2)/737280. + (618480793*CST_LN7)/737280.) + EILN(11)*((-285311670611*CST_LN2)/1.8579456e8 + (285311670611*CST_LN11)/1.8579456e8))))))))))));

                    
    COMPLEX16 PN5_tail = eta*(CST_2PI*((17*EILN(2))/21. + e*((3*EILN(1))/28. + (9*EILN(3))/4. + e*((41*EILN(2))/42. + (88*EILN(4))/21. + e*((139*EIL(1))/1008. - (53*EILN(1))/168. + (387*EILN(3))/112. + (3125*EILN(5))/504. + e*((107*EIL(2))/336. - (605*EILN(2))/336. + (68*EILN(4))/7. + (405*EILN(6))/56. + e*((2449*EIL(1))/32256. + (10377*EIL(3))/17920. - (949*EILN(1))/16128. - (11997*EILN(3))/1792. + (784375*EILN(5))/32256. + (117649*EILN(7))/23040. + e*((173*EIL(2))/1440. + (1811*EIL(4))/1890. + (215*EILN(2))/864. - (2525*EILN(4))/126. + (31509*EILN(6))/560. - (4096*EILN(8))/945. + e*((4867*EIL(1))/80640. + (7767*EIL(3))/71680. + (340625*EIL(5))/225792. - (7867*EILN(1))/129024. + (8469*EILN(3))/3584. - (6790625*EILN(5))/129024. + (117649*EILN(7))/960. - (14348907*EILN(9))/501760. + e*((25979*EIL(2))/241920. - (67*EIL(4))/3780. + (18063*EIL(6))/7840. - (7123*EILN(2))/48384. + (39811*EILN(4))/3780. - (566757*EILN(6))/4480. + (1681408*EILN(8))/6615. - (1953125*EILN(10))/24192. + e*((4717987*EIL(1))/9.289728e7 + (257049*EIL(3))/1.605632e6 - (61428125*EIL(5))/1.73408256e8 + (549067883*EIL(7))/1.5925248e8 - (628919*EILN(1))/1.2386304e7 - (1452069*EILN(3))/2.29376e6 + (94221875*EILN(5))/2.654208e6 - (252121807*EILN(7))/884736. + (16257311631*EILN(9))/3.211264e7 - (285311670611*EILN(11))/1.560674304e9 + e*((295669*EIL(2))/3.38688e6 + (51109*EIL(4))/211680. - (33129*EIL(6))/31360. + (18688*EIL(8))/3675. - (29333*EILN(2))/345600. - (95731*EILN(4))/30240. + (6406371*EILN(6))/62720. - (4037888*EILN(8))/6615. + (330078125*EILN(10))/338688. - (454896*EILN(12))/1225.)))))))))))
        +I*((34*lgxOx0p*EILN(2))/7. + e*(lgxOx0p*((9*EILN(1))/14. + (27*EILN(3))/2.) - (3*EILN(1)*CST_LN2)/7. + EILN(3)*(-9*CST_LN2 + 9*CST_LN3) + e*(lgxOx0p*((41*EILN(2))/7. + (176*EILN(4))/7.) + (352*EILN(4)*CST_LN2)/21. + e*(lgxOx0p*((-139*EIL(1))/168. - (53*EILN(1))/28. + (1161*EILN(3))/56. + (3125*EILN(5))/84.) + (139*EIL(1)*CST_LN2)/252. + (53*EILN(1)*CST_LN2)/42. + EILN(3)*((-387*CST_LN2)/28. + (387*CST_LN3)/28.) + EILN(5)*((-3125*CST_LN2)/126. + (3125*CST_LN5)/126.) + e*(lgxOx0p*((-107*EIL(2))/56. - (605*EILN(2))/56. + (408*EILN(4))/7. + (1215*EILN(6))/28.) + (272*EILN(4)*CST_LN2)/7. + (405*EILN(6)*CST_LN3)/14. + e*(lgxOx0p*((-2449*EIL(1))/5376. - (31131*EIL(3))/8960. - (949*EILN(1))/2688. - (35991*EILN(3))/896. + (784375*EILN(5))/5376. + (117649*EILN(7))/3840.) + (2449*EIL(1)*CST_LN2)/8064. + (949*EILN(1)*CST_LN2)/4032. + EILN(3)*((11997*CST_LN2)/448. - (11997*CST_LN3)/448.) + EIL(3)*((10377*CST_LN2)/4480. - (10377*CST_LN3)/4480.) + EILN(5)*((-784375*CST_LN2)/8064. + (784375*CST_LN5)/8064.) + EILN(7)*((-117649*CST_LN2)/5760. + (117649*CST_LN7)/5760.) + e*(lgxOx0p*((-173*EIL(2))/240. - (1811*EIL(4))/315. + (215*EILN(2))/144. - (2525*EILN(4))/21. + (94527*EILN(6))/280. - (8192*EILN(8))/315.) - (3622*EIL(4)*CST_LN2)/945. - (5050*EILN(4)*CST_LN2)/63. - (32768*EILN(8)*CST_LN2)/945. + (31509*EILN(6)*CST_LN3)/140. + e*(lgxOx0p*((-4867*EIL(1))/13440. - (23301*EIL(3))/35840. - (340625*EIL(5))/37632. - (7867*EILN(1))/21504. + (25407*EILN(3))/1792. - (6790625*EILN(5))/21504. + (117649*EILN(7))/160. - (43046721*EILN(9))/250880.) + (4867*EIL(1)*CST_LN2)/20160. + (7867*EILN(1)*CST_LN2)/32256. + EILN(9)*((14348907*CST_LN2)/125440. - (14348907*CST_LN3)/62720.) + EIL(3)*((7767*CST_LN2)/17920. - (7767*CST_LN3)/17920.) + EILN(3)*((-8469*CST_LN2)/896. + (8469*CST_LN3)/896.) + EILN(5)*((6790625*CST_LN2)/32256. - (6790625*CST_LN5)/32256.) + EIL(5)*((340625*CST_LN2)/56448. - (340625*CST_LN5)/56448.) + EILN(7)*((-117649*CST_LN2)/240. + (117649*CST_LN7)/240.) + e*(lgxOx0p*((-25979*EIL(2))/40320. + (67*EIL(4))/630. - (54189*EIL(6))/3920. - (7123*EILN(2))/8064. + (39811*EILN(4))/630. - (1700271*EILN(6))/2240. + (3362816*EILN(8))/2205. - (1953125*EILN(10))/4032.) + (67*EIL(4)*CST_LN2)/945. + (39811*EILN(4)*CST_LN2)/945. + (13451264*EILN(8)*CST_LN2)/6615. - (18063*EIL(6)*CST_LN3)/1960. - (566757*EILN(6)*CST_LN3)/1120. - (1953125*EILN(10)*CST_LN5)/6048. + e*(lgxOx0p*((-4717987*EIL(1))/1.548288e7 - (771147*EIL(3))/802816. + (61428125*EIL(5))/2.8901376e7 - (549067883*EIL(7))/2.654208e7 - (628919*EILN(1))/2.064384e6 - (4356207*EILN(3))/1.14688e6 + (94221875*EILN(5))/442368. - (252121807*EILN(7))/147456. + (48771934893*EILN(9))/1.605632e7 - (285311670611*EILN(11))/2.60112384e8) + (4717987*EIL(1)*CST_LN2)/2.322432e7 + (628919*EILN(1)*CST_LN2)/3.096576e6 + EILN(3)*((1452069*CST_LN2)/573440. - (1452069*CST_LN3)/573440.) + EIL(3)*((257049*CST_LN2)/401408. - (257049*CST_LN3)/401408.) + EILN(9)*((-16257311631*CST_LN2)/8.02816e6 + (16257311631*CST_LN3)/4.01408e6) + EIL(5)*((-61428125*CST_LN2)/4.3352064e7 + (61428125*CST_LN5)/4.3352064e7) + EILN(5)*((-94221875*CST_LN2)/663552. + (94221875*CST_LN5)/663552.) + e*(lgxOx0p*((-295669*EIL(2))/564480. - (51109*EIL(4))/35280. + (99387*EIL(6))/15680. - (37376*EIL(8))/1225. - (29333*EILN(2))/57600. - (95731*EILN(4))/5040. + (19219113*EILN(6))/31360. - (8075776*EILN(8))/2205. + (330078125*EILN(10))/56448. - (2729376*EILN(12))/1225.) - (51109*EIL(4)*CST_LN2)/52920. - (149504*EIL(8)*CST_LN2)/3675. - (95731*EILN(4)*CST_LN2)/7560. - (32303104*EILN(8)*CST_LN2)/6615. + EILN(12)*((-1819584*CST_LN2)/1225. - (1819584*CST_LN3)/1225.) + (33129*EIL(6)*CST_LN3)/7840. + (6406371*EILN(6)*CST_LN3)/15680. + (330078125*EILN(10)*CST_LN5)/84672.) + EILN(7)*((252121807*CST_LN2)/221184. - (252121807*CST_LN7)/221184.) + EIL(7)*((549067883*CST_LN2)/3.981312e7 - (549067883*CST_LN7)/3.981312e7) + EILN(11)*((285311670611*CST_LN2)/3.90168576e8 - (285311670611*CST_LN11)/3.90168576e8))))))))))))
        +(CST_2PI*((-107*EILN(2))/42. + e*((-473*EILN(1))/112. + (27*EILN(3))/16. + e*((-615*EILN(2))/28. + (470*EILN(4))/21. + e*((-3439*EIL(1))/8064. + (1601*EILN(1))/2688. - (61155*EILN(3))/896. + (613625*EILN(5))/8064. + e*((-97*EIL(2))/84. + (4895*EILN(2))/672. - (1236*EILN(4))/7. + (42741*EILN(6))/224. + e*((-6193*EIL(1))/43008. - (158841*EIL(3))/71680. - (9587*EILN(1))/64512. + (281583*EILN(3))/7168. - (17652125*EILN(5))/43008. + (38035613*EILN(7))/92160. + e*((-503*EIL(2))/1440. - (28387*EIL(4))/7560. - (4283*EILN(2))/1728. + (70037*EILN(4))/504. - (1982853*EILN(6))/2240. + (771671*EILN(8))/945. + e*((-79187*EIL(1))/5.16096e6 - (218511*EIL(3))/573440. - (14407375*EIL(5))/2.408448e6 - (90025*EILN(1))/1.032192e6 - (1453383*EILN(3))/114688. + (410416625*EILN(5))/1.032192e6 - (1329989017*EILN(7))/737280. + (6107613759*EILN(9))/4.01408e6 + e*((-21103*EIL(2))/120960. + EIL(4)/90. - (2306259*EIL(6))/250880. - (2977*EILN(2))/2688. - (203479*EILN(4))/3780. + (17926893*EILN(6))/17920. - (7751524*EILN(8))/2205. + (75062125*EILN(10))/27648. + e*((5870273*EIL(1))/7.4317824e7 - (13101939*EIL(3))/3.211264e7 + (836124125*EIL(5))/6.93633024e8 - (8791122781*EIL(7))/6.3700992e8 - (15352091*EILN(1))/2.4772608e8 - (8015787*EILN(3))/9.17504e6 - (2005202125*EILN(5))/1.0616832e7 + (40803447923*EILN(7))/1.769472e7 - (169746271461*EILN(9))/2.5690112e7 + (146328958845311*EILN(11))/3.121348608e10 + e*((47753*EIL(2))/3.38688e6 - (643151*EIL(4))/846720. + (1910133*EIL(6))/501760. - (64726667*EIL(8))/3.1752e6 - (389417*EILN(2))/345600. + (31149*EILN(4))/4480. - (142812909*EILN(6))/250880. + (1055270873*EILN(8))/211680. - (97981054625*EILN(10))/8.128512e6 + (618544449*EILN(12))/78400.)))))))))))
        +I*((-107*lgxOx0p*EILN(2))/7. + e*(lgxOx0p*((-1419*EILN(1))/56. + (81*EILN(3))/8.) + EILN(1)*(-4.5 + (473*CST_LN2)/28.) + EILN(3)*(-13.5 - (27*CST_LN2)/4. + (27*CST_LN3)/4.) + e*(lgxOx0p*((-1845*EILN(2))/14. + (940*EILN(4))/7.) + EILN(4)*(-48 + (1880*CST_LN2)/21.) + e*(lgxOx0p*((1601*EILN(1))/448. - (183465*EILN(3))/448. + (613625*EILN(5))/1344.) + EILN(1)*(-2.0625 - (1601*CST_LN2)/672.) + EIL(1)*(-1.3125 + (3439*lgxOx0p)/1344. - (3439*CST_LN2)/2016.) + EILN(3)*(18.5625 + (61155*CST_LN2)/224. - (61155*CST_LN3)/224.) + EILN(5)*(-117.1875 - (613625*CST_LN2)/2016. + (613625*CST_LN5)/2016.) + e*((-3*EIL(2))/2. + lgxOx0p*((97*EIL(2))/14. + (4895*EILN(2))/112. - (7416*EILN(4))/7. + (128223*EILN(6))/112.) + EILN(4)*(72 - (4944*CST_LN2)/7.) + EILN(6)*(-243 + (42741*CST_LN3)/56.) + e*(lgxOx0p*((476523*EIL(3))/35840. - (9587*EILN(1))/10752. + (844749*EILN(3))/3584. - (17652125*EILN(5))/7168. + (38035613*EILN(7))/15360.) + EIL(1)*(-1.86328125 + (6193*lgxOx0p)/7168. - (6193*CST_LN2)/10752.) + EILN(1)*(-2.0234375 + (9587*CST_LN2)/16128.) + EIL(3)*(-1.79296875 - (158841*CST_LN2)/17920. + (158841*CST_LN3)/17920.) + EILN(3)*(-4.0078125 - (281583*CST_LN2)/1792. + (281583*CST_LN3)/1792.) + EILN(5)*(197.75390625 + (17652125*CST_LN2)/10752. - (17652125*CST_LN5)/10752.) + EILN(7)*(-459.56640625 - (38035613*CST_LN2)/23040. + (38035613*CST_LN7)/23040.) + e*((-41*EIL(2))/20. + lgxOx0p*((503*EIL(2))/240. + (28387*EIL(4))/1260. - (4283*EILN(2))/288. + (70037*EILN(4))/84. - (5948559*EILN(6))/1120. + (1543342*EILN(8))/315.) + EIL(4)*(-2.2 + (28387*CST_LN2)/1890.) + EILN(4)*(-29 + (70037*CST_LN2)/126.) + EILN(8)*(-819.2 + (6173368*CST_LN2)/945.) + EILN(6)*(461.7 - (1982853*CST_LN3)/560.) + e*(lgxOx0p*((655533*EIL(3))/286720. + (14407375*EIL(5))/401408. - (90025*EILN(1))/172032. - (4360149*EILN(3))/57344. + (410416625*EILN(5))/172032. - (1329989017*EILN(7))/122880. + (18322841277*EILN(9))/2.00704e6) + EIL(1)*(-2.18291015625 + (79187*lgxOx0p)/860160. - (79187*CST_LN2)/1.29024e6) + EILN(1)*(-1.91259765625 + (90025*CST_LN2)/258048.) + EILN(3)*(2.30712890625 + (1453383*CST_LN2)/28672. - (1453383*CST_LN3)/28672.) + EIL(3)*(-2.25439453125 - (218511*CST_LN2)/143360. + (218511*CST_LN3)/143360.) + EILN(9)*(-1401.26044921875 - (6107613759*CST_LN2)/1.00352e6 + (6107613759*CST_LN3)/501760.) + EIL(5)*(-2.74658203125 - (14407375*CST_LN2)/602112. + (14407375*CST_LN5)/602112.) + EILN(5)*(-110.77880859375 - (410416625*CST_LN2)/258048. + (410416625*CST_LN5)/258048.) + EILN(7)*(976.57861328125 + (1329989017*CST_LN2)/184320. - (1329989017*CST_LN7)/184320.) + e*((-1157*EIL(2))/480. + lgxOx0p*((21103*EIL(2))/20160. - EIL(4)/15. + (6918777*EIL(6))/125440. - (2977*EILN(2))/448. - (203479*EILN(4))/630. + (53780679*EILN(6))/8960. - (15503048*EILN(8))/735. + (75062125*EILN(10))/4608.) + EILN(8)*(1930.9714285714285 - (62012192*CST_LN2)/2205.) + EILN(4)*(10.233333333333333 - (203479*CST_LN2)/945.) + EIL(4)*(-2.442857142857143 - (2*CST_LN2)/45.) + EIL(6)*(-3.4714285714285715 + (2306259*CST_LN3)/62720.) + EILN(6)*(-325.0125 + (17926893*CST_LN3)/4480.) + EILN(10)*(-2325.1488095238096 + (75062125*CST_LN5)/6912.) + e*(lgxOx0p*((39305817*EIL(3))/1.605632e7 - (836124125*EIL(5))/1.15605504e8 + (8791122781*EIL(7))/1.0616832e8 - (15352091*EILN(1))/4.128768e7 - (24047361*EILN(3))/4.58752e6 - (2005202125*EILN(5))/1.769472e6 + (40803447923*EILN(7))/2.94912e6 - (509238814383*EILN(9))/1.2845056e7 + (146328958845311*EILN(11))/5.20224768e9) + EILN(1)*(-1.8427144368489583 + (15352091*CST_LN2)/6.193152e7) + EIL(1)*(-2.392769368489583 - (5870273*lgxOx0p)/1.2386304e7 + (5870273*CST_LN2)/1.8579456e7) + EILN(9)*(3634.5192901611326 + (169746271461*CST_LN2)/6.422528e6 - (169746271461*CST_LN3)/3.211264e6) + EILN(3)*(1.232501220703125 + (8015787*CST_LN2)/2.29376e6 - (8015787*CST_LN3)/2.29376e6) + EIL(3)*(-2.6600254603794644 - (13101939*CST_LN2)/8.02816e6 + (13101939*CST_LN3)/8.02816e6) + e*((-17791*EIL(2))/6720. + lgxOx0p*((-47753*EIL(2))/564480. + (643151*EIL(4))/141120. - (5730399*EIL(6))/250880. + (64726667*EIL(8))/529200. - (389417*EILN(2))/57600. + (93447*EILN(4))/2240. - (428438727*EILN(6))/125440. + (1055270873*EILN(8))/35280. - (97981054625*EILN(10))/1.354752e6 + (1855633347*EILN(12))/39200.) + EIL(4)*(-2.9160714285714286 + (643151*CST_LN2)/211680.) + EILN(4)*(1.2375 + (31149*CST_LN2)/1120.) + EIL(8)*(-5.688888888888889 + (64726667*CST_LN2)/396900.) + EILN(8)*(-1879.7714285714285 + (1055270873*CST_LN2)/26460.) + EILN(6)*(138.20625 - (142812909*CST_LN3)/62720.) + EIL(6)*(-2.6035714285714286 - (1910133*CST_LN3)/125440.) + EILN(12)*(-5998.628571428571 + (618544449*CST_LN2)/19600. + (618544449*CST_LN3)/19600.) + EILN(10)*(6587.921626984127 - (97981054625*CST_LN5)/2.032128e6)) + EILN(5)*(40.575663248697914 + (2005202125*CST_LN2)/2.654208e6 - (2005202125*CST_LN5)/2.654208e6) + EIL(5)*(-2.578099568684896 + (836124125*CST_LN2)/1.73408256e8 - (836124125*CST_LN5)/1.73408256e8) + EIL(7)*(-4.428113810221354 - (8791122781*CST_LN2)/1.5925248e8 + (8791122781*CST_LN7)/1.5925248e8) + EILN(7)*(-820.9962361653646 - (40803447923*CST_LN2)/4.42368e6 + (40803447923*CST_LN7)/4.42368e6) + EILN(11)*(-3769.273245820545 - (146328958845311*CST_LN2)/7.80337152e9 + (146328958845311*CST_LN11)/7.80337152e9))))))))))));

    COMPLEX16 PN6_tail = (((-23.358866213151927 + (2*CST_PISQ)/3.)*EILN(2) + lgxOx0p*((-428*EILN(2))/35. - 18*lgxOx0p*EILN(2)) + e*(EILN(1)*(4.379787414965986 - CST_PISQ/8. + CST_LN2*(-1.5285714285714285 + (3*CST_LN2)/2.)) + lgxOx0p*(lgxOx0p*((27*EILN(1))/8. - (729*EILN(3))/8.) + EILN(1)*(2.2928571428571427 - (9*CST_LN2)/2.) + EILN(3)*(-61.90714285714286 + (243*CST_LN2)/2. - (243*CST_LN3)/2.)) + EILN(3)*(-118.25426020408163 + (27*CST_PISQ)/8. + (-41.27142857142857 - (81*CST_LN3)/2.)*CST_LN3 + CST_LN2*(41.27142857142857 - (81*CST_LN2)/2. + 81*CST_LN3)) + e*((58.39716553287982 - (5*CST_PISQ)/3.)*EILN(2) + lgxOx0p*((214*EILN(2))/7. + lgxOx0p*(45*EILN(2) - 288*EILN(4)) + EILN(4)*(-195.65714285714284 - 384*CST_LN2)) + EILN(4)*(-373.74185941043083 + (32*CST_PISQ)/3. + (-130.43809523809523 - 128*CST_LN2)*CST_LN2) + e*(EILN(1)*(-2.3723848497732427 + (13*CST_PISQ)/192. + (0.8279761904761904 - (13*CST_LN2)/16.)*CST_LN2) + EIL(1)*(0.4258126653439153 - (7*CST_PISQ)/576. + lgxOx0p*(0.22291666666666668 + (21*lgxOx0p)/64. - (7*CST_LN2)/16.) + (-0.1486111111111111 + (7*CST_LN2)/48.)*CST_LN2) + EILN(3)*(280.85386798469386 - (513*CST_PISQ)/64. + CST_LN2*(-98.01964285714286 + (1539*CST_LN2)/16. - (1539*CST_LN3)/8.) + CST_LN3*(98.01964285714286 + (1539*CST_LN3)/16.)) + lgxOx0p*(lgxOx0p*((-117*EILN(1))/64. + (13851*EILN(3))/64. - (46875*EILN(5))/64.) + EILN(1)*(-1.2419642857142856 + (39*CST_LN2)/16.) + EILN(3)*(147.02946428571428 - (4617*CST_LN2)/16. + (4617*CST_LN3)/16.) + EILN(5)*(-497.58184523809524 + (15625*CST_LN2)/16. - (15625*CST_LN5)/16.)) + EILN(5)*(-950.4746994283825 + (15625*CST_PISQ)/576. + (-331.7212301587302 - (15625*CST_LN5)/48.)*CST_LN5 + CST_LN2*(331.7212301587302 - (15625*CST_LN2)/48. + (15625*CST_LN5)/24.)) + e*((1.4599291383219954 - CST_PISQ/24.)*EIL(2) + (-33.578370181405894 + (23*CST_PISQ)/24.)*EILN(2) + EILN(4)*(934.3546485260771 - (80*CST_PISQ)/3. + CST_LN2*(326.0952380952381 + 320*CST_LN2)) + lgxOx0p*((107*EIL(2))/140. - (2461*EILN(2))/140. + lgxOx0p*((9*EIL(2))/8. - (207*EILN(2))/8. + 720*EILN(4) - (6561*EILN(6))/4.) + EILN(4)*(489.14285714285717 + 960*CST_LN2) + EILN(6)*(-1114.3285714285714 - 2187*CST_LN3)) + EILN(6)*(-2128.5766836734692 + (243*CST_PISQ)/4. + (-742.8857142857142 - 729*CST_LN3)*CST_LN3) + e*(EILN(1)*(-0.0380189879771353 + (5*CST_PISQ)/4608. + (0.013268849206349206 - (5*CST_LN2)/384.)*CST_LN2) + EIL(1)*(0.1786892434925359 - (47*CST_PISQ)/9216. + lgxOx0p*(0.0935453869047619 + (141*lgxOx0p)/1024. - (47*CST_LN2)/256.) + (-0.06236359126984127 + (47*CST_LN2)/768.)*CST_LN2) + EIL(3)*(3.1411287866709183 - (459*CST_PISQ)/5120. + CST_LN2*(-1.0962723214285715 + (1377*CST_LN2)/1280. - (1377*CST_LN3)/640.) + CST_LN3*(1.0962723214285715 + (1377*CST_LN3)/1280.)) + EILN(3)*(-197.70634127869897 + (2889*CST_PISQ)/512. + (-69.00066964285715 - (8667*CST_LN3)/128.)*CST_LN3 + CST_LN2*(69.00066964285715 - (8667*CST_LN2)/128. + (8667*CST_LN3)/64.)) + EILN(5)*(2554.400754713778 - (671875*CST_PISQ)/9216. + CST_LN2*(-891.5008060515873 + (671875*CST_LN2)/768. - (671875*CST_LN5)/384.) + CST_LN5*(891.5008060515873 + (671875*CST_LN5)/768.)) + lgxOx0p*(lgxOx0p*((12393*EIL(3))/5120. - (15*EILN(1))/512. - (78003*EILN(3))/512. + (2015625*EILN(5))/1024. - (17294403*EILN(7))/5120.) + EILN(1)*(-0.019903273809523808 + (5*CST_LN2)/128.) + EILN(3)*(-103.50100446428571 + (26001*CST_LN2)/128. - (26001*CST_LN3)/128.) + EIL(3)*(1.644408482142857 - (4131*CST_LN2)/1280. + (4131*CST_LN3)/1280.) + EILN(5)*(1337.251209077381 - (671875*CST_LN2)/256. + (671875*CST_LN5)/256.) + EILN(7)*(-2294.7682552083334 + (5764801*CST_LN2)/1280. - (5764801*CST_LN7)/1280.)) + EILN(7)*(-4383.437998191551 + (5764801*CST_PISQ)/46080. + (-1529.8455034722222 - (5764801*CST_LN7)/3840.)*CST_LN7 + CST_LN2*(1529.8455034722222 - (5764801*CST_LN2)/3840. + (5764801*CST_LN7)/1920.)) + e*((0.535307350718065 - (11*CST_PISQ)/720.)*EIL(2) + (5.271966332829428 - (65*CST_PISQ)/432.)*EILN(2) + EILN(8)*(-8504.703645250693 + (32768*CST_PISQ)/135. + (-5936.382645502646 - (524288*CST_LN2)/45.)*CST_LN2) + EILN(4)*(-786.4151625094482 + (202*CST_PISQ)/9. + (-274.4634920634921 - (808*CST_LN2)/3.)*CST_LN2) + EIL(4)*(5.7099450743260265 - (22*CST_PISQ)/135. + CST_LN2*(1.9928042328042328 + (88*CST_LN2)/45.)) + EILN(6)*(6172.872382653061 - (7047*CST_PISQ)/40. + CST_LN3*(2154.3685714285716 + (21141*CST_LN3)/10.)) + lgxOx0p*((1177*EIL(2))/4200. + (1391*EILN(2))/504. + lgxOx0p*((33*EIL(2))/80. + (22*EIL(4))/5. + (65*EILN(2))/16. - 606*EILN(4) + (190269*EILN(6))/40. - (32768*EILN(8))/5.) + EILN(8)*(-4452.286984126984 - (262144*CST_LN2)/15.) + EILN(4)*(-411.6952380952381 - 808*CST_LN2) + EIL(4)*(2.9892063492063494 + (88*CST_LN2)/15.) + EILN(6)*(3231.552857142857 + (63423*CST_LN3)/10.)) + e*(EILN(1)*(-0.1078788783851214 + (227*CST_PISQ)/73728. + (0.037650359623015876 - (227*CST_LN2)/6144.)*CST_LN2) + EIL(1)*(0.10369678970763653 - (1091*CST_PISQ)/368640. + lgxOx0p*(0.05428617931547619 + (3273*lgxOx0p)/40960. - (1091*CST_LN2)/10240.) + (-0.03619078621031746 + (1091*CST_LN2)/30720.)*CST_LN2) + EIL(3)*(0.8083787318638392 - (189*CST_PISQ)/8192. + CST_LN2*(-0.28212890625 + (567*CST_LN2)/2048. - (567*CST_LN3)/1024.) + (0.28212890625 + (567*CST_LN3)/2048.)*CST_LN3) + EILN(3)*(55.31620179468272 - (12933*CST_PISQ)/8192. + CST_LN2*(-19.305678013392857 + (38799*CST_LN2)/2048. - (38799*CST_LN3)/1024.) + CST_LN3*(19.305678013392857 + (38799*CST_LN3)/2048.)) + EILN(9)*(-15781.430264343084 + (129140163*CST_PISQ)/286720. + (-11015.622959980867 - (387420489*CST_LN3)/17920.)*CST_LN3 + CST_LN2*(5507.811479990433 - (387420489*CST_LN2)/71680. + (387420489*CST_LN3)/17920.)) + EIL(5)*(9.547178900508307 - (15625*CST_PISQ)/57344. + CST_LN2*(-3.3320212850765305 + (46875*CST_LN2)/14336. - (46875*CST_LN5)/7168.) + CST_LN5*(3.3320212850765305 + (46875*CST_LN5)/14336.)) + EILN(5)*(-2502.421669588788 + (5265625*CST_PISQ)/73728. + (-873.3598012772817 - (5265625*CST_LN5)/6144.)*CST_LN5 + CST_LN2*(873.3598012772817 - (5265625*CST_LN2)/6144. + (5265625*CST_LN5)/3072.)) + EILN(7)*(13698.243744348596 - (28824005*CST_PISQ)/73728. + CST_LN2*(-4780.767198350694 + (28824005*CST_LN2)/6144. - (28824005*CST_LN7)/3072.) + CST_LN7*(4780.767198350694 + (28824005*CST_LN7)/6144.)) + lgxOx0p*(lgxOx0p*((5103*EIL(3))/8192. + (421875*EIL(5))/57344. - (681*EILN(1))/8192. + (349191*EILN(3))/8192. - (15796875*EILN(5))/8192. + (86472015*EILN(7))/8192. - (3486784401*EILN(9))/286720.) + EILN(1)*(-0.05647553943452381 + (227*CST_LN2)/2048.) + EILN(9)*(-8261.71721998565 + (1162261467*CST_LN2)/71680. - (1162261467*CST_LN3)/35840.) + EIL(3)*(0.423193359375 - (1701*CST_LN2)/2048. + (1701*CST_LN3)/2048.) + EILN(3)*(28.958517020089285 - (116397*CST_LN2)/2048. + (116397*CST_LN3)/2048.) + EILN(5)*(-1310.0397019159227 + (5265625*CST_LN2)/2048. - (5265625*CST_LN5)/2048.) + EIL(5)*(4.998031927614796 - (140625*CST_LN2)/14336. + (140625*CST_LN5)/14336.) + EILN(7)*(7171.150797526042 - (28824005*CST_LN2)/2048. + (28824005*CST_LN7)/2048.)) + e*((0.350788529069035 - (173*CST_PISQ)/17280.)*EIL(2) + (-0.8617637274817335 + (85*CST_PISQ)/3456.)*EILN(2) + EIL(4)*(0.630318612100925 - (17*CST_PISQ)/945. + (0.2199848828420257 + (68*CST_LN2)/315.)*CST_LN2) + EILN(4)*(305.4820614764424 - (1177*CST_PISQ)/135. + CST_LN2*(106.61502645502645 + (4708*CST_LN2)/45.)) + EILN(8)*(28551.505094770182 - (770048*CST_PISQ)/945. + CST_LN2*(19929.284595616024 + (12320768*CST_LN2)/315.)) + EILN(6)*(-6891.267013392857 + (62937*CST_PISQ)/320. + (-2405.0925 - (188811*CST_LN3)/80.)*CST_LN3) + EIL(6)*(15.20411916909621 - (243*CST_PISQ)/560. + CST_LN3*(5.306326530612245 + (729*CST_LN3)/140.)) + lgxOx0p*((18511*EIL(2))/100800. - (1819*EILN(2))/4032. + lgxOx0p*((173*EIL(2))/640. + (17*EIL(4))/35. + (6561*EIL(6))/560. - (85*EILN(2))/128. + (1177*EILN(4))/5. - (1699299*EILN(6))/320. + (770048*EILN(8))/35. - (9765625*EILN(10))/448.) + EIL(4)*(0.32997732426303855 + (68*CST_LN2)/105.) + EILN(4)*(159.92253968253968 + (4708*CST_LN2)/15.) + EILN(8)*(14946.963446712018 + (6160384*CST_LN2)/105.) + EILN(6)*(-3607.63875 - (566433*CST_LN3)/80.) + EIL(6)*(7.959489795918367 + (2187*CST_LN3)/140.) + EILN(10)*(-14808.983489229026 - (9765625*CST_LN5)/336.)) + EILN(10)*(-28287.937482987574 + (9765625*CST_PISQ)/12096. + (-9872.655659486016 - (9765625*CST_LN5)/1008.)*CST_LN5) + e*(EILN(1)*(-0.0680163655222198 + (34349*CST_PISQ)/1.769472e7 + (0.023738109447337962 - (34349*CST_LN2)/1.47456e6)*CST_LN2) + EIL(1)*(0.06808435086530393 - (10315*CST_PISQ)/5.308416e6 + lgxOx0p*(0.03564275509465939 + (10315*lgxOx0p)/196608. - (10315*CST_LN2)/147456.) + (-0.023761836729772928 + (10315*CST_LN2)/442368.)*CST_LN2) + EIL(3)*(0.710630895408875 - (46521*CST_PISQ)/2.29376e6 + CST_LN2*(-0.24801433952487245 + (139563*CST_LN2)/573440. - (139563*CST_LN3)/286720.) + (0.24801433952487245 + (139563*CST_LN3)/573440.)*CST_LN3) + EILN(3)*(-9.413281625864457 + (176067*CST_PISQ)/655360. + (-3.285290352957589 - (528201*CST_LN3)/163840.)*CST_LN3 + CST_LN2*(3.285290352957589 - (528201*CST_LN2)/163840. + (528201*CST_LN3)/81920.)) + EILN(9)*(56714.51501248296 - (2970223749*CST_PISQ)/1.835008e6 + CST_LN2*(-19793.69750621562 + (8910671247*CST_LN2)/458752. - (8910671247*CST_LN3)/114688.) + CST_LN3*(39587.39501243124 + (8910671247*CST_LN3)/114688.)) + EIL(5)*(-0.5856487057024771 + (828125*CST_PISQ)/4.9545216e7 + CST_LN2*(0.20439482420029645 - (828125*CST_LN2)/4.128768e6 + (828125*CST_LN5)/2.064384e6) + (-0.20439482420029645 - (828125*CST_LN5)/4.128768e6)*CST_LN5) + EILN(5)*(1227.593353656254 - (185984375*CST_PISQ)/5.308416e6 + CST_LN2*(-428.4372615646013 + (185984375*CST_LN2)/442368. - (185984375*CST_LN5)/221184.) + CST_LN5*(428.4372615646013 + (185984375*CST_LN5)/442368.)) + e*((0.23072094418124392 - (59*CST_PISQ)/8960.)*EIL(2) + (-0.20418731142920132 + (1007*CST_PISQ)/172800.)*EILN(2) + EILN(8)*(-39562.05892121081 + (1067008*CST_PISQ)/945. + (-27614.779984882844 - (17072128*CST_LN2)/315.)*CST_LN2) + EILN(4)*(-70.04415576971529 + (2159*CST_PISQ)/1080. + (-24.445820105820104 - (2159*CST_LN2)/90.)*CST_LN2) + EIL(4)*(1.2281943544613612 - (53*CST_PISQ)/1512. + (0.4286470143613001 + (53*CST_LN2)/126.)*CST_LN2) + EIL(8)*(35.43626518854455 - (2048*CST_PISQ)/2025. + CST_LN2*(24.734927689594358 + (32768*CST_LN2)/675.)) + EILN(12)*(-84072.6973574344 + (419904*CST_PISQ)/175. + CST_LN2*(-29341.863183673468 - (5038848*CST_LN2)/175. - (10077696*CST_LN3)/175.) + (-29341.863183673468 - (5038848*CST_LN3)/175.)*CST_LN3) + EIL(6)*(-3.8010297922740524 + (243*CST_PISQ)/2240. + (-1.3265816326530613 - (729*CST_LN3)/560.)*CST_LN3) + EILN(6)*(4057.599303252551 - (14823*CST_PISQ)/128. + CST_LN3*(1416.1258928571428 + (44469*CST_LN3)/32.)) + EILN(10)*(108437.0936847857 - (224609375*CST_PISQ)/72576. + CST_LN5*(37845.18002802973 + (224609375*CST_LN5)/6048.)) + lgxOx0p*((18939*EIL(2))/156800. - (107749*EILN(2))/1.008e6 + lgxOx0p*((1593*EIL(2))/8960. + (53*EIL(4))/56. - (6561*EIL(6))/2240. + (2048*EIL(8))/75. - (1007*EILN(2))/6400. - (2159*EILN(4))/40. + (400221*EILN(6))/128. - (1067008*EILN(8))/35. + (224609375*EILN(10))/2688. - (11337408*EILN(12))/175.) + EILN(8)*(-20711.084988662133 - (8536064*CST_LN2)/105.) + EILN(4)*(-36.668730158730156 - (2159*CST_LN2)/30.) + EIL(4)*(0.6429705215419501 + (53*CST_LN2)/42.) + EIL(8)*(18.551195767195768 + (16384*CST_LN2)/225.) + EILN(12)*(-44012.79477551021 - (15116544*CST_LN2)/175. - (15116544*CST_LN3)/175.) + EIL(6)*(-1.9898724489795918 - (2187*CST_LN3)/560.) + EILN(6)*(2124.1888392857145 + (133407*CST_LN3)/32.) + EILN(10)*(56767.7700420446 + (224609375*CST_LN5)/2016.))) + EIL(7)*(23.464584191708244 - (213297637*CST_PISQ)/3.1850496e8 + CST_LN2*(-8.189277145322949 + (213297637*CST_LN2)/2.654208e7 - (213297637*CST_LN7)/1.327104e7) + CST_LN7*(8.189277145322949 + (213297637*CST_LN7)/2.654208e7)) + EILN(7)*(-17145.635086676328 + (4329365551*CST_PISQ)/8.84736e6 + (-5983.9269432689525 - (4329365551*CST_LN7)/737280.)*CST_LN7 + CST_LN2*(5983.9269432689525 - (4329365551*CST_LN2)/737280. + (4329365551*CST_LN7)/368640.)) + lgxOx0p*(lgxOx0p*((1256067*EIL(3))/2.29376e6 - (828125*EIL(5))/1.835008e6 + (213297637*EIL(7))/1.179648e7 - (34349*EILN(1))/655360. - (4753809*EILN(3))/655360. + (185984375*EILN(5))/196608. - (4329365551*EILN(7))/327680. + (80196041223*EILN(9))/1.835008e6 - (3138428376721*EILN(11))/8.257536e7) + EILN(1)*(-0.03560716417100694 + (34349*CST_LN2)/491520.) + EILN(3)*(-4.927935529436384 + (1584603*CST_LN2)/163840. - (1584603*CST_LN3)/163840.) + EIL(3)*(0.37202150928730865 - (418689*CST_LN2)/573440. + (418689*CST_LN3)/573440.) + EILN(9)*(29690.546259323433 - (26732013741*CST_LN2)/458752. + (26732013741*CST_LN3)/229376.) + EIL(5)*(-0.30659223630044463 + (828125*CST_LN2)/1.376256e6 - (828125*CST_LN5)/1.376256e6) + EILN(5)*(642.6558923469019 - (185984375*CST_LN2)/147456. + (185984375*CST_LN5)/147456.) + EILN(7)*(-8975.89041490343 + (4329365551*CST_LN2)/245760. - (4329365551*CST_LN7)/245760.) + EIL(7)*(12.283915717984422 - (213297637*CST_LN2)/8.84736e6 + (213297637*CST_LN7)/8.84736e6) + EILN(11)*(-25820.518896105077 + (3138428376721*CST_LN2)/6.193152e7 - (3138428376721*CST_LN11)/6.193152e7)) + EILN(11)*(-49322.036508621015 + (3138428376721*CST_PISQ)/2.22953472e9 + (-17213.67926407005 - (3138428376721*CST_LN11)/1.8579456e8)*CST_LN11 + CST_LN2*(17213.67926407005 - (3138428376721*CST_LN2)/1.8579456e8 + (3138428376721*CST_LN11)/9.289728e7))))))))))))
        +I*((428*CST_PI*EILN(2))/105. + 12*lgxOx0p*CST_PI*EILN(2) + e*(lgxOx0p*((-9*CST_PI*EILN(1))/4. + (243*CST_PI*EILN(3))/4.) + EILN(1)*((-107*CST_PI)/140. + (3*CST_PI*CST_LN2)/2.) + EILN(3)*((2889*CST_PI)/140. - (81*CST_PI*CST_LN2)/2. + (81*CST_PI*CST_LN3)/2.) + e*((-214*CST_PI*EILN(2))/21. + lgxOx0p*(-30*CST_PI*EILN(2) + 192*CST_PI*EILN(4)) + EILN(4)*((6848*CST_PI)/105. + 128*CST_PI*CST_LN2) + e*(lgxOx0p*((39*CST_PI*EILN(1))/32. - (4617*CST_PI*EILN(3))/32. + (15625*CST_PI*EILN(5))/32.) + EILN(1)*((1391*CST_PI)/3360. - (13*CST_PI*CST_LN2)/16.) + EIL(1)*((107*CST_PI)/1440. + (7*lgxOx0p*CST_PI)/32. - (7*CST_PI*CST_LN2)/48.) + EILN(3)*((-54891*CST_PI)/1120. + (1539*CST_PI*CST_LN2)/16. - (1539*CST_PI*CST_LN3)/16.) + EILN(5)*((334375*CST_PI)/2016. - (15625*CST_PI*CST_LN2)/48. + (15625*CST_PI*CST_LN5)/48.) + e*((107*CST_PI*EIL(2))/420. + (2461*CST_PI*EILN(2))/420. + lgxOx0p*((3*CST_PI*EIL(2))/4. + (69*CST_PI*EILN(2))/4. - 480*CST_PI*EILN(4) + (2187*CST_PI*EILN(6))/2.) + EILN(4)*((-3424*CST_PI)/21. - 320*CST_PI*CST_LN2) + EILN(6)*((26001*CST_PI)/70. + 729*CST_PI*CST_LN3) + e*(lgxOx0p*((4131*CST_PI*EIL(3))/2560. + (5*CST_PI*EILN(1))/256. + (26001*CST_PI*EILN(3))/256. - (671875*CST_PI*EILN(5))/512. + (5764801*CST_PI*EILN(7))/2560.) + EIL(1)*((5029*CST_PI)/161280. + (47*lgxOx0p*CST_PI)/512. - (47*CST_PI*CST_LN2)/768.) + EILN(1)*((107*CST_PI)/16128. - (5*CST_PI*CST_LN2)/384.) + EIL(3)*((49113*CST_PI)/89600. - (1377*CST_PI*CST_LN2)/1280. + (1377*CST_PI*CST_LN3)/1280.) + EILN(3)*((309123*CST_PI)/8960. - (8667*CST_PI*CST_LN2)/128. + (8667*CST_PI*CST_LN3)/128.) + EILN(5)*((-14378125*CST_PI)/32256. + (671875*CST_PI*CST_LN2)/768. - (671875*CST_PI*CST_LN5)/768.) + EILN(7)*((88119101*CST_PI)/115200. - (5764801*CST_PI*CST_LN2)/3840. + (5764801*CST_PI*CST_LN7)/3840.) + e*((1177*CST_PI*EIL(2))/12600. - (1391*CST_PI*EILN(2))/1512. + lgxOx0p*((11*CST_PI*EIL(2))/40. + (44*CST_PI*EIL(4))/15. - (65*CST_PI*EILN(2))/24. + 404*CST_PI*EILN(4) - (63423*CST_PI*EILN(6))/20. + (65536*CST_PI*EILN(8))/15.) + EIL(4)*((4708*CST_PI)/4725. + (88*CST_PI*CST_LN2)/45.) + EILN(4)*((43228*CST_PI)/315. + (808*CST_PI*CST_LN2)/3.) + EILN(8)*((7012352*CST_PI)/4725. + (262144*CST_PI*CST_LN2)/45.) + EILN(6)*((-754029*CST_PI)/700. - (21141*CST_PI*CST_LN3)/10.) + e*(lgxOx0p*((1701*CST_PI*EIL(3))/4096. + (140625*CST_PI*EIL(5))/28672. + (227*CST_PI*EILN(1))/4096. - (116397*CST_PI*EILN(3))/4096. + (5265625*CST_PI*EILN(5))/4096. - (28824005*CST_PI*EILN(7))/4096. + (1162261467*CST_PI*EILN(9))/143360.) + EILN(1)*((24289*CST_PI)/1.29024e6 - (227*CST_PI*CST_LN2)/6144.) + EIL(1)*((116737*CST_PI)/6.4512e6 + (1091*lgxOx0p*CST_PI)/20480. - (1091*CST_PI*CST_LN2)/30720.) + EILN(3)*((-1383831*CST_PI)/143360. + (38799*CST_PI*CST_LN2)/2048. - (38799*CST_PI*CST_LN3)/2048.) + EIL(3)*((2889*CST_PI)/20480. - (567*CST_PI*CST_LN2)/2048. + (567*CST_PI*CST_LN3)/2048.) + EILN(9)*((13817997441*CST_PI)/5.0176e6 - (387420489*CST_PI*CST_LN2)/71680. + (387420489*CST_PI*CST_LN3)/35840.) + EIL(5)*((334375*CST_PI)/200704. - (46875*CST_PI*CST_LN2)/14336. + (46875*CST_PI*CST_LN5)/14336.) + EILN(5)*((112684375*CST_PI)/258048. - (5265625*CST_PI*CST_LN2)/6144. + (5265625*CST_PI*CST_LN5)/6144.) + EILN(7)*((-88119101*CST_PI)/36864. + (28824005*CST_PI*CST_LN2)/6144. - (28824005*CST_PI*CST_LN7)/6144.) + e*((18511*CST_PI*EIL(2))/302400. + (1819*CST_PI*EILN(2))/12096. + lgxOx0p*((173*CST_PI*EIL(2))/960. + (34*CST_PI*EIL(4))/105. + (2187*CST_PI*EIL(6))/280. + (85*CST_PI*EILN(2))/192. - (2354*CST_PI*EILN(4))/15. + (566433*CST_PI*EILN(6))/160. - (1540096*CST_PI*EILN(8))/105. + (9765625*CST_PI*EILN(10))/672.) + EILN(8)*((-164790272*CST_PI)/33075. - (6160384*CST_PI*CST_LN2)/315.) + EILN(4)*((-251878*CST_PI)/4725. - (4708*CST_PI*CST_LN2)/45.) + EIL(4)*((3638*CST_PI)/33075. + (68*CST_PI*CST_LN2)/315.) + EIL(6)*((26001*CST_PI)/9800. + (729*CST_PI*CST_LN3)/140.) + EILN(6)*((962037*CST_PI)/800. + (188811*CST_PI*CST_LN3)/80.) + EILN(10)*((208984375*CST_PI)/42336. + (9765625*CST_PI*CST_LN5)/1008.) + e*(lgxOx0p*((418689*CST_PI*EIL(3))/1.14688e6 - (828125*CST_PI*EIL(5))/2.752512e6 + (213297637*CST_PI*EIL(7))/1.769472e7 + (34349*CST_PI*EILN(1))/983040. + (1584603*CST_PI*EILN(3))/327680. - (185984375*CST_PI*EILN(5))/294912. + (4329365551*CST_PI*EILN(7))/491520. - (26732013741*CST_PI*EILN(9))/917504. + (3138428376721*CST_PI*EILN(11))/1.2386304e8) + EIL(1)*((220741*CST_PI)/1.8579456e7 + (10315*lgxOx0p*CST_PI)/294912. - (10315*CST_PI*CST_LN2)/442368.) + EILN(1)*((525049*CST_PI)/4.42368e7 - (34349*CST_PI*CST_LN2)/1.47456e6) + EILN(9)*((-317813941143*CST_PI)/3.211264e7 + (8910671247*CST_PI*CST_LN2)/458752. - (8910671247*CST_PI*CST_LN3)/229376.) + EIL(3)*((4977747*CST_PI)/4.01408e7 - (139563*CST_PI*CST_LN2)/573440. + (139563*CST_PI*CST_LN3)/573440.) + EILN(3)*((18839169*CST_PI)/1.14688e7 - (528201*CST_PI*CST_LN2)/163840. + (528201*CST_PI*CST_LN3)/163840.) + EILN(5)*((-3980065625*CST_PI)/1.8579456e7 + (185984375*CST_PI*CST_LN2)/442368. - (185984375*CST_PI*CST_LN5)/442368.) + EIL(5)*((-17721875*CST_PI)/1.73408256e8 + (828125*CST_PI*CST_LN2)/4.128768e6 - (828125*CST_PI*CST_LN5)/4.128768e6) + e*((6313*CST_PI*EIL(2))/156800. + (107749*CST_PI*EILN(2))/3.024e6 + lgxOx0p*((531*CST_PI*EIL(2))/4480. + (53*CST_PI*EIL(4))/84. - (2187*CST_PI*EIL(6))/1120. + (4096*CST_PI*EIL(8))/225. + (1007*CST_PI*EILN(2))/9600. + (2159*CST_PI*EILN(4))/60. - (133407*CST_PI*EILN(6))/64. + (2134016*CST_PI*EILN(8))/105. - (224609375*CST_PI*EILN(10))/4032. + (7558272*CST_PI*EILN(12))/175.) + EIL(4)*((5671*CST_PI)/26460. + (53*CST_PI*CST_LN2)/126.) + EILN(4)*((231013*CST_PI)/18900. + (2159*CST_PI*CST_LN2)/90.) + EIL(8)*((438272*CST_PI)/70875. + (16384*CST_PI*CST_LN2)/675.) + EILN(8)*((228339712*CST_PI)/33075. + (8536064*CST_PI*CST_LN2)/315.) + EILN(6)*((-1586061*CST_PI)/2240. - (44469*CST_PI*CST_LN3)/32.) + EIL(6)*((-26001*CST_PI)/39200. - (729*CST_PI*CST_LN3)/560.) + EILN(12)*((89859456*CST_PI)/6125. + (5038848*CST_PI*CST_LN2)/175. + (5038848*CST_PI*CST_LN3)/175.) + EILN(10)*((-4806640625*CST_PI)/254016. - (224609375*CST_PI*CST_LN5)/6048.)) + EIL(7)*((3260406737*CST_PI)/7.962624e8 - (213297637*CST_PI*CST_LN2)/2.654208e7 + (213297637*CST_PI*CST_LN7)/2.654208e7) + EILN(7)*((66177444851*CST_PI)/2.21184e7 - (4329365551*CST_PI*CST_LN2)/737280. + (4329365551*CST_PI*CST_LN7)/737280.) + EILN(11)*((335811836309147*CST_PI)/3.90168576e10 - (3138428376721*CST_PI*CST_LN2)/1.8579456e8 + (3138428376721*CST_PI*CST_LN11)/1.8579456e8))))))))))));


#undef EIL
#undef EILN
    return PN0 + v2*(PN2 + v*(PN3_tail + v*(PN4 + v*(I*ImPN5 + PN5_tail + v*(PN6 + PN6_tail)))));
}



// h21
COMPLEX16 calculate_h21FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe)
{
    COMPLEX16 ret;
    INT pm = p - 1;
    if (pm==0) return 0.0;
    INT pm2 = pm*pm;
    INT pm3 = pm2*pm;
    INT pm4 = pm2*pm2;
    INT pm5 = pm4*pm;
    REAL8 *ePow = dc->ePow;
    REAL8 *sq1m2Pow = dc->sq1mesqPow;
    REAL8 m1p3eta = -1. + 3.*eta;
    REAL8 y = sq1m2Pow[1];
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 dm = sqrt(1. - 4.*eta);
#define EVALEJ(pp, qq, aa, bb) evaluate_Jpqab((pp), (qq), (aa), (bb), pe)
// #define EVALEJ(p, q, a, b) 1.0
#define SQ1M2(n) sq1m2Pow[n]
#define EPOW(n) ePow[n]
    REAL8 ImPN1, ImPN3, ImPN5, PN6;
    ImPN1 = -EVALEJ(pm,pm,1,0) - pm*EVALEJ(pm,pm,0,0)*SQ1M2(1) + EVALEJ(pm,pm,2,0)*SQ1M2(2);
    ImPN1 *= 2.*sq1m2Pow[1]*dm/(3.*ePow[1]);

    ImPN3 = EPOW(1)*EVALEJ(-1 + pm,pm,1,0) - 2*EVALEJ(pm,pm,1,1) - 
        EPOW(1)*EVALEJ(1 + pm,pm,1,0) - 2*pm*EVALEJ(pm,pm,0,1)*SQ1M2(1) + 
        2*EVALEJ(pm,pm,2,1)*SQ1M2(2) + 
        EVALEJ(pm,pm,0,0)*(2*SQ1M2(1) - 
        (pm*SQ1M2(1)*(28*(-4 + eta) + 42*SQ1M2(1) + (33 - 10*eta)*SQ1M2(2)))/21.) + 
        (EVALEJ(pm,pm,2,0)*SQ1M2(2)*(28*(-13 + eta) + 84*SQ1M2(1) - 
        4*(25 + 6*eta)*SQ1M2(2) + (27 + 16*eta)*pm*SQ1M2(3)))/42. + 
        (EVALEJ(pm,pm,1,0)*(-28*(-4 + eta) - 84*SQ1M2(1) + (46 - 8*eta)*SQ1M2(2) - 
        (111 + 16*eta)*pm*SQ1M2(3)))/42. + 
        ((303 + 80*eta)*EVALEJ(pm,pm,3,0)*SQ1M2(4))/42. - 
        ((27 + 16*eta)*EVALEJ(pm,pm,4,0)*SQ1M2(6))/14.;    
    ImPN3 *= dm/(ePow[1]*y);

    ImPN5 = 3*EPOW(1)*EVALEJ(-1 + pm,pm,1,1) - 3*EVALEJ(pm,pm,1,2) - 
        3*EPOW(1)*EVALEJ(1 + pm,pm,1,1) - 3*EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*SQ1M2(1) - 
        3*pm*EVALEJ(pm,pm,0,2)*SQ1M2(1) + 3*EVALEJ(pm,pm,2,2)*SQ1M2(2) + 
        (EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*
        (399 - 294*eta + 252*SQ1M2(1) + (-495 + 206*eta)*SQ1M2(2)))/84. - 
        (EPOW(1)*EVALEJ(1 + pm,pm,1,0)*
        (399 - 294*eta + 252*SQ1M2(1) + (-495 + 206*eta)*SQ1M2(2)))/84. + 
        EVALEJ(pm,pm,1,1)*(-9.5 + 7*eta - 6*SQ1M2(1) - 
        ((-495 + 206*eta)*SQ1M2(2))/42. + ((-181 + 12*eta)*pm*SQ1M2(3))/14.) - 
        (EVALEJ(pm,pm,2,1)*SQ1M2(2)*(357 + 294*eta - 252*SQ1M2(1) + 
        (657 - 110*eta)*SQ1M2(2) + 3*(-97 + 12*eta)*pm*SQ1M2(3)))/42. + 
        ((303 + 80*eta)*EVALEJ(pm,pm,3,1)*SQ1M2(4))/14. + 
        (EVALEJ(pm,pm,3,0)*(-46322 - 8905*eta + 1245*eta2 + 378*(71 + 8*eta)*SQ1M2(1) - 
        3*(4441 + 383*eta + 139*eta2)*SQ1M2(2) + 
        (3938 + 1636*eta + 15*eta2)*pm*SQ1M2(3))*SQ1M2(4))/756. + 
        EVALEJ(pm,pm,0,1)*(6*SQ1M2(1) - 
        (pm*SQ1M2(1)*(63 - 210*eta + 252*SQ1M2(1) + (-159 + 122*eta)*SQ1M2(2)))/
        42. + (-5 + 2*eta)*pm2*SQ1M2(4)) + 
        (EVALEJ(pm,pm,0,0)*(1512 + 126*(9 - 10*eta + 18*pm)*SQ1M2(1) + 
        126*(6 + (25 + 6*eta)*pm)*SQ1M2(2) - 
        2*(477 - 366*eta + (435 + 506*eta + 76*eta2)*pm)*SQ1M2(3) + 
        3*pm*(1314 - 7*eta2*pm + 5*eta*(-44 + 21*pm))*SQ1M2(4) + 
        (-551 + 614*eta - 16*eta2)*pm*SQ1M2(5)))/252. - 
        (EVALEJ(pm,pm,2,0)*SQ1M2(2)*(252*(49 - 50*eta + eta2) + 
        108*(338 + 3*eta)*SQ1M2(1) + (-30974 + 6314*eta - 234*eta2)*SQ1M2(2) - 
        18*(-402 - 22*eta + (-933 - 152*eta + 32*eta2)*pm)*SQ1M2(3) - 
        3*(-911 + 610*eta + 36*eta2 + 18*(27 + 16*eta)*pm)*SQ1M2(4) + 
        3*(-233 - 358*eta + 104*eta2)*pm*SQ1M2(5)))/756. - 
        (3*(27 + 16*eta)*EVALEJ(pm,pm,4,1)*SQ1M2(6))/14. - 
        (EVALEJ(pm,pm,4,0)*(-56510 - 4216*eta + 264*eta2 + 162*(27 + 16*eta)*SQ1M2(1) - 
        3*(-785 + 1088*eta + 54*eta2)*SQ1M2(2) + 
        3*(10 + 68*eta + 3*eta2)*pm*SQ1M2(3))*SQ1M2(6))/756. + 
        (EVALEJ(pm,pm,1,0)*(252*(35 - 10*eta + 2*eta2) + 756*(7 + 6*eta)*SQ1M2(1) - 
        24*(246 - 22*eta + 61*eta2)*SQ1M2(2) - 
        36*(5*(-207 + 22*eta) + (-462 - 240*eta + 32*eta2)*pm)*SQ1M2(3) - 
        6*(-361 + 66*eta + 8*eta2 + 9*(303 + 80*eta)*pm)*SQ1M2(4) + 
        2*(37 - 706*eta + 306*eta2)*pm*SQ1M2(5) + 63*(-15 + eta)*eta*pm2*SQ1M2(6)))
        /1512. - ((15902 + 3784*eta + 357*eta2)*EVALEJ(pm,pm,5,0)*SQ1M2(8))/756. + 
        (5*(10 + 68*eta + 3*eta2)*EVALEJ(pm,pm,6,0)*SQ1M2(10))/252.;   
    ImPN5 *= dm/(ePow[1]*sq1m2Pow[3]);


    PN6 = (3*eta*pm*EVALEJ(pm,pm,1,0))/5. + 
        (eta*EVALEJ(pm,pm,3,0)*SQ1M2(1)*(-5786 + 3695*pm*SQ1M2(1)))/315. + 
        (2*eta*EVALEJ(pm,pm,2,0)*(567*SQ1M2(1) + pm*(95 - 2799*SQ1M2(2))))/945. - 
        (5431*eta*EVALEJ(pm,pm,5,0)*SQ1M2(3))/63. - 
        (eta*EVALEJ(pm,pm,4,0)*SQ1M2(1)*(-9865 - 10323*SQ1M2(2) + 2475*pm*SQ1M2(3)))/
        315. + (275*eta*EVALEJ(pm,pm,6,0)*SQ1M2(5))/7.;
    PN6 *= dm * y / ePow[1];

    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];

    // tail part
    COMPLEX16 PN4_tail, PN6_tail;
    COMPLEX16 lnIpvO2 =  pm>0 ? (log(0.5*pm) + 0.5*I*CST_PI) : (log(-0.5*pm) - 0.5*I*CST_PI);
    REAL8 lnxOx0p = log(v2/CONST_X0P);
    COMPLEX16 pref_tail;
    pref_tail = -1. + 6.*lnxOx0p + 4*lnIpvO2;
    PN4_tail = -(pm*EVALEJ(pm,pm,2,0)*SQ1M2(1)) + EVALEJ(pm,pm,3,0)*(-3 + (5*pm*SQ1M2(1))/2.) - 
        (35*EVALEJ(pm,pm,5,0)*SQ1M2(2))/2. + 
        EVALEJ(pm,pm,4,0)*(7 + 6*SQ1M2(2) - (3*pm*SQ1M2(3))/2.) + 
        (15*EVALEJ(pm,pm,6,0)*SQ1M2(4))/2.;
    PN4_tail *= dm * pref_tail * y / (3.*ePow[1]*pm);

    PN6_tail = pref_tail * (
        (3*pm2*EVALEJ(pm,pm,1,0))/2. - pm2*EVALEJ(pm,pm,2,1)*SQ1M2(1) + 
        (pm*EVALEJ(pm,pm,3,1)*(-6 + 5*pm*SQ1M2(1)))/2. - 
        (35*pm*EVALEJ(pm,pm,5,1)*SQ1M2(2))/2. + 
        EVALEJ(pm,pm,3,0)*(-3 + (pm*(-14*(-7 + eta) - 175*SQ1M2(1) + 
        (23 + 3*eta)*SQ1M2(2)))/14. - 
        (pm2*SQ1M2(1)*(202 - 40*eta - 105*SQ1M2(1) + (48 + 38*eta)*SQ1M2(2)))/24.) + 
        EVALEJ(pm,pm,4,1)*(pm*(7 + 6*SQ1M2(2)) - (3*pm2*SQ1M2(3))/2.) - 
        (EVALEJ(pm,pm,5,0)*SQ1M2(2)*(420 - 
        2*pm*(760 - 70*eta - 282*SQ1M2(1) + (652 + 191*eta)*SQ1M2(2)) + 
        (231 + 80*eta)*pm2*SQ1M2(3)))/24. + 
        (pm*EVALEJ(pm,pm,2,0)*(84*SQ1M2(1) + 
        pm*(-98 + (154 - 28*eta)*SQ1M2(1) - 84*SQ1M2(2) + (-33 + 17*eta)*SQ1M2(3)))
        )/42. + (15*pm*EVALEJ(pm,pm,6,1)*SQ1M2(4))/2. + 
        (EVALEJ(pm,pm,6,0)*(1260 - 2*pm*
        (9186 + 1550*eta - 630*SQ1M2(1) + 15*(131 + 67*eta)*SQ1M2(2)) + 
        15*(27 + 16*eta)*pm2*SQ1M2(3))*SQ1M2(4))/168. + 
        EVALEJ(pm,pm,4,0)*(7 + 6*SQ1M2(2) + 
        (pm2*(837 + 20*eta - 90*SQ1M2(1) + 15*(3 + 7*eta)*SQ1M2(2))*SQ1M2(3))/60. + 
        (pm*(98*(-7 + eta) + 735*SQ1M2(1) + 7*(-167 + 9*eta)*SQ1M2(2) + 
        378*SQ1M2(3) - 6*(50 + 19*eta)*SQ1M2(4)))/42.) + 
        (3*(199 + 80*eta)*pm*EVALEJ(pm,pm,7,0)*SQ1M2(6))/8. - 
        (5*(27 + 16*eta)*pm*EVALEJ(pm,pm,8,0)*SQ1M2(8))/8.
    ) + (
        4*pm*(1 + pm)*EVALEJ(pm,pm,2,0)*SQ1M2(1) - 
        2*(1 + pm)*EVALEJ(pm,pm,3,0)*(-6 + 5*pm*SQ1M2(1)) + 
        70*(1 + pm)*EVALEJ(pm,pm,5,0)*SQ1M2(2) + 
        2*(1 + pm)*EVALEJ(pm,pm,4,0)*(-14 - 12*SQ1M2(2) + 3*pm*SQ1M2(3)) - 
        30*(1 + pm)*EVALEJ(pm,pm,6,0)*SQ1M2(4) - 
        (eta*pm*EVALEJ(pm,pm,7,0)*(3832*CST_LN2 - 960*(4*CST_LN2) + (8*CST_LN2))*SQ1M2(6))/32.   
    );
    PN6_tail *= dm / (ePow[1]*pm2*y);

    // print_debug("PN1 = %.16e\n", v*ImPN1);
    // print_debug("PN3 = %.16e\n", dc->vPow[3]*ImPN3);
    // print_debug("PN4_tail = %.16e, %.16e\n", dc->vPow[4]*creal(PN4_tail), dc->vPow[4]*cimag(PN4_tail) );
    // print_debug("ImPN5 = %.16e\n", dc->vPow[5]*ImPN5);
    // print_debug("PN6 = %.16e\n", dc->vPow[6]*PN6);
    // print_debug("PN6_tail = %.16e, %.16e\n", dc->vPow[6]*creal(PN6_tail), dc->vPow[6]*cimag(PN6_tail) );
    // return v2*PN0;
    return  v*(I*ImPN1 + 
                v2 * (I*ImPN3 + 
                    v*(PN4_tail + 
                        v*(I*ImPN5 + 
                            v*(PN6 + PN6_tail)))));
#undef EVALEJ
#undef SQ1M2
#undef EPOW
}


COMPLEX16 calculate_h21_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 dm = sqrt(1. - 4.*eta);
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];
    REAL8 e = dc->ePow[1], e2 = dc->ePow[2];
    COMPLEX16 eilvec[50], eilnvec[50];
    eilvec[0] = 1.0;
    eilnvec[0] = 1.0;
    for (size_t p=1; p<50; p++) {
        eilvec[p] = cexp(I*p*l);
        eilnvec[p] = cexp(-I*p*l);
    }
#define EIL(x) eilvec[x]
#define EILN(x) eilnvec[x]
    COMPLEX16 ImPN1 = ((2*EILN(1))/3. + e*((4*EILN(2))/3. + e*(EIL(1)/12. - (2*EILN(1))/3. + (9*EILN(3))/4. + e*(EIL(2)/9. - (5*EILN(2))/3. + (32*EILN(4))/9. + e*(-0.013888888888888888*EIL(1) + (9*EIL(3))/64. + (7*EILN(1))/96. - (27*EILN(3))/8. + (3125*EILN(5))/576. + e*(-0.027777777777777776*EIL(2) + (8*EIL(4))/45. + (17*EILN(2))/36. - (56*EILN(4))/9. + (81*EILN(6))/10. + e*((-37*EIL(1))/4608. - (9*EIL(3))/160. + (3125*EIL(5))/13824. - (5*EILN(1))/432. + (765*EILN(3))/512. - (3125*EILN(5))/288. + (823543*EILN(7))/69120. + e*(-0.005555555555555556*EIL(2) - (14*EIL(4))/135. + (81*EIL(6))/280. - (7*EILN(2))/108. + (166*EILN(4))/45. - (729*EILN(6))/40. + (16384*EILN(8))/945. + e*((-857*EIL(1))/138240. + (9*EIL(3))/4096. - (34375*EIL(5))/193536. + (823543*EIL(7))/2.21184e6 - (1355*EILN(1))/221184. - (1629*EILN(3))/5120. + (884375*EILN(5))/110592. - (823543*EILN(7))/27648. + (14348907*EILN(9))/573440. + e*((-169*EIL(2))/25920. + (19*EIL(4))/945. - (81*EIL(6))/280. + (4096*EIL(8))/8505. - (23*EILN(2))/8640. - (449*EILN(4))/405. + (17901*EILN(6))/1120. - (45056*EILN(8))/945. + (1953125*EILN(10))/54432. + e*((-125923*EIL(1))/2.654208e7 - (981*EIL(3))/143360. + (1409375*EIL(5))/2.4772608e7 - (9058973*EIL(7))/1.990656e7 + (14348907*EIL(9))/2.29376e7 - (26239*EILN(1))/5.5296e6 + (10359*EILN(3))/327680. - (7253125*EILN(5))/2.322432e6 + (534479407*EILN(7))/1.769472e7 - (43046721*EILN(9))/573440. + (285311670611*EILN(11))/5.5738368e9)))))))))));

    COMPLEX16 ImPN3 =eta*((10*EILN(1))/21. + e*((4*EILN(2))/21. + e*((-5*EIL(1))/28. + (5*EILN(1))/7. - (27*EILN(3))/28. + e*((-25*EIL(2))/63. + (53*EILN(2))/21. - (32*EILN(4))/9. + e*((19*EIL(1))/252. - (45*EIL(3))/64. - (89*EILN(1))/672. + (27*EILN(3))/4. - (34375*EILN(5))/4032. + e*((17*EIL(2))/84. - (8*EIL(4))/7. - (5*EILN(2))/4. + (328*EILN(4))/21. - (243*EILN(6))/14. + e*((1639*EIL(1))/32256. + (2169*EIL(3))/4480. - (171875*EIL(5))/96768. + (767*EILN(1))/12096. - (17235*EILN(3))/3584. + (265625*EILN(5))/8064. - (2235331*EILN(7))/69120. + e*((53*EIL(2))/1260. + (194*EIL(4))/189. - (1053*EIL(6))/392. + (215*EILN(2))/756. - (874*EILN(4))/63. + (3645*EILN(6))/56. - (376832*EILN(8))/6615. + e*((1013*EIL(1))/20160. - (4779*EIL(3))/143360. + (3125*EIL(5))/1568. - (588245*EIL(7))/147456. + (8585*EILN(1))/172032. + (81*EILN(3))/56. - (8871875*EILN(5))/258048. + (117649*EILN(7))/960. - (387420489*EILN(9))/4.01408e6 + e*((1633*EIL(2))/25920. - (1723*EIL(4))/6615. + (891*EIL(6))/245. - (69632*EIL(8))/11907. + (517*EILN(2))/12096. + (15737*EILN(4))/2835. - (606447*EILN(6))/7840. + (1470464*EILN(8))/6615. - (60546875*EILN(10))/381024. + e*((8803877*EIL(1))/1.8579456e8 + (4527*EIL(3))/57344. - (139478125*EIL(5))/1.73408256e8 + (252592403*EIL(7))/3.981312e7 - (272629233*EIL(9))/3.211264e7 + (733739*EILN(1))/1.548288e7 - (329949*EILN(3))/2.29376e6 + (561128125*EILN(5))/3.2514048e7 - (575421259*EILN(7))/3.538944e6 + (3142410633*EILN(9))/8.02816e6 - (285311670611*EILN(11))/1.11476736e9)))))))))))
        +((-17*EILN(1))/42. + e*((40*EILN(2))/21. + e*((737*EIL(1))/336. - (275*EILN(1))/84. + (621*EILN(3))/112. + e*((515*EIL(2))/126. - (277*EILN(2))/42. + (98*EILN(4))/9. + e*((235*EIL(1))/1008. + (1701*EIL(3))/256. - (3175*EILN(1))/2688. - (189*EILN(3))/16. + (298175*EILN(5))/16128. + e*((-571*EIL(2))/1008. + (12907*EIL(4))/1260. - (25*EILN(2))/18. - (5027*EILN(4))/252. + (16263*EILN(6))/560. + e*((43307*EIL(1))/129024. - (46233*EIL(3))/17920. + (5941265*EIL(5))/387072. - (42239*EILN(1))/48384. - (14607*EILN(3))/14336. - (1037945*EILN(5))/32256. + (11998777*EILN(7))/276480. + e*((10741*EIL(2))/20160. - (9911*EIL(4))/1512. + (354267*EIL(6))/15680. - (12115*EILN(2))/12096. + (127*EILN(4))/2520. - (22329*EILN(6))/448. + (827597*EILN(8))/13230. + e*((129601*EIL(1))/483840. + (644031*EIL(3))/573440. - (9276215*EIL(5))/677376. + (290661623*EIL(7))/8.84736e6 - (4612277*EILN(1))/6.193152e6 - (10881*EILN(3))/17920. + (5450965*EILN(5))/3.096576e6 - (5136103*EILN(7))/69120. + (1406749437*EILN(9))/1.605632e7 + e*((33791*EIL(2))/103680. + (36527*EIL(4))/13230. - (648819*EIL(6))/25088. + (90201113*EIL(8))/1.90512e6 - (239453*EILN(2))/241920. + (4069*EILN(4))/4536. + (213417*EILN(6))/62720. - (22574521*EILN(8))/211680. + (729739175*EILN(10))/6.096384e6 + e*((170384141*EIL(1))/7.4317824e8 + (66933*EIL(3))/229376. + (4607313955*EIL(5))/6.93633024e8 - (1461904115*EIL(7))/3.1850496e7 + (8705357847*EIL(9))/1.2845056e8 - (40569527*EILN(1))/6.193152e7 - (2348217*EILN(3))/1.835008e6 + (648154735*EILN(5))/1.30056192e8 + (38599981*EILN(7))/1.4155776e7 - (4719397581*EILN(9))/3.211264e7 + (712498890301*EILN(11))/4.45906944e9)))))))))));

    COMPLEX16 ImPN5 = eta2*((79*EILN(1))/252. + e*((-4*EILN(2))/63. + e*((307*EIL(1))/2016. - (23*EILN(1))/84. - (279*EILN(3))/224. + e*((17*EIL(2))/27. - (43*EILN(2))/63. - (704*EILN(4))/189. + e*((-2035*EIL(1))/12096. + (5853*EIL(3))/3584. - (151*EILN(1))/16128. - (387*EILN(3))/448. - (803125*EILN(5))/96768. + e*((-53*EIL(2))/84. + (472*EIL(4))/135. + (319*EILN(2))/756. - (8*EILN(4))/63. - (81*EILN(6))/5. + e*((-115567*EIL(1))/774144. - (4209*EIL(3))/2240. + (15659375*EIL(5))/2.322432e6 - (13*EILN(1))/81. + (45069*EILN(3))/28672. + (34375*EILN(5))/12096. - (48589037*EILN(7))/1.65888e6 + e*((-1073*EIL(2))/7560. - (13432*EIL(4))/2835. + (5967*EIL(6))/490. - (1591*EILN(2))/4536. + (3544*EILN(4))/945. + (729*EILN(6))/70. - (999424*EILN(8))/19845. + e*((-4518683*EIL(1))/2.322432e7 + (265329*EIL(3))/1.14688e6 - (115459375*EIL(5))/1.0838016e7 + (1116136063*EIL(7))/5.308416e7 - (801949*EILN(1))/4.128768e6 - (278181*EILN(3))/286720. + (132415625*EILN(5))/1.8579456e7 + (29529899*EILN(7))/1.10592e6 - (2683245609*EILN(9))/3.211264e7 + e*((-158849*EIL(2))/544320. + (32561*EIL(4))/19845. - (34641*EIL(6))/1568. + (6275072*EIL(8))/178605. - (49471*EILN(2))/181440. - (665*EILN(4))/243. + (44739*EILN(6))/3920. + (1163264*EILN(8))/19845. - (154296875*EILN(10))/1.143072e6 + e*((-1019268601*EIL(1))/4.45906944e9 - (3529227*EIL(3))/8.02816e6 + (23090228125*EIL(5))/4.161798144e9 - (20627516819*EIL(7))/4.7775744e8 + (10508182893*EIL(9))/1.835008e8 - (212352671*EILN(1))/9.289728e8 - (4339497*EILN(3))/1.835008e7 - (55365625*EILN(5))/7.962624e6 + (6472930331*EILN(7))/4.2467328e8 + (3773762541*EILN(9))/3.211264e7 - (200003481098311*EILN(11))/9.364045824e11)))))))))))
        +eta*((-509*EILN(1))/189. + e*((-233*EILN(2))/21. + e*((-1115*EIL(1))/189. + (269*EILN(1))/126. - (369*EILN(3))/14. + e*((-32083*EIL(2))/2268. + (1739*EILN(2))/756. - (29852*EILN(4))/567. + e*((-12953*EIL(1))/3024. - (25623*EIL(3))/896. + (56747*EILN(1))/12096. + (1179*EILN(3))/112. - (2332825*EILN(5))/24192. + e*((-1643*EIL(2))/1008. - (120233*EIL(4))/2268. + (13235*EILN(2))/9072. + (27695*EILN(4))/756. - (9381*EILN(6))/56. + e*((-2005781*EIL(1))/290304. + (17193*EIL(3))/1792. - (81019685*EIL(5))/870912. + (1597291*EILN(1))/217728. - (23913*EILN(3))/3584. + (14192485*EILN(5))/145152. - (34923623*EILN(7))/124416. + e*((-895421*EIL(2))/90720. + (213313*EIL(4))/5670. - (246375*EIL(6))/1568. + (1133*EILN(2))/168. - (31679*EILN(4))/1134. + (62523*EILN(6))/280. - (9090479*EILN(8))/19845. + e*((-29215715*EIL(1))/3.483648e6 - (4544721*EIL(3))/286720. + (788588105*EIL(5))/8.128512e6 - (10271781457*EIL(7))/3.981312e7 + (86272961*EILN(1))/9.289728e6 + (530073*EILN(3))/71680. - (1122677555*EILN(5))/1.3934592e7 + (127771469*EILN(7))/276480. - (1176027687*EILN(9))/1.605632e6 + e*((-70694333*EIL(2))/6.53184e6 - (15101813*EIL(4))/476280. + (3331011*EIL(6))/15680. - (29009965*EIL(8))/69984. + (3891301*EILN(2))/435456. + (1392833*EILN(4))/102060. - (16077*EILN(6))/80. + (1703744929*EILN(8))/1.90512e6 - (7902264775*EILN(10))/6.858432e6 + e*((-2710123391*EIL(1))/2.7869184e8 - (9662535*EIL(3))/802816. - (27947445835*EIL(5))/3.90168576e8 + (33815406611*EIL(7))/7.962624e7 - (52553432277*EIL(9))/8.02816e7 + (15312762383*EILN(1))/1.3934592e9 + (2264067*EILN(3))/286720. + (7057642855*EILN(5))/1.95084288e8 - (7291525703*EILN(7))/1.5925248e7 + (26467185033*EILN(9))/1.605632e7 - (5814788384471*EILN(11))/3.2514048e9)))))))))))
        +((-43*EILN(1))/189. + e*((1801*EILN(2))/126. + e*((-929*EIL(1))/216. - (12349*EILN(1))/504. + (1035*EILN(3))/28. + e*((-23605*EIL(2))/4536. - (36223*EILN(2))/1512. + (35186*EILN(4))/567. + e*((-67247*EIL(1))/12096. - (1341*EIL(3))/224. - (255235*EILN(1))/6048. + (81*EILN(3))/448. + (635*EILN(5))/8. + e*((-45445*EIL(2))/6048. - (17378*EIL(4))/2835. - (1320691*EILN(2))/18144. + (9497*EILN(4))/126. + (693*EILN(6))/10. + e*((-487901*EIL(1))/82944. - (406359*EIL(3))/35840. - (307847*EIL(5))/62208. - (41731213*EILN(1))/870912. - (485361*EILN(3))/3584. + (150949193*EILN(5))/580608. - (5086879*EILN(7))/1.24416e6 + e*((-432071*EIL(2))/72576. - (551207*EIL(4))/30240. - (8541*EIL(6))/6272. - (4659727*EILN(2))/72576. - (5149783*EILN(4))/18144. + (2983329*EILN(6))/4480. - (3233257*EILN(8))/15876. + e*((-470468147*EIL(1))/6.967296e7 - (1553661*EIL(3))/286720. - (331525141*EIL(5))/1.0838016e7 + (252372007*EIL(7))/3.981312e7 - (504014671*EILN(1))/9.289728e6 - (18193527*EILN(3))/286720. - (8756028007*EILN(5))/1.3934592e7 + (4969571443*EILN(7))/3.31776e6 - (5104107243*EILN(9))/8.02816e6 + e*((-89964809*EIL(2))/1.306368e7 - (613931*EIL(4))/238140. - (26387199*EIL(6))/501760. + (713712583*EIL(8))/3.429216e7 - (92475863*EILN(2))/1.24416e6 - (11333683*EILN(4))/816480. - (21806163*EILN(6))/15680. + (11881441387*EILN(8))/3.81024e6 - (323912653165*EILN(10))/2.19469824e8 + e*((-8545821089*EIL(1))/1.11476736e9 - (24052185*EIL(3))/3.211264e6 + (89206955*EIL(5))/1.5925248e7 - (87144460621*EIL(7))/9.5551488e8 + (14818419009*EIL(9))/3.211264e8 - (333133833491*EILN(1))/5.5738368e9 - (210803193*EILN(3))/2.29376e6 + (44735618315*EILN(5))/2.60112384e8 - (191103902893*EILN(7))/6.3700992e7 + (56456249301*EILN(9))/9.17504e6 - (22003051499627*EILN(11))/7.3156608e9)))))))))));


    COMPLEX16 PN6 = eta*((25*EILN(1))/3. + e*((13049*EILN(2))/315. + e*((33577*EIL(1))/2520. + (3617*EILN(1))/315. + (328687*EILN(3))/2520. + e*((3131*EIL(2))/140. - (7597*EILN(2))/420. + (208703*EILN(4))/630. + e*((12655*EIL(1))/432. + (1491859*EIL(3))/40320. + (232273*EILN(1))/6720. - (768119*EILN(3))/5040. + (2561959*EILN(5))/3456. + e*((38813*EIL(2))/1080. + (32477*EIL(4))/540. + (195275*EILN(2))/3024. - (2059321*EILN(4))/3780. + (7688777*EILN(6))/5040. + e*((3227593*EIL(1))/64512. + (4054571*EIL(3))/100800. + (93041831*EIL(5))/967680. + (1624277*EILN(1))/30240. + (3676417*EILN(3))/21504. - (30245053*EILN(5))/20160. + (227161799*EILN(7))/76800. + e*((6016309*EIL(2))/100800. + (17155843*EIL(4))/453600. + (2181401*EIL(6))/14400. + (11629519*EILN(2))/181440. + (25915189*EILN(4))/50400. - (40269151*EILN(6))/11200. + (1243946051*EILN(8))/226800. + e*((2142077761*EIL(1))/2.90304e7 + (914702587*EIL(3))/1.29024e7 + (846200363*EIL(5))/4.064256e7 + (15638909569*EIL(7))/6.63552e7 + (720039545*EILN(1))/9.289728e6 + (165594083*EILN(3))/3.2256e6 + (34262255411*EILN(5))/2.322432e7 - (32684158571*EILN(7))/4.1472e6 + (1183080744779*EILN(9))/1.204224e8 + e*((1899453*EIL(2))/22400. + (2840603*EIL(4))/33075. - (6934769*EIL(6))/282240. + (13820924471*EIL(8))/3.81024e7 + (14004287*EILN(2))/151200. - (3016529*EILN(4))/45360. + (10958452709*EILN(6))/2.8224e6 - (22885180819*EILN(8))/1.4112e6 + (52177448933*EILN(10))/3.048192e6 + e*((559261387421*EIL(1))/5.5738368e9 + (1739627627*EIL(3))/1.806336e7 + (38183596189*EIL(5))/3.46816512e8 - (508563239267*EIL(7))/4.1803776e9 + (1598430340453*EIL(9))/2.8901376e9 + (1610851709*EILN(1))/1.548288e7 + (23202580073*EILN(3))/2.064384e8 - (52292343767*EILN(5))/9.7542144e7 + (2336445461417*EILN(7))/2.4772608e8 - (11495867067709*EILN(9))/3.612672e8 + (6825819824023387*EILN(11))/2.341011456e11)))))))))));

    REAL8 lgxOx0p = log(v2/CONST_X0P);
                    
    COMPLEX16 ImPN4_tail = (CST_2PI*(EILN(1)/6. + e*((2*EILN(2))/3. + e*(EIL(1)/48. - EILN(1)/6. + (27*EILN(3))/16. + e*(EIL(2)/18. - (5*EILN(2))/6. + (32*EILN(4))/9. + e*(-0.003472222222222222*EIL(1) + (27*EIL(3))/256. + (7*EILN(1))/384. - (81*EILN(3))/32. + (15625*EILN(5))/2304. + e*(-0.013888888888888888*EIL(2) + (8*EIL(4))/45. + (17*EILN(2))/72. - (56*EILN(4))/9. + (243*EILN(6))/20. + e*((-37*EIL(1))/18432. - (27*EIL(3))/640. + (15625*EIL(5))/55296. - (5*EILN(1))/1728. + (2295*EILN(3))/2048. - (15625*EILN(5))/1152. + (5764801*EILN(7))/276480. + e*(-0.002777777777777778*EIL(2) - (14*EIL(4))/135. + (243*EIL(6))/560. - (7*EILN(2))/216. + (166*EILN(4))/45. - (2187*EILN(6))/80. + (32768*EILN(8))/945. + e*((-857*EIL(1))/552960. + (27*EIL(3))/16384. - (171875*EIL(5))/774144. + (5764801*EIL(7))/8.84736e6 - (1355*EILN(1))/884736. - (4887*EILN(3))/20480. + (4421875*EILN(5))/442368. - (5764801*EILN(7))/110592. + (129140163*EILN(9))/2.29376e6 + e*((-169*EIL(2))/51840. + (19*EIL(4))/945. - (243*EIL(6))/560. + (8192*EIL(8))/8505. - (23*EILN(2))/17280. - (449*EILN(4))/405. + (53703*EILN(6))/2240. - (90112*EILN(8))/945. + (9765625*EILN(10))/108864. + e*((-125923*EIL(1))/1.0616832e8 - (2943*EIL(3))/573440. + (7046875*EIL(5))/9.9090432e7 - (63412811*EIL(7))/7.962624e7 + (129140163*EIL(9))/9.17504e7 - (26239*EILN(1))/2.21184e7 + (31077*EILN(3))/1.31072e6 - (36265625*EILN(5))/9.289728e6 + (3741355849*EILN(7))/7.077888e7 - (387420489*EILN(9))/2.29376e6 + (3138428376721*EILN(11))/2.22953472e10)))))))))))
        +I*(lgxOx0p*EILN(1) + EILN(1)*(-0.16666666666666666 - (2*CST_LN2)/3.) + e*((-2*EILN(2))/3. + 4*lgxOx0p*EILN(2) + e*(lgxOx0p*(-EILN(1) + (81*EILN(3))/8.) + EIL(1)*(0.020833333333333332 - lgxOx0p/8. + CST_LN2/12.) + EILN(1)*(0.16666666666666666 + (2*CST_LN2)/3.) + EILN(3)*(-1.6875 - (27*CST_LN2)/4. + (27*CST_LN3)/4.) + e*(EIL(2)/18. + (5*EILN(2))/6. + lgxOx0p*(-0.3333333333333333*EIL(2) - 5*EILN(2) + (64*EILN(4))/3.) + EILN(4)*(-3.5555555555555554 + (128*CST_LN2)/9.) + e*(lgxOx0p*((-81*EIL(3))/128. + (7*EILN(1))/64. - (243*EILN(3))/16. + (15625*EILN(5))/384.) + EILN(1)*(-0.018229166666666668 - (7*CST_LN2)/96.) + EIL(1)*(-0.003472222222222222 + lgxOx0p/48. - CST_LN2/72.) + EILN(3)*(2.53125 + (81*CST_LN2)/8. - (81*CST_LN3)/8.) + EIL(3)*(0.10546875 + (27*CST_LN2)/64. - (27*CST_LN3)/64.) + EILN(5)*(-6.781684027777778 - (15625*CST_LN2)/576. + (15625*CST_LN5)/576.) + e*(-0.013888888888888888*EIL(2) - (17*EILN(2))/72. + lgxOx0p*(EIL(2)/12. - (16*EIL(4))/15. + (17*EILN(2))/12. - (112*EILN(4))/3. + (729*EILN(6))/10.) + EILN(4)*(6.222222222222222 - (224*CST_LN2)/9.) + EIL(4)*(0.17777777777777778 - (32*CST_LN2)/45.) + EILN(6)*(-12.15 + (243*CST_LN3)/5.) + e*(lgxOx0p*((81*EIL(3))/320. - (15625*EIL(5))/9216. - (5*EILN(1))/288. + (6885*EILN(3))/1024. - (15625*EILN(5))/192. + (5764801*EILN(7))/46080.) + EIL(1)*(-0.002007378472222222 + (37*lgxOx0p)/3072. - (37*CST_LN2)/4608.) + EILN(1)*(0.0028935185185185184 + (5*CST_LN2)/432.) + EIL(3)*(-0.0421875 - (27*CST_LN2)/160. + (27*CST_LN3)/160.) + EILN(3)*(-1.12060546875 - (2295*CST_LN2)/512. + (2295*CST_LN3)/512.) + EILN(5)*(13.563368055555555 + (15625*CST_LN2)/288. - (15625*CST_LN5)/288.) + EIL(5)*(0.28257016782407407 + (15625*CST_LN2)/13824. - (15625*CST_LN5)/13824.) + EILN(7)*(-20.85069806134259 - (5764801*CST_LN2)/69120. + (5764801*CST_LN7)/69120.) + e*(-0.002777777777777778*EIL(2) + (7*EILN(2))/216. + lgxOx0p*(EIL(2)/60. + (28*EIL(4))/45. - (729*EIL(6))/280. - (7*EILN(2))/36. + (332*EILN(4))/15. - (6561*EILN(6))/40. + (65536*EILN(8))/315.) + EIL(4)*(-0.1037037037037037 + (56*CST_LN2)/135.) + EILN(4)*(-3.688888888888889 + (664*CST_LN2)/45.) + EILN(8)*(-34.67513227513228 + (262144*CST_LN2)/945.) + EILN(6)*(27.3375 - (2187*CST_LN3)/20.) + EIL(6)*(0.43392857142857144 - (243*CST_LN3)/140.) + e*(lgxOx0p*((-81*EIL(3))/8192. + (171875*EIL(5))/129024. - (5764801*EIL(7))/1.47456e6 - (1355*EILN(1))/147456. - (14661*EILN(3))/10240. + (4421875*EILN(5))/73728. - (5764801*EILN(7))/18432. + (387420489*EILN(9))/1.14688e6) + EIL(1)*(-0.0015498408564814815 + (857*lgxOx0p)/92160. - (857*CST_LN2)/138240.) + EILN(1)*(0.0015315303096064814 + (1355*CST_LN2)/221184.) + EILN(3)*(0.238623046875 + (4887*CST_LN2)/5120. - (4887*CST_LN3)/5120.) + EIL(3)*(0.00164794921875 + (27*CST_LN2)/4096. - (27*CST_LN3)/4096.) + EILN(9)*(-56.300643048967636 - (129140163*CST_LN2)/573440. + (129140163*CST_LN3)/286720.) + EIL(5)*(-0.2220194175760582 - (171875*CST_LN2)/193536. + (171875*CST_LN5)/193536.) + EILN(5)*(-9.99591968677662 - (4421875*CST_LN2)/110592. + (4421875*CST_LN5)/110592.) + EILN(7)*(52.12674515335648 + (5764801*CST_LN2)/27648. - (5764801*CST_LN7)/27648.) + EIL(7)*(0.651584314416956 + (5764801*CST_LN2)/2.21184e6 - (5764801*CST_LN7)/2.21184e6) + e*((-169*EIL(2))/51840. + (23*EILN(2))/17280. + lgxOx0p*((169*EIL(2))/8640. - (38*EIL(4))/315. + (729*EIL(6))/280. - (16384*EIL(8))/2835. - (23*EILN(2))/2880. - (898*EILN(4))/135. + (161109*EILN(6))/1120. - (180224*EILN(8))/315. + (9765625*EILN(10))/18144.) + EILN(8)*(95.35661375661375 - (720896*CST_LN2)/945.) + EIL(8)*(0.9631981187536743 - (65536*CST_LN2)/8505.) + EILN(4)*(1.108641975308642 - (1796*CST_LN2)/405.) + EIL(4)*(0.020105820105820106 - (76*CST_LN2)/945.) + EIL(6)*(-0.43392857142857144 + (243*CST_LN3)/140.) + EILN(6)*(-23.974553571428572 + (53703*CST_LN3)/560.) + EILN(10)*(-89.70481518224574 + (9765625*CST_LN5)/27216.) + e*(lgxOx0p*((8829*EIL(3))/286720. - (7046875*EIL(5))/1.6515072e7 + (63412811*EIL(7))/1.327104e7 - (387420489*EIL(9))/4.58752e7 - (26239*EILN(1))/3.6864e6 + (93231*EILN(3))/655360. - (36265625*EILN(5))/1.548288e6 + (3741355849*EILN(7))/1.179648e7 - (1162261467*EILN(9))/1.14688e6 + (3138428376721*EILN(11))/3.7158912e9) + EIL(1)*(-0.0011860694414303626 + (125923*lgxOx0p)/1.769472e7 - (125923*CST_LN2)/2.654208e7) + EILN(1)*(0.0011862973813657406 + (26239*CST_LN2)/5.5296e6) + EILN(9)*(168.90192914690292 + (387420489*CST_LN2)/573440. - (387420489*CST_LN3)/286720.) + EIL(9)*(1.407516076224191 + (129140163*CST_LN2)/2.29376e7 - (129140163*CST_LN3)/1.14688e7) + EIL(3)*(-0.005132184709821429 - (2943*CST_LN2)/143360. + (2943*CST_LN3)/143360.) + EILN(3)*(-0.023709869384765624 - (31077*CST_LN2)/327680. + (31077*CST_LN3)/327680.) + EILN(5)*(3.903841425712357 + (36265625*CST_LN2)/2.322432e6 - (36265625*CST_LN5)/2.322432e6) + EIL(5)*(0.07111559469233114 + (7046875*CST_LN2)/2.4772608e7 - (7046875*CST_LN5)/2.4772608e7) + EIL(7)*(-0.7963808287318351 - (63412811*CST_LN2)/1.990656e7 + (63412811*CST_LN7)/1.990656e7) + EILN(7)*(-52.859777507075556 - (3741355849*CST_LN2)/1.769472e7 + (3741355849*CST_LN7)/1.769472e7) + EILN(11)*(-140.7660687482364 - (3138428376721*CST_LN2)/5.5738368e9 + (3138428376721*CST_LN11)/5.5738368e9)))))))))))));

    COMPLEX16 ImPN6_tail = eta*(CST_2PI*(EILN(1)/28. + e*((-5*EILN(2))/21. + e*((-37*EIL(1))/672. + (11*EILN(1))/42. - (351*EILN(3))/224. + e*((-19*EIL(2))/84. + (47*EILN(2))/28. - (16*EILN(4))/3. + e*((83*EIL(1))/4032. - (297*EIL(3))/512. - (227*EILN(1))/5376. + (405*EILN(3))/64. - (453125*EILN(5))/32256. + e*((109*EIL(2))/1008. - (388*EIL(4))/315. - (107*EILN(2))/144. + (1180*EILN(4))/63. - (8991*EILN(6))/280. + e*((393*EIL(1))/28672. + (1377*EIL(3))/3584. - (203125*EIL(5))/86016. + (31*EILN(1))/1792. - (119475*EILN(3))/28672. + (171875*EILN(5))/3584. - (823543*EILN(7))/12288. + e*((113*EIL(2))/5040. + (1019*EIL(4))/945. - (33291*EIL(6))/7840. + (479*EILN(2))/3024. - (4951*EILN(4))/315. + (124659*EILN(6))/1120. - (868352*EILN(8))/6615. + e*((103247*EIL(1))/7.74144e6 - (29619*EIL(3))/1.14688e6 + (28203125*EIL(5))/1.0838016e7 - (129296251*EIL(7))/1.769472e7 + (164015*EILN(1))/1.2386304e7 + (345249*EILN(3))/286720. - (297109375*EILN(5))/6.193152e6 + (266004389*EILN(7))/1.10592e6 - (7877549943*EILN(9))/3.211264e7 + e*((229*EIL(2))/6912. - (1193*EIL(4))/4410. + (44469*EIL(6))/7840. - (241664*EIL(8))/19845. + (1777*EILN(2))/80640. + (11539*EILN(4))/1890. - (4014603*EILN(6))/31360. + (217088*EILN(8))/441. - (224609375*EILN(10))/508032. + e*((3697843*EIL(1))/2.97271296e8 + (1107*EIL(3))/17920. - (1444109375*EIL(5))/1.387266048e9 + (114472477*EIL(7))/9.95328e6 - (25440612111*EIL(9))/1.2845056e9 + (240773*EILN(1))/1.93536e7 - (2197233*EILN(3))/1.835008e7 + (95609375*EILN(5))/4.064256e6 - (44020843979*EILN(7))/1.4155776e8 + (387420489*EILN(9))/401408. - (34522712143931*EILN(11))/4.45906944e10)))))))))))
        +I*((3*lgxOx0p*EILN(1))/14. + EILN(1)*(-0.03571428571428571 - CST_LN2/7.) + e*((5*EILN(2))/21. - (10*lgxOx0p*EILN(2))/7. + e*(lgxOx0p*((11*EILN(1))/7. - (1053*EILN(3))/112.) + EILN(1)*(-0.2619047619047619 - (22*CST_LN2)/21.) + EIL(1)*(-0.05505952380952381 + (37*lgxOx0p)/112. - (37*CST_LN2)/168.) + EILN(3)*(1.5669642857142858 + (351*CST_LN2)/56. - (351*CST_LN3)/56.) + e*((-19*EIL(2))/84. - (47*EILN(2))/28. + lgxOx0p*((19*EIL(2))/14. + (141*EILN(2))/14. - 32*EILN(4)) + EILN(4)*(5.333333333333333 - (64*CST_LN2)/3.) + e*(lgxOx0p*((891*EIL(3))/256. - (227*EILN(1))/896. + (1215*EILN(3))/32. - (453125*EILN(5))/5376.) + EIL(1)*(0.02058531746031746 - (83*lgxOx0p)/672. + (83*CST_LN2)/1008.) + EILN(1)*(0.042224702380952384 + (227*CST_LN2)/1344.) + EIL(3)*(-0.580078125 - (297*CST_LN2)/128. + (297*CST_LN3)/128.) + EILN(3)*(-6.328125 - (405*CST_LN2)/16. + (405*CST_LN3)/16.) + EILN(5)*(14.047774057539682 + (453125*CST_LN2)/8064. - (453125*CST_LN5)/8064.) + e*((109*EIL(2))/1008. + (107*EILN(2))/144. + lgxOx0p*((-109*EIL(2))/168. + (776*EIL(4))/105. - (107*EILN(2))/24. + (2360*EILN(4))/21. - (26973*EILN(6))/140.) + EIL(4)*(-1.2317460317460318 + (1552*CST_LN2)/315.) + EILN(4)*(-18.73015873015873 + (4720*CST_LN2)/63.) + EILN(6)*(32.11071428571429 - (8991*CST_LN3)/70.) + e*(lgxOx0p*((-4131*EIL(3))/1792. + (203125*EIL(5))/14336. + (93*EILN(1))/896. - (358425*EILN(3))/14336. + (515625*EILN(5))/1792. - (823543*EILN(7))/2048.) + EILN(1)*(-0.017299107142857144 - (31*CST_LN2)/448.) + EIL(1)*(0.013706752232142858 - (1179*lgxOx0p)/14336. + (393*CST_LN2)/7168.) + EILN(3)*(4.166957310267857 + (119475*CST_LN2)/7168. - (119475*CST_LN3)/7168.) + EIL(3)*(0.3842075892857143 + (1377*CST_LN2)/896. - (1377*CST_LN3)/896.) + EIL(5)*(-2.361479259672619 - (203125*CST_LN2)/21504. + (203125*CST_LN5)/21504.) + EILN(5)*(-47.95619419642857 - (171875*CST_LN2)/896. + (171875*CST_LN5)/896.) + EILN(7)*(67.02010091145833 + (823543*CST_LN2)/3072. - (823543*CST_LN7)/3072.) + e*((113*EIL(2))/5040. - (479*EILN(2))/3024. + lgxOx0p*((-113*EIL(2))/840. - (2038*EIL(4))/315. + (99873*EIL(6))/3920. + (479*EILN(2))/504. - (9902*EILN(4))/105. + (373977*EILN(6))/560. - (1736704*EILN(8))/2205.) + EILN(8)*(131.27014361300076 - (6946816*CST_LN2)/6615.) + EILN(4)*(15.717460317460317 - (19804*CST_LN2)/315.) + EIL(4)*(1.0783068783068783 - (4076*CST_LN2)/945.) + EIL(6)*(-4.246301020408163 + (33291*CST_LN3)/1960.) + EILN(6)*(-111.30267857142857 + (124659*CST_LN3)/280.) + e*(lgxOx0p*((88857*EIL(3))/573440. - (28203125*EIL(5))/1.806336e6 + (129296251*EIL(7))/2.94912e6 + (164015*EILN(1))/2.064384e6 + (1035747*EILN(3))/143360. - (297109375*EILN(5))/1.032192e6 + (266004389*EILN(7))/184320. - (23632649829*EILN(9))/1.605632e7) + EILN(1)*(-0.013241641735904431 - (164015*CST_LN2)/3.096576e6) + EIL(1)*(0.013336924396494709 - (103247*lgxOx0p)/1.29024e6 + (103247*CST_LN2)/1.93536e6) + EILN(9)*(245.309944713359 + (7877549943*CST_LN2)/8.02816e6 - (7877549943*CST_LN3)/4.01408e6) + EIL(3)*(-0.025825718470982143 - (29619*CST_LN2)/286720. + (29619*CST_LN3)/286720.) + EILN(3)*(-1.2041329520089286 - (345249*CST_LN2)/71680. + (345249*CST_LN3)/71680.) + EILN(5)*(47.97385483191758 + (297109375*CST_LN2)/1.548288e6 - (297109375*CST_LN5)/1.548288e6) + EIL(5)*(2.602240576134968 + (28203125*CST_LN2)/2.709504e6 - (28203125*CST_LN5)/2.709504e6) + EIL(7)*(-7.307052668818721 - (129296251*CST_LN2)/4.42368e6 + (129296251*CST_LN7)/4.42368e6) + EILN(7)*(-240.5276954933449 - (266004389*CST_LN2)/276480. + (266004389*CST_LN7)/276480.) + e*((229*EIL(2))/6912. - (1777*EILN(2))/80640. + lgxOx0p*((-229*EIL(2))/1152. + (1193*EIL(4))/735. - (133407*EIL(6))/3920. + (483328*EIL(8))/6615. + (1777*EILN(2))/13440. + (11539*EILN(4))/315. - (12043809*EILN(6))/15680. + (434176*EILN(8))/147. - (224609375*EILN(10))/84672.) + EIL(4)*(-0.27052154195011335 + (2386*CST_LN2)/2205.) + EILN(4)*(-6.105291005291005 + (23078*CST_LN2)/945.) + EIL(8)*(-12.177576215671454 + (1933312*CST_LN2)/19845.) + EILN(8)*(-492.26303854875283 + (1736704*CST_LN2)/441.) + EILN(6)*(128.01667729591836 - (4014603*CST_LN3)/7840.) + EIL(6)*(5.672066326530612 - (44469*CST_LN3)/1960.) + EILN(10)*(442.11658911249685 - (224609375*CST_LN5)/127008.) + e*(lgxOx0p*((-3321*EIL(3))/8960. + (1444109375*EIL(5))/2.31211008e8 - (114472477*EIL(7))/1.65888e6 + (76321836333*EIL(9))/6.422528e8 + (240773*EILN(1))/3.2256e6 - (6591699*EILN(3))/9.17504e6 + (95609375*EILN(5))/677376. - (44020843979*EILN(7))/2.359296e7 + (1162261467*EILN(9))/200704. - (34522712143931*EILN(11))/7.4317824e9) + EILN(1)*(-0.012440734540343916 - (240773*CST_LN2)/4.8384e6) + EIL(1)*(0.012439287108298542 - (3697843*lgxOx0p)/4.9545216e7 + (3697843*CST_LN2)/7.4317824e7) + EILN(3)*(0.11973969595772879 + (2197233*CST_LN2)/4.58752e6 - (2197233*CST_LN3)/4.58752e6) + EIL(3)*(0.061774553571428574 + (1107*CST_LN2)/4480. - (1107*CST_LN3)/4480.) + EIL(9)*(-19.805761929726113 - (25440612111*CST_LN2)/3.211264e8 + (25440612111*CST_LN3)/1.605632e8) + EILN(9)*(-965.1538808394452 - (387420489*CST_LN2)/100352. + (387420489*CST_LN3)/50176.) + EIL(5)*(-1.040975072576706 - (1444109375*CST_LN2)/3.46816512e8 + (1444109375*CST_LN5)/3.46816512e8) + EILN(5)*(-23.52444703286407 - (95609375*CST_LN2)/1.016064e6 + (95609375*CST_LN5)/1.016064e6) + EILN(7)*(310.97443177258526 + (44020843979*CST_LN2)/3.538944e7 - (44020843979*CST_LN7)/3.538944e7) + EIL(7)*(11.500980279867541 + (114472477*CST_LN2)/2.48832e6 - (114472477*CST_LN7)/2.48832e6) + EILN(11)*(774.2133781153002 + (34522712143931*CST_LN2)/1.11476736e10 - (34522712143931*CST_LN11)/1.11476736e10)))))))))))))
        +(CST_2PI*((-17*EILN(1))/168. + e*(-0.047619047619047616*EILN(2) + e*((569*EIL(1))/1344. - (275*EILN(1))/336. + (351*EILN(3))/448. + e*((113*EIL(2))/63. - (64*EILN(2))/21. + (26*EILN(4))/9. + e*((-185*EIL(1))/4032. + (4671*EIL(3))/1024. - (3175*EILN(1))/10752. - (459*EILN(3))/64. + (440875*EILN(5))/64512. + e*((-949*EIL(2))/2016. + (12067*EIL(4))/1260. - (115*EILN(2))/144. - (3515*EILN(4))/252. + (14769*EILN(6))/1120. + e*((-4237*EIL(1))/516096. - (156843*EIL(3))/71680. + (28131325*EIL(5))/1.548288e6 - (42239*EILN(1))/193536. - (75573*EILN(3))/57344. - (3089725*EILN(5))/129024. + (24696343*EILN(7))/1.10592e6 + e*((737*EIL(2))/8064. - (10331*EIL(4))/1512. + (1015173*EIL(6))/31360. - (13459*EILN(2))/24192. - (5669*EILN(4))/2520. - (32967*EILN(6))/896. + (45097*EILN(8))/1323. + e*((-7673*EIL(1))/483840. + (267273*EIL(3))/458752. - (11742925*EIL(5))/677376. + (1955571233*EIL(7))/3.538944e7 - (4612277*EILN(1))/2.4772608e7 - (9531*EILN(3))/17920. - (68295175*EILN(5))/1.2386304e7 - (685853*EILN(7))/13824. + (3018279429*EILN(9))/6.422528e7 + e*((109*EIL(2))/41472. + (63709*EIL(4))/26460. - (1946457*EIL(6))/50176. + (87104537*EIL(8))/952560. - (265367*EILN(2))/483840. + (4951*EILN(4))/4536. - (2115369*EILN(6))/125440. - (1143077*EILN(8))/21168. + (695570875*EILN(10))/1.2192768e7 + e*((-54670843*EIL(1))/2.97271296e9 - (18657*EIL(3))/917504. + (21721444775*EIL(5))/2.774532096e9 - (50850403513*EIL(7))/6.3700992e8 + (75937604247*EIL(9))/5.1380224e8 - (40569527*EILN(1))/2.4772608e8 - (39769839*EILN(3))/3.670016e7 + (4101773675*EILN(5))/5.20224768e8 - (2872440221*EILN(7))/5.6623104e7 - (3904716213*EILN(9))/1.2845056e8 + (990007698647*EILN(11))/1.783627776e10)))))))))))
        +I*((-17*lgxOx0p*EILN(1))/28. + EILN(1)*(0.10119047619047619 + (17*CST_LN2)/42.) + e*((-83*EILN(2))/21. - (2*lgxOx0p*EILN(2))/7. + e*(lgxOx0p*((-275*EILN(1))/56. + (1053*EILN(3))/224.) + EIL(1)*(0.9233630952380952 - (569*lgxOx0p)/224. + (569*CST_LN2)/336.) + EILN(1)*(0.8184523809523809 + (275*CST_LN2)/84.) + EILN(3)*(-14.283482142857142 - (351*CST_LN2)/112. + (351*CST_LN3)/112.) + e*((176*EIL(2))/63. + (85*EILN(2))/21. + lgxOx0p*((-226*EIL(2))/21. - (128*EILN(2))/7. + (52*EILN(4))/3.) + EILN(4)*(-34.888888888888886 + (104*CST_LN2)/9.) + e*(lgxOx0p*((-14013*EIL(3))/512. - (3175*EILN(1))/1792. - (1377*EILN(3))/32. + (440875*EILN(5))/10752.) + EIL(1)*(0.3707837301587302 + (185*lgxOx0p)/672. - (185*CST_LN2)/1008.) + EILN(1)*(0.29529389880952384 + (3175*CST_LN2)/2688.) + EILN(3)*(13.921875 + (459*CST_LN2)/16. - (459*CST_LN3)/16.) + EIL(3)*(6.2490234375 + (4671*CST_LN2)/256. - (4671*CST_LN3)/256.) + EILN(5)*(-71.93816654265873 - (440875*CST_LN2)/16128. + (440875*CST_LN5)/16128.) + e*((563*EIL(2))/2016. + (55*EILN(2))/144. + lgxOx0p*((949*EIL(2))/336. - (12067*EIL(4))/210. - (115*EILN(2))/24. - (3515*EILN(4))/42. + (44307*EILN(6))/560.) + EILN(4)*(37.948412698412696 - (3515*CST_LN2)/63.) + EIL(4)*(12.243650793650794 - (12067*CST_LN2)/315.) + EILN(6)*(-134.68660714285716 + (14769*CST_LN3)/280.) + e*(lgxOx0p*((470529*EIL(3))/35840. - (28131325*EIL(5))/258048. - (42239*EILN(1))/32256. - (226719*EILN(3))/28672. - (3089725*EILN(5))/21504. + (24696343*EILN(7))/184320.) + EIL(1)*(0.36027987041170634 + (4237*lgxOx0p)/86016. - (4237*CST_LN2)/129024.) + EILN(1)*(0.21824880125661375 + (42239*CST_LN2)/48384.) + EILN(3)*(-0.8969552176339286 + (75573*CST_LN2)/14336. - (75573*CST_LN3)/14336.) + EIL(3)*(-1.1755998883928571 - (156843*CST_LN2)/17920. + (156843*CST_LN3)/17920.) + EILN(5)*(89.05106801835318 + (3089725*CST_LN2)/32256. - (3089725*CST_LN5)/32256.) + EIL(5)*(22.23832064835483 + (28131325*CST_LN2)/387072. - (28131325*CST_LN5)/387072.) + EILN(7)*(-236.79536223234953 - (24696343*CST_LN2)/276480. + (24696343*CST_LN7)/276480.) + e*((31909*EIL(2))/40320. + (8083*EILN(2))/24192. + lgxOx0p*((-737*EIL(2))/1344. + (10331*EIL(4))/252. - (3045519*EIL(6))/15680. - (13459*EILN(2))/4032. - (5669*EILN(4))/420. - (98901*EILN(6))/448. + (90194*EILN(8))/441.) + EILN(4)*(-6.950396825396825 - (5669*CST_LN2)/630.) + EIL(4)*(-5.7215608465608465 + (10331*CST_LN2)/378.) + EILN(8)*(-398.1758125472411 + (360776*CST_LN2)/1323.) + EILN(6)*(188.66852678571428 - (32967*CST_LN3)/224.) + EIL(6)*(38.44658801020408 - (1015173*CST_LN3)/7840.) + e*(lgxOx0p*((-801819*EIL(3))/229376. + (11742925*EIL(5))/112896. - (1955571233*EIL(7))/5.89824e6 - (4612277*EILN(1))/4.128768e6 - (28593*EILN(3))/8960. - (68295175*EILN(5))/2.064384e6 - (685853*EILN(7))/2304. + (9054838287*EILN(9))/3.211264e7) + EIL(1)*(0.3154348544973545 + (7673*lgxOx0p)/80640. - (7673*CST_LN2)/120960.) + EILN(1)*(0.18618455513444526 + (4612277*CST_LN2)/6.193152e6) + EIL(3)*(1.621476091657366 + (267273*CST_LN2)/114688. - (267273*CST_LN3)/114688.) + EILN(3)*(0.22600446428571427 + (9531*CST_LN2)/4480. - (9531*CST_LN3)/4480.) + EILN(9)*(-647.5353855210908 - (3018279429*CST_LN2)/1.605632e7 + (3018279429*CST_LN3)/8.02816e6) + EILN(5)*(-25.342896880296173 + (68295175*CST_LN2)/3.096576e6 - (68295175*CST_LN5)/3.096576e6) + EIL(5)*(-16.463972741874528 - (11742925*CST_LN2)/169344. + (11742925*CST_LN5)/169344.) + EIL(7)*(64.19462260493526 + (1955571233*CST_LN2)/8.84736e6 - (1955571233*CST_LN7)/8.84736e6) + EILN(7)*(371.309693287037 + (685853*CST_LN2)/3456. - (685853*CST_LN7)/3456.) + e*((133529*EIL(2))/207360. + (161711*EILN(2))/483840. + lgxOx0p*((-109*EIL(2))/6912. - (63709*EIL(4))/4410. + (5839371*EIL(6))/25088. - (87104537*EIL(8))/158760. - (265367*EILN(2))/80640. + (4951*EILN(4))/756. - (6346107*EILN(6))/62720. - (1143077*EILN(8))/3528. + (695570875*EILN(10))/2.032128e6) + EIL(8)*(104.4457472495171 - (87104537*CST_LN2)/119070.) + EILN(8)*(691.1557917611489 - (1143077*CST_LN2)/2646.) + EIL(4)*(3.820445956160242 - (63709*CST_LN2)/6615.) + EILN(4)*(-0.31371252204585537 + (4951*CST_LN2)/1134.) + EILN(6)*(-71.00694355867347 - (2115369*CST_LN3)/31360.) + EIL(6)*(-38.79259008290816 + (1946457*CST_LN3)/12544.) + EILN(10)*(-1025.8598273173081 + (695570875*CST_LN5)/3.048192e6) + e*(lgxOx0p*((55971*EIL(3))/458752. - (21721444775*EIL(5))/4.62422016e8 + (50850403513*EIL(7))/1.0616832e8 - (227812812741*EIL(9))/2.5690112e8 - (40569527*EILN(1))/4.128768e7 - (119309517*EILN(3))/1.835008e7 + (4101773675*EILN(5))/8.6704128e7 - (2872440221*EILN(7))/9.437184e6 - (11714148639*EILN(9))/6.422528e7 + (990007698647*EILN(11))/2.97271296e9) + EIL(1)*(0.2844368441815519 + (54670843*lgxOx0p)/4.9545216e8 - (54670843*CST_LN2)/7.4317824e8) + EILN(1)*(0.16376768647047577 + (40569527*CST_LN2)/6.193152e7) + EIL(9)*(166.56227452609005 + (75937604247*CST_LN2)/1.2845056e8 - (75937604247*CST_LN3)/6.422528e7) + EILN(9)*(1231.4789774135668 + (3904716213*CST_LN2)/3.211264e7 - (3904716213*CST_LN3)/1.605632e7) + EILN(3)*(0.5881037848336356 + (39769839*CST_LN2)/9.17504e6 - (39769839*CST_LN3)/9.17504e6) + EIL(3)*(0.936417715890067 - (18657*CST_LN2)/229376. + (18657*CST_LN3)/229376.) + EIL(5)*(9.724863090933225 + (21721444775*CST_LN2)/6.93633024e8 - (21721444775*CST_LN5)/6.93633024e8) + EILN(5)*(-1.264402841734748 - (4101773675*CST_LN2)/1.30056192e8 + (4101773675*CST_LN5)/1.30056192e8) + EILN(7)*(-171.27496456216883 + (2872440221*CST_LN2)/1.4155776e7 - (2872440221*CST_LN7)/1.4155776e7) + EIL(7)*(-81.81248662658189 - (50850403513*CST_LN2)/1.5925248e8 + (50850403513*CST_LN7)/1.5925248e8) + EILN(11)*(-1591.1351269124327 - (990007698647*CST_LN2)/4.45906944e9 + (990007698647*CST_LN11)/4.45906944e9)))))))))))));


#undef EIL
#undef EILN

    return dm*v*(I*ImPN1 + v2*(I*ImPN3 + v*(I*ImPN4_tail + v*(I*ImPN5 + v*(PN6 + I*ImPN6_tail)))));
}


COMPLEX16 calculate_h21FourierMode_zeroe(REAL8 eta, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 dm = sqrt(1. - 4.*eta);
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];

    REAL8 ImPN1 = 2./3.;
    REAL8 ImPN3 = (-17. + 20.*eta)/42.;
    REAL8 ImPN5 = (-172. - 2036.*eta + 237.*eta2)/756.;
    REAL8 PN6 = (25.*eta)/3.;

    REAL8 lnxOx0p = log(v2/CONST_X0P);
    COMPLEX16 lgpref = 1. - 6.*lnxOx0p + 2.*I*CST_PI + 4.*CST_LN2;
    COMPLEX16 PN4_tail = lgpref / 6.;
    COMPLEX16 PN6_tail = (-17. + 6.*eta)*lgpref/168.;
    return dm*v*(I*ImPN1 + 
            v2*(I*ImPN3 + 
                v*(PN4_tail + 
                    v*(I*ImPN5 + 
                        v*(PN6 + PN6_tail)))));
}



#define SQ2O3 0.816496580927726032732428024901963797321982494
COMPLEX16 calculate_h20FourierModep(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe)
{
    COMPLEX16 ret;
    INT pm = p;
    if (pm==0) return 0.0;
    INT pm2 = pm*pm;
    INT pm3 = pm2*pm;
    INT pm4 = pm2*pm2;
    INT pm5 = pm4*pm;
    REAL8 *ePow = dc->ePow;
    REAL8 *sq1m2Pow = dc->sq1mesqPow;
    REAL8 m1p3eta = -1. + 3.*eta;
    REAL8 y = sq1m2Pow[1];
    REAL8 yp = y+1.0, myp = -y+1.0;
    REAL8 yp2 = yp*yp, myp2 = myp*myp;
    REAL8 yp3 = yp2*yp, myp3 = myp2*myp;
    REAL8 yp4 = yp3*yp, myp4 = myp3*myp;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
#define EVALEJ(pp, qq, aa, bb) evaluate_Jpqab((pp), (qq), (aa), (bb), pe)
// #define EVALEJ(p, q, a, b) 1.0
#define SQ1M2(n) sq1m2Pow[n]
#define EPOW(n) ePow[n]
    REAL8 PN0, PN2, PN4, ImPN5, PN6;
    PN0 = EVALEJ(pm,pm,0,0);
    PN0 *= SQ2O3;

    PN2 = -((-26 + eta)*EVALEJ(pm,pm,1,0)*SQ1M2(2)) + 
        EVALEJ(pm,pm,0,0)*(-14 - ((-57 + 17*eta)*SQ1M2(2))/6.) + 
        (-26 + eta)*EVALEJ(pm,pm,2,0)*SQ1M2(4);    
    PN2 *= SQ2O3 / (sq1m2Pow[2] * 7);

    PN4 = (-21*(-5 + 2*eta)*pm*EPOW(1)*EVALEJ(-1 + pm,pm,0,1)*SQ1M2(3))/4. - 
        (21*(-5 + 2*eta)*pm*EPOW(1)*EVALEJ(1 + pm,pm,0,1)*SQ1M2(3))/4. - 
        (7*(-15 + eta)*eta*pm*EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*SQ1M2(4))/8. + 
        (7*(-15 + eta)*eta*pm*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*SQ1M2(4))/16. - 
        (7*(-15 + eta)*eta*pm*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*SQ1M2(4))/16. - 
        (EVALEJ(pm,pm,2,0)*(-5604 + 1377*eta - 63*eta2 + 
        (188 + 939*eta + 107*eta2)*SQ1M2(2))*SQ1M2(4))/24. + 
        EVALEJ(pm,pm,1,0)*(-((76 + 45*eta)*SQ1M2(2)) - 
        ((1391 - 1921*eta + 187*eta2)*SQ1M2(4))/36.) + 
        EVALEJ(pm,pm,0,0)*(-87.5 + 63*eta + (21.5 - 19*eta)*SQ1M2(2) - 
        7*(-5 + 2*eta)*SQ1M2(3) + ((155 + 3*eta2*(-3 + 7*pm) - 7*eta*(37 + 45*pm))*
        SQ1M2(4))/24.) + ((-8645 + 8494*eta + 1454*eta2)*EVALEJ(pm,pm,3,0)*
        SQ1M2(6))/72. - ((-217 + 716*eta + 316*eta2)*EVALEJ(pm,pm,4,0)*SQ1M2(8))/24.;   
    PN4 *= SQ2O3/(7.*sq1m2Pow[4]);

    ImPN5 = (122*eta*EVALEJ(pm,pm,0,0))/5. - (73*eta*EVALEJ(pm,pm,1,0))/3. + 
        29*eta*EVALEJ(pm,pm,2,0)*SQ1M2(2);
    ImPN5 *= SQ2O3 * pm / 7.;

    PN6 = (-5*pm*EPOW(1)*EVALEJ(-1 + pm,pm,0,1)*
        (-7*(2880 + 2016*eta2 + eta*(-12512 + 123*CST_PISQ)) + 
        144*(45 - 69*eta + 40*eta2)*SQ1M2(2))*SQ1M2(3))/128. - 
        (5*pm*EPOW(1)*EVALEJ(1 + pm,pm,0,1)*
        (-7*(2880 + 2016*eta2 + eta*(-12512 + 123*CST_PISQ)) + 
        144*(45 - 69*eta + 40*eta2)*SQ1M2(2))*SQ1M2(3))/128. - 
        (pm*EPOW(1)*EVALEJ(-1 + pm,pm,0,0)*
        (-268800 + eta*(208512 + 4305*CST_PISQ - 18220*SQ1M2(2)) + 
        1480*eta3*SQ1M2(2) - 20*eta2*(-84 + 269*SQ1M2(2)))*SQ1M2(4))/64. + 
        (3*pm*EPOW(1)*EVALEJ(-1 + pm,pm,1,0)*
        (-89600 + 400*eta3*SQ1M2(2) - 20*eta2*(-28 + 701*SQ1M2(2)) + 
        eta*(69504 + 1435*CST_PISQ + 41060*SQ1M2(2)))*SQ1M2(4))/128. - 
        (3*pm*EPOW(1)*EVALEJ(1 + pm,pm,1,0)*
        (-89600 + 400*eta3*SQ1M2(2) - 20*eta2*(-28 + 701*SQ1M2(2)) + 
        eta*(69504 + 1435*CST_PISQ + 41060*SQ1M2(2)))*SQ1M2(4))/128. + 
        (EVALEJ(pm,pm,2,0)*SQ1M2(4)*(990*
        (-43680 + 39720*eta2 + eta*(-293984 + 4305*CST_PISQ)) + 
        2138400*(130 - 57*eta + 2*eta2)*SQ1M2(1) + 
        (-146613232 - 125712300*eta2 - 3250000*eta3 + 
        15*eta*(34834276 + 426195*CST_PISQ))*SQ1M2(2) - 
        237600*(130 - 57*eta + 2*eta2)*SQ1M2(3) + 
        300*(22104 - 145521*eta + 107771*eta2 + 3914*eta3)*SQ1M2(4)))/31680. + 
        (45*(130 - 57*eta + 2*eta2)*pm*EVALEJ(pm,pm,1,1)*SQ1M2(5))/2. - 
        (EVALEJ(pm,pm,3,0)*(-32337248 - 27638880*eta2 - 2516000*eta3 + 
        3*eta*(34241272 + 1278585*CST_PISQ) + 
        285120*(130 - 57*eta + 2*eta2)*SQ1M2(1) + 
        12*(-1464553 + 5483195*eta + 1601000*eta2 + 16470*eta3)*SQ1M2(2))*SQ1M2(6))
        /6336. - (EVALEJ(pm,pm,1,0)*SQ1M2(2)*
        (-132*(-104160 + 34320*eta2 + eta*(-25064 + 4305*CST_PISQ)) + 
        33*(-414560 + 200320*eta2 + 21*eta*(-13856 + 615*CST_PISQ))*SQ1M2(2) + 
        31680*(130 - 57*eta + 2*eta2)*SQ1M2(3) + 
        8*(-73189 + 603425*eta - 259145*eta2 + 14210*eta3)*SQ1M2(4) + 
        330*eta*(-2179 - 265*eta + 85*eta2)*pm2*SQ1M2(6)))/2112. - 
        (45*(130 - 57*eta + 2*eta2)*pm*EVALEJ(pm,pm,2,1)*SQ1M2(7))/2. + 
        (EVALEJ(pm,pm,4,0)*(-21812592 + 10734100*eta2 - 931300*eta3 + 
        175*eta*(629020 + 4059*CST_PISQ) + 
        100*(7371 - 64060*eta + 156071*eta2 + 4921*eta3)*SQ1M2(2))*SQ1M2(8))/3520.\
        + EVALEJ(pm,pm,0,0)*((-35*(13312 + 7056*eta2 + 
        eta*(-35968 + 615*CST_PISQ)))/64. + 630*(-5 + 2*eta)*SQ1M2(1) + 
        (5*(67872 + 40608*eta2 + eta*(-111200 + 861*CST_PISQ))*SQ1M2(2))/64. + 
        (35*(8640 + 2016*eta2 + eta*(-14816 + 123*CST_PISQ))*SQ1M2(3))/96. + 
        ((-80*(-1143 + 57*eta + 400*eta2) + 
        9*(-89600 + 69504*eta + 560*eta2 + 1435*eta*CST_PISQ)*pm)*SQ1M2(4))/
        192. - 5*(-75 - 4*eta + 43*eta2)*SQ1M2(5) + 
        (5*(-54135 + eta3*(2195 + 21978*pm + 2376*pm2) + 
        9*eta2*(1832 - 8877*pm + 68640*pm2) - 
        27*eta*(9725 + 10021*pm + 103510*pm2))*SQ1M2(6))/4752. + 
        (35*eta*(23 - 73*eta + 13*eta2)*pm2*SQ1M2(8))/16.) - 
        ((-990343 + 2378115*eta + 2895480*eta2 + 77300*eta3)*EVALEJ(pm,pm,5,0)*
        SQ1M2(10))/264. + (25*(-361 - 1653*eta + 13602*eta2 + 602*eta3)*
        EVALEJ(pm,pm,6,0)*SQ1M2(12))/88.;
    PN6 *= SQ2O3 / (105. * sq1m2Pow[6]);

    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];

    // tail part
    COMPLEX16 PN3_tail, PN5_tail, PN6_tail;
    COMPLEX16 lnIpvO2 =  pm>0 ? (log(0.5*pm) + 0.5*I*CST_PI) : (log(-0.5*pm) - 0.5*I*CST_PI);
    REAL8 lnxOx0p = log(v2/CONST_X0P);
    COMPLEX16 pref_tail, pref_tail_PN6;
    pref_tail = 3.*lnxOx0p + 2*lnIpvO2;
    PN3_tail = -EVALEJ(pm,pm,2,0) + (5*EVALEJ(pm,pm,3,0))/2. - (3*EVALEJ(pm,pm,4,0)*SQ1M2(2))/2.;
    PN3_tail *= pref_tail * SQ2O3 * I / (1.*pm);

    PN5_tail = pref_tail * (
        EVALEJ(pm,pm,2,0)*(5 + (19*(-3 + 2*eta)*SQ1M2(2))/42.) + 
        (EVALEJ(pm,pm,4,0)*SQ1M2(2)*(134 - 4*eta + (81 + 2*eta)*SQ1M2(2)))/4. + 
        EVALEJ(pm,pm,3,0)*(-12.5 - ((93 + 22*eta)*SQ1M2(2))/12.) + 
        (5*(-26 + eta)*EVALEJ(pm,pm,5,0)*SQ1M2(4))/2. - 
        (15*(-26 + eta)*EVALEJ(pm,pm,6,0)*SQ1M2(6))/14.
    ) + (
        6*EVALEJ(pm,pm,2,0) - 15*EVALEJ(pm,pm,3,0) - 
        (EVALEJ(pm,pm,4,0)*(42*CST_LN2 - 7*(54 + 6*CST_LN2))*SQ1M2(2))/42.
    );
    PN5_tail *= I * SQ2O3 / (sq1m2Pow[2]*pm);

    pref_tail_PN6 = 515063 + 396900*lnxOx0p*lnxOx0p + 179760*lnIpvO2 + 
        176400*lnIpvO2*lnIpvO2 + 2520*lnxOx0p*(107 + 210*lnIpvO2) + 
        29400*CST_PISQ;
    PN6_tail = -0.000011337868480725624*EVALEJ(pm,pm,2,0) + EVALEJ(pm,pm,3,0)/35280. - 
        (EVALEJ(pm,pm,4,0)*SQ1M2(2))/58800.;
    PN6_tail *= pref_tail_PN6 * SQ2O3;

    // print_debug("PN0 = %.16e\n", PN0);
    // print_debug("PN2 = %.16e\n", dc->vPow[2]*PN2);
    // print_debug("PN3_tail = %.16e, %.16e\n", dc->vPow[3]*creal(PN3_tail), dc->vPow[3]*cimag(PN3_tail) );
    // print_debug("PN4 = %.16e\n", dc->vPow[4]*PN4);
    // print_debug("ImPN5 = %.16e\n", dc->vPow[5]*ImPN5);
    // print_debug("PN5_tail = %.16e, %.16e\n", dc->vPow[5]*creal(PN5_tail), dc->vPow[5]*cimag(PN5_tail) );
    // print_debug("PN6 = %.16e\n", dc->vPow[6]*PN6);
    // print_debug("PN6_tail = %.16e, %.16e\n", dc->vPow[6]*creal(PN6_tail), dc->vPow[6]*cimag(PN6_tail) );
    // return v2*PN0;
    return (PN0 + 
        v2 * (PN2 + 
            v*(PN3_tail + 
                v*(PN4 + 
                    v*( I*ImPN5 + PN5_tail +
                        v*(PN6 + PN6_tail) )))));
#undef EVALEJ
#undef SQ1M2
#undef EPOW
}


COMPLEX16 calculate_h20FourierMode_zeroe(REAL8 eta, DynParamsCache *dc)
{
    return 0.0;;
}


COMPLEX16 calculate_h20_pc_e10(REAL8 eta, REAL8 l, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 v = dc->vPow[1], v2 = dc->vPow[2];
    REAL8 e = dc->ePow[1], e2 = dc->ePow[2];
    COMPLEX16 eilvec[50], eilnvec[50];
    eilvec[0] = 1.0;
    eilnvec[0] = 1.0;
    for (size_t p=1; p<50; p++) {
        eilvec[p] = cexp(I*p*l);
        eilnvec[p] = cexp(-I*p*l);
    }
#define EIL(x) eilvec[x]
#define EILN(x) eilnvec[x]
    COMPLEX16 PN0 = (e*(EIL(1)/2. + EILN(1)/2. + e*(EIL(2)/2. + EILN(2)/2. + e*(-0.0625*EIL(1) + (9*EIL(3))/16. - EILN(1)/16. + (9*EILN(3))/16. + e*(-0.16666666666666666*EIL(2) + (2*EIL(4))/3. - EILN(2)/6. + (2*EILN(4))/3. + e*(EIL(1)/384. - (81*EIL(3))/256. + (625*EIL(5))/768. + EILN(1)/384. - (81*EILN(3))/256. + (625*EILN(5))/768. + e*(EIL(2)/48. - (8*EIL(4))/15. + (81*EIL(6))/80. + EILN(2)/48. - (8*EILN(4))/15. + (81*EILN(6))/80. + e*(-0.00005425347222222222*EIL(1) + (729*EIL(3))/10240. - (15625*EIL(5))/18432. + (117649*EIL(7))/92160. - EILN(1)/18432. + (729*EILN(3))/10240. - (15625*EILN(5))/18432. + (117649*EILN(7))/92160. + e*(-0.001388888888888889*EIL(2) + (8*EIL(4))/45. - (729*EIL(6))/560. + (512*EIL(8))/315. - EILN(2)/720. + (8*EILN(4))/45. - (729*EILN(6))/560. + (512*EILN(8))/315. + e*(EIL(1)/1.47456e6 - (729*EIL(3))/81920. + (390625*EIL(5))/1.032192e6 - (5764801*EIL(7))/2.94912e6 + (4782969*EIL(9))/2.29376e6 + EILN(1)/1.47456e6 - (729*EILN(3))/81920. + (390625*EILN(5))/1.032192e6 - (5764801*EILN(7))/2.94912e6 + (4782969*EILN(9))/2.29376e6 + e*(EIL(2)/17280. - (32*EIL(4))/945. + (6561*EIL(6))/8960. - (8192*EIL(8))/2835. + (390625*EIL(10))/145152. + EILN(2)/17280. - (32*EILN(4))/945. + (6561*EILN(6))/8960. - (8192*EILN(8))/2835. + (390625*EILN(10))/145152.)))))))))));

    COMPLEX16 PN2 = eta*(e*((-11*EIL(1))/84. - (11*EILN(1))/84. + e*((-5*EIL(2))/84. - (5*EILN(2))/84. + e*(-0.001488095238095238*EIL(1) + (3*EIL(3))/224. - EILN(1)/672. + (3*EILN(3))/224. + e*(-0.027777777777777776*EIL(2) + EIL(4)/9. - EILN(2)/36. + EILN(4)/9. + e*((13*EIL(1))/16128. - (351*EIL(3))/3584. + (8125*EIL(5))/32256. + (13*EILN(1))/16128. - (351*EILN(3))/3584. + (8125*EILN(5))/32256. + e*((19*EIL(2))/2016. - (76*EIL(4))/315. + (513*EIL(6))/1120. + (19*EILN(2))/2016. - (76*EILN(4))/315. + (513*EILN(6))/1120. + e*((-25*EIL(1))/774144. + (1215*EIL(3))/28672. - (390625*EIL(5))/774144. + (84035*EIL(7))/110592. - (25*EILN(1))/774144. + (1215*EILN(3))/28672. - (390625*EILN(5))/774144. + (84035*EILN(7))/110592. + e*((-31*EIL(2))/30240. + (124*EIL(4))/945. - (7533*EIL(6))/7840. + (7936*EIL(8))/6615. - (31*EILN(2))/30240. + (124*EILN(4))/945. - (7533*EILN(6))/7840. + (7936*EILN(8))/6615. + e*((37*EIL(1))/6.193152e7 - (8991*EIL(3))/1.14688e6 + (14453125*EIL(5))/4.3352064e7 - (30471091*EIL(7))/1.769472e7 + (58989951*EIL(9))/3.211264e7 + (37*EILN(1))/6.193152e7 - (8991*EILN(3))/1.14688e6 + (14453125*EILN(5))/4.3352064e7 - (30471091*EILN(7))/1.769472e7 + (58989951*EILN(9))/3.211264e7 + e*((43*EIL(2))/725760. - (688*EIL(4))/19845. + (94041*EIL(6))/125440. - (176128*EIL(8))/59535. + (16796875*EIL(10))/6.096384e6 + (43*EILN(2))/725760. - (688*EILN(4))/19845. + (94041*EILN(6))/125440. - (176128*EILN(8))/59535. + (16796875*EILN(10))/6.096384e6)))))))))))
        +(e*((-61*EIL(1))/28. - (61*EILN(1))/28. + e*((-113*EIL(2))/28. - (113*EILN(2))/28. + e*((-59*EIL(1))/224. - (1485*EIL(3))/224. - (59*EILN(1))/224. - (1485*EILN(3))/224. + e*((19*EIL(2))/12. - (31*EIL(4))/3. + (19*EILN(2))/12. - (31*EILN(4))/3. + e*((-4973*EIL(1))/5376. + (17757*EIL(3))/3584. - (168125*EIL(5))/10752. - (4973*EILN(1))/5376. + (17757*EILN(3))/3584. - (168125*EILN(5))/10752. + e*((-769*EIL(2))/672. + (1144*EIL(4))/105. - (26001*EIL(6))/1120. - (769*EILN(2))/672. + (1144*EILN(4))/105. - (26001*EILN(6))/1120. + e*((-226763*EIL(1))/258048. - (342477*EIL(3))/143360. + (5408125*EIL(5))/258048. - (6269011*EIL(7))/184320. - (226763*EILN(1))/258048. - (342477*EILN(3))/143360. + (5408125*EILN(5))/258048. - (6269011*EILN(7))/184320. + e*((-1343*EIL(2))/2016. - (1784*EIL(4))/315. + (293949*EIL(6))/7840. - (21760*EIL(8))/441. - (1343*EILN(2))/2016. - (1784*EILN(4))/315. + (293949*EILN(6))/7840. - (21760*EILN(8))/441. + e*((-18169117*EIL(1))/2.064384e7 - (380043*EIL(3))/1.14688e6 - (185348125*EIL(5))/1.4450688e7 + (377770939*EIL(7))/5.89824e6 - (2281476213*EIL(9))/3.211264e7 - (18169117*EILN(1))/2.064384e7 - (380043*EILN(3))/1.14688e6 - (185348125*EILN(5))/1.4450688e7 + (377770939*EILN(7))/5.89824e6 - (2281476213*EILN(9))/3.211264e7 + e*((-171217*EIL(2))/241920. + (4348*EIL(4))/6615. - (3398193*EIL(6))/125440. + (2102272*EIL(8))/19845. - (206640625*EIL(10))/2.032128e6 - (171217*EILN(2))/241920. + (4348*EILN(4))/6615. - (3398193*EILN(6))/125440. + (2102272*EILN(8))/19845. - (206640625*EILN(10))/2.032128e6)))))))))));

    COMPLEX16 PN4 = eta2*(e*((-95*EIL(1))/1008. - (95*EILN(1))/1008. + e*((-265*EIL(2))/336. - (265*EILN(2))/336. + e*((709*EIL(1))/2688. - (2127*EIL(3))/896. + (709*EILN(1))/2688. - (2127*EILN(3))/896. + e*((4091*EIL(2))/3024. - (4091*EIL(4))/756. + (4091*EILN(2))/3024. - (4091*EILN(4))/756. + e*((-743*EIL(1))/21504. + (60183*EIL(3))/14336. - (464375*EIL(5))/43008. - (743*EILN(1))/21504. + (60183*EILN(3))/14336. - (464375*EILN(5))/43008. + e*((-3305*EIL(2))/8064. + (661*EIL(4))/63. - (17847*EIL(6))/896. - (3305*EILN(2))/8064. + (661*EILN(4))/63. - (17847*EILN(6))/896. + e*((13775*EIL(1))/9.289728e6 - (223155*EIL(3))/114688. + (215234375*EIL(5))/9.289728e6 - (46303285*EIL(7))/1.327104e6 + (13775*EILN(1))/9.289728e6 - (223155*EILN(3))/114688. + (215234375*EILN(5))/9.289728e6 - (46303285*EILN(7))/1.327104e6 + e*((6089*EIL(2))/120960. - (6089*EIL(4))/945. + (1479627*EIL(6))/31360. - (389696*EIL(8))/6615. + (6089*EILN(2))/120960. - (6089*EILN(4))/945. + (1479627*EILN(6))/31360. - (389696*EILN(8))/6615. + e*((-2599*EIL(1))/8.257536e7 + (1894671*EIL(3))/4.58752e6 - (1015234375*EIL(5))/5.7802752e7 + (2140388257*EIL(7))/2.359296e7 - (12430936431*EIL(9))/1.2845056e8 - (2599*EILN(1))/8.257536e7 + (1894671*EILN(3))/4.58752e6 - (1015234375*EILN(5))/5.7802752e7 + (2140388257*EILN(7))/2.359296e7 - (12430936431*EILN(9))/1.2845056e8 + e*((-29147*EIL(2))/8.70912e6 + (116588*EIL(4))/59535. - (21248163*EIL(6))/501760. + (29846528*EIL(8))/178605. - (11385546875*EIL(10))/7.3156608e7 - (29147*EILN(2))/8.70912e6 + (116588*EILN(4))/59535. - (21248163*EILN(6))/501760. + (29846528*EILN(8))/178605. - (11385546875*EILN(10))/7.3156608e7)))))))))))
        +eta*(e*((4979*EIL(1))/1008. + (4979*EILN(1))/1008. + e*((3389*EIL(2))/336. + (3389*EILN(2))/336. + e*((8237*EIL(1))/896. + (2205*EIL(3))/128. + (8237*EILN(1))/896. + (2205*EILN(3))/128. + e*((30841*EIL(2))/3024. + (2915*EIL(4))/108. + (30841*EILN(2))/3024. + (2915*EILN(4))/108. + e*((819745*EIL(1))/64512. + (158805*EIL(3))/14336. + (5183425*EIL(5))/129024. + (819745*EILN(1))/64512. + (158805*EILN(3))/14336. + (5183425*EILN(5))/129024. + e*((93595*EIL(2))/8064. + (7057*EIL(4))/630. + (258633*EIL(6))/4480. + (93595*EILN(2))/8064. + (7057*EILN(4))/630. + (258633*EILN(6))/4480. + e*((156575221*EIL(1))/9.289728e6 + (1135935*EIL(3))/114688. + (95199325*EIL(5))/9.289728e6 + (107186569*EIL(7))/1.327104e6 + (156575221*EILN(1))/9.289728e6 + (1135935*EILN(3))/114688. + (95199325*EILN(5))/9.289728e6 + (107186569*EILN(7))/1.327104e6 + e*((158377*EIL(2))/10080. + (14641*EIL(4))/2520. + (64503*EIL(6))/7840. + (1949389*EIL(8))/17640. + (158377*EILN(2))/10080. + (14641*EILN(4))/2520. + (64503*EILN(6))/7840. + (1949389*EILN(8))/17640. + e*((5169784241*EIL(1))/2.4772608e8 + (14882121*EIL(3))/917504. - (549831775*EIL(5))/1.73408256e8 + (83032901*EIL(7))/1.4155776e7 + (19037961603*EIL(9))/1.2845056e8 + (5169784241*EILN(1))/2.4772608e8 + (14882121*EILN(3))/917504. - (549831775*EILN(5))/1.73408256e8 + (83032901*EILN(7))/1.4155776e7 + (19037961603*EILN(9))/1.2845056e8 + e*((163927487*EIL(2))/8.70912e6 + (639601*EIL(4))/34020. - (2697057*EIL(6))/125440. + (1907729*EIL(8))/357210. + (1018866925*EIL(10))/5.225472e6 + (163927487*EILN(2))/8.70912e6 + (639601*EILN(4))/34020. - (2697057*EILN(6))/125440. + (1907729*EILN(8))/357210. + (1018866925*EILN(10))/5.225472e6)))))))))))
        +(e*((-107*EIL(1))/126. - (107*EILN(1))/126. + e*((-2357*EIL(2))/336. - (2357*EILN(2))/336. + e*((-459*EIL(1))/112. - (297*EIL(3))/16. - (459*EILN(1))/112. - (297*EILN(3))/16. + e*((6535*EIL(2))/3024. - (14297*EIL(4))/378. + (6535*EILN(2))/3024. - (14297*EILN(4))/378. + e*((-89*EIL(1))/8. + (3717*EIL(3))/256. - (122875*EIL(5))/1792. - (89*EILN(1))/8. + (3717*EILN(3))/256. - (122875*EILN(5))/1792. + e*((-9109*EIL(2))/1008. + (50333*EIL(4))/1260. - (32589*EIL(6))/280. - (9109*EILN(2))/1008. + (50333*EILN(4))/1260. - (32589*EILN(6))/280. + e*((-2791723*EIL(1))/165888. - (118701*EIL(3))/10240. + (103245875*EIL(5))/1.161216e6 - (157118549*EIL(7))/829440. - (2791723*EILN(1))/165888. - (118701*EILN(3))/10240. + (103245875*EILN(5))/1.161216e6 - (157118549*EILN(7))/829440. + e*((-3100627*EIL(2))/241920. - (341479*EIL(4))/15120. + (11169279*EIL(6))/62720. - (31694851*EIL(8))/105840. - (3100627*EILN(2))/241920. - (341479*EILN(4))/15120. + (11169279*EILN(6))/62720. - (31694851*EILN(8))/105840. + e*((-695932631*EIL(1))/3.096576e7 - (610209*EIL(3))/57344. - (69864875*EIL(5))/1.354752e6 + (590750419*EIL(7))/1.769472e6 - (7441151103*EIL(9))/1.605632e7 - (695932631*EILN(1))/3.096576e7 - (610209*EILN(3))/57344. - (69864875*EILN(5))/1.354752e6 + (590750419*EILN(7))/1.769472e6 - (7441151103*EILN(9))/1.605632e7 + e*((-151518403*EIL(2))/8.70912e6 - (493477*EIL(4))/68040. - (117475839*EIL(6))/1.00352e6 + (426933163*EIL(8))/714420. - (14749997375*EIL(10))/2.0901888e7 - (151518403*EILN(2))/8.70912e6 - (493477*EILN(4))/68040. - (117475839*EILN(6))/1.00352e6 + (426933163*EILN(8))/714420. - (14749997375*EILN(10))/2.0901888e7)))))))))));

    COMPLEX16 ImPN5 = eta*(e*((941*EIL(1))/210. - (941*EILN(1))/210. + e*((1411*EIL(2))/105. - (1411*EILN(2))/105. + e*((-417*EIL(1))/560. + (16859*EIL(3))/560. + (417*EILN(1))/560. - (16859*EILN(3))/560. + e*((-2036*EIL(2))/315. + (37301*EIL(4))/630. + (2036*EILN(2))/315. - (37301*EILN(4))/630. + e*((1663*EIL(1))/5760. - (29903*EIL(3))/1280. + (248815*EIL(5))/2304. - (1663*EILN(1))/5760. + (29903*EILN(3))/1280. - (248815*EILN(5))/2304. + e*((2719*EIL(2))/1680. - (131491*EIL(4))/2100. + (524917*EIL(6))/2800. - (2719*EILN(2))/1680. + (131491*EILN(4))/2100. - (524917*EILN(6))/2800. + e*((371089*EIL(1))/1.93536e6 + (2761719*EIL(3))/358400. - (55623815*EIL(5))/387072. + (434391839*EIL(7))/1.3824e6 - (371089*EILN(1))/1.93536e6 - (2761719*EILN(3))/358400. + (55623815*EILN(5))/387072. - (434391839*EILN(7))/1.3824e6 + e*((19573*EIL(2))/75600. + (1043131*EIL(4))/37800. - (5873083*EIL(6))/19600. + (271525561*EIL(8))/529200. - (19573*EILN(2))/75600. - (1043131*EILN(4))/37800. + (5873083*EILN(6))/19600. - (271525561*EILN(8))/529200. + e*((8891287*EIL(1))/5.16096e7 - (1714269*EIL(3))/2.8672e6 + (582785815*EIL(5))/7.225344e6 - (8629153687*EIL(7))/1.47456e7 + (65905707927*EIL(9))/8.02816e7 - (8891287*EILN(1))/5.16096e7 + (1714269*EILN(3))/2.8672e6 - (582785815*EILN(5))/7.225344e6 + (8629153687*EILN(7))/1.47456e7 - (65905707927*EILN(9))/8.02816e7 + e*((602857*EIL(2))/1.8144e6 - (8854243*EIL(4))/1.5876e6 + (64426977*EIL(6))/313600. - (10379269507*EIL(8))/9.5256e6 + (3939230545*EIL(10))/3.048192e6 - (602857*EILN(2))/1.8144e6 + (8854243*EILN(4))/1.5876e6 - (64426977*EILN(6))/313600. + (10379269507*EILN(8))/9.5256e6 - (3939230545*EILN(10))/3.048192e6)))))))))));

    COMPLEX16 PN6 = eta3*(e*((-20743*EIL(1))/199584. - (20743*EILN(1))/199584. + e*((-9847*EIL(2))/28512. - (9847*EILN(2))/28512. + e*((-334109*EIL(1))/1.596672e6 - (77347*EIL(3))/177408. - (334109*EILN(1))/1.596672e6 - (77347*EILN(3))/177408. + e*((-636203*EIL(2))/598752. + (247*EIL(4))/1944. - (636203*EILN(2))/598752. + (247*EILN(4))/1944. + e*((217079*EIL(1))/5.474304e6 - (1108097*EIL(3))/315392. + (178240625*EIL(5))/7.6640256e7 + (217079*EILN(1))/5.474304e6 - (1108097*EILN(3))/315392. + (178240625*EILN(5))/7.6640256e7 + e*((2431991*EIL(2))/4.790016e6 - (1814807*EIL(4))/187110. + (786167*EIL(6))/98560. + (2431991*EILN(2))/4.790016e6 - (1814807*EILN(4))/187110. + (786167*EILN(6))/98560. + e*((-4055717*EIL(1))/1.839366144e9 + (32798349*EIL(3))/1.261568e7 - (44083578125*EIL(5))/1.839366144e9 + (2424762697*EIL(7))/1.1943936e8 - (4055717*EILN(1))/1.839366144e9 + (32798349*EILN(3))/1.261568e7 - (44083578125*EILN(5))/1.839366144e9 + (2424762697*EILN(7))/1.1943936e8 + e*((-5838563*EIL(2))/7.185024e7 + (5221379*EIL(4))/561330. - (37734651*EIL(6))/689920. + (88084576*EIL(8))/1.964655e6 - (5838563*EILN(2))/7.185024e7 + (5221379*EILN(4))/561330. - (37734651*EILN(6))/689920. + (88084576*EILN(8))/1.964655e6 + e*((769339*EIL(1))/1.337720832e10 - (72461457*EIL(3))/1.0092544e8 + (403368359375*EIL(5))/1.4714929152e10 - (4936310977199*EIL(7))/4.204265472e10 + (36679290183*EIL(9))/4.0370176e8 + (769339*EILN(1))/1.337720832e10 - (72461457*EILN(3))/1.0092544e8 + (403368359375*EILN(5))/1.4714929152e10 - (4936310977199*EILN(7))/4.204265472e10 + (36679290183*EILN(9))/4.0370176e8 + e*((11376047*EIL(2))/1.72440576e9 - (21517726*EIL(4))/5.893965e6 + (788148063*EIL(6))/1.103872e7 - (4244545024*EIL(8))/1.7681895e7 + (228642578125*EIL(10))/1.316818944e9 + (11376047*EILN(2))/1.72440576e9 - (21517726*EILN(4))/5.893965e6 + (788148063*EILN(6))/1.103872e7 - (4244545024*EILN(8))/1.7681895e7 + (228642578125*EILN(10))/1.316818944e9)))))))))))
        +eta2*(e*((-80359*EIL(1))/11088. - (80359*EILN(1))/11088. + e*((-59069*EIL(2))/2772. - (59069*EILN(2))/2772. + e*((-3314083*EIL(1))/88704. - (387357*EIL(3))/9856. - (3314083*EILN(1))/88704. - (387357*EILN(3))/9856. + e*((-574271*EIL(2))/8316. - (5951*EIL(4))/108. - (574271*EILN(2))/8316. - (5951*EILN(4))/108. + e*((-23028965*EIL(1))/304128. - (21035629*EIL(3))/157696. - (234662075*EIL(5))/4.257792e6 - (23028965*EILN(1))/304128. - (21035629*EILN(3))/157696. - (234662075*EILN(5))/4.257792e6 + e*((-22321997*EIL(2))/266112. - (11089033*EIL(4))/41580. - (592677*EIL(6))/49280. - (22321997*EILN(2))/266112. - (11089033*EILN(4))/41580. - (592677*EILN(6))/49280. + e*((-13800164447*EIL(1))/1.02187008e8 - (435494401*EIL(3))/6.30784e6 - (55701643175*EIL(5))/1.02187008e8 + (831135403*EIL(7))/6.63552e6 - (13800164447*EILN(1))/1.02187008e8 - (435494401*EILN(3))/6.30784e6 - (55701643175*EILN(5))/1.02187008e8 + (831135403*EILN(7))/6.63552e6 + e*((-1184511509*EIL(2))/7.98336e6 + (13292857*EIL(4))/498960. - (770767749*EIL(6))/689920. + (1564330753*EIL(8))/3.49272e6 - (1184511509*EILN(2))/7.98336e6 + (13292857*EILN(4))/498960. - (770767749*EILN(6))/689920. + (1564330753*EILN(8))/3.49272e6 + e*((-156112829149*EIL(1))/7.4317824e8 - (9170154761*EIL(3))/5.046272e7 + (1887160455625*EIL(5))/5.722472448e9 - (5281636747303*EIL(7))/2.33570304e9 + (1572640485921*EIL(9))/1.41295616e9 - (156112829149*EILN(1))/7.4317824e8 - (9170154761*EILN(3))/5.046272e7 + (1887160455625*EILN(5))/5.722472448e9 - (5281636747303*EILN(7))/2.33570304e9 + (1572640485921*EILN(9))/1.41295616e9 + e*((-5134793437*EIL(2))/2.395008e7 - (1371006013*EIL(4))/5.23908e6 + (12250602927*EIL(6))/1.103872e7 - (70454537149*EIL(8))/1.571724e7 + (348922018175*EIL(10))/1.46313216e8 - (5134793437*EILN(2))/2.395008e7 - (1371006013*EILN(4))/5.23908e6 + (12250602927*EILN(6))/1.103872e7 - (70454537149*EILN(8))/1.571724e7 + (348922018175*EILN(10))/1.46313216e8)))))))))))
        +eta*(e*((52.29204545454545 - (451*CST_PISQ)/192.)*EIL(1) + (52.29204545454545 - (451*CST_PISQ)/192.)*EILN(1) + e*((116.28940295815296 - (697*CST_PISQ)/192.)*EIL(2) + (116.28940295815296 - (697*CST_PISQ)/192.)*EILN(2) + e*((200.45017812049062 - (18409*CST_PISQ)/3072.)*EIL(1) + (200.87972808441558 - (5289*CST_PISQ)/1024.)*EIL(3) + (200.45017812049062 - (18409*CST_PISQ)/3072.)*EILN(1) + (200.87972808441558 - (5289*CST_PISQ)/1024.)*EILN(3) + e*((297.6412067099567 - (287*CST_PISQ)/36.)*EIL(2) + (295.1253968253968 - (7913*CST_PISQ)/1152.)*EIL(4) + (297.6412067099567 - (287*CST_PISQ)/36.)*EILN(2) + (295.1253968253968 - (7913*CST_PISQ)/1152.)*EILN(4) + e*((414.9309319947992 - (199301*CST_PISQ)/18432.)*EIL(1) + (484.3738598315747 - (89667*CST_PISQ)/8192.)*EIL(3) + (378.22438836843133 - (632425*CST_PISQ)/73728.)*EIL(5) + (414.9309319947992 - (199301*CST_PISQ)/18432.)*EILN(1) + (484.3738598315747 - (89667*CST_PISQ)/8192.)*EILN(3) + (378.22438836843133 - (632425*CST_PISQ)/73728.)*EILN(5) + e*((439.13330552549303 - (452435*CST_PISQ)/36864.)*EIL(2) + (815.9063900913901 - (22099*CST_PISQ)/1440.)*EIL(4) + (408.96030438311686 - (204057*CST_PISQ)/20480.)*EIL(6) + (439.13330552549303 - (452435*CST_PISQ)/36864.)*EILN(2) + (815.9063900913901 - (22099*CST_PISQ)/1440.)*EILN(4) + (408.96030438311686 - (204057*CST_PISQ)/20480.)*EILN(6) + e*((725.6506822843859 - (60521617*CST_PISQ)/3.538944e6)*EIL(1) + (462.3532783647017 - (9268173*CST_PISQ)/655360.)*EIL(3) + (1415.791510159491 - (78466825*CST_PISQ)/3.538944e6)*EIL(5) + (311.07054425274885 - (184488479*CST_PISQ)/1.769472e7)*EIL(7) + (725.6506822843859 - (60521617*CST_PISQ)/3.538944e6)*EILN(1) + (462.3532783647017 - (9268173*CST_PISQ)/655360.)*EILN(3) + (1415.791510159491 - (78466825*CST_PISQ)/3.538944e6)*EILN(5) + (311.07054425274885 - (184488479*CST_PISQ)/1.769472e7)*EILN(7) + e*((721.73402643749 - (4012219*CST_PISQ)/221184.)*EIL(2) + (376.3265129469296 - (1073831*CST_PISQ)/69120.)*EIL(4) + (2516.2526607432746 - (9527949*CST_PISQ)/286720.)*EIL(6) + (-49.47952859089764 - (863911*CST_PISQ)/96768.)*EIL(8) + (721.73402643749 - (4012219*CST_PISQ)/221184.)*EILN(2) + (376.3265129469296 - (1073831*CST_PISQ)/69120.)*EILN(4) + (2516.2526607432746 - (9527949*CST_PISQ)/286720.)*EILN(6) + (-49.47952859089764 - (863911*CST_PISQ)/96768.)*EILN(8) + e*((1121.1707269225212 - (3512140237*CST_PISQ)/1.4155776e8)*EIL(1) + (799.7854805012493 - (13286091*CST_PISQ)/655360.)*EIL(3) + (-8.785092961447143 - (751580225*CST_PISQ)/4.9545216e7)*EIL(5) + (4535.186929810991 - (14627000633*CST_PISQ)/2.8311552e8)*EIL(7) + (-899.8478939418758 - (274716441*CST_PISQ)/7.340032e7)*EIL(9) + (1121.1707269225212 - (3512140237*CST_PISQ)/1.4155776e8)*EILN(1) + (799.7854805012493 - (13286091*CST_PISQ)/655360.)*EILN(3) + (-8.785092961447143 - (751580225*CST_PISQ)/4.9545216e7)*EILN(5) + (4535.186929810991 - (14627000633*CST_PISQ)/2.8311552e8)*EILN(7) + (-899.8478939418758 - (274716441*CST_PISQ)/7.340032e7)*EILN(9) + e*((1058.5699417882947 - (41628079*CST_PISQ)/1.65888e6)*EIL(2) + (953.5550198126389 - (33435541*CST_PISQ)/1.45152e6)*EIL(4) + (-1072.1630275339894 - (48472209*CST_PISQ)/4.58752e6)*EIL(6) + (8204.791972101972 - (360628087*CST_PISQ)/4.35456e6)*EIL(8) + (-2614.863656011088 + (1765991975*CST_PISQ)/2.22953472e8)*EIL(10) + (1058.5699417882947 - (41628079*CST_PISQ)/1.65888e6)*EILN(2) + (953.5550198126389 - (33435541*CST_PISQ)/1.45152e6)*EILN(4) + (-1072.1630275339894 - (48472209*CST_PISQ)/4.58752e6)*EILN(6) + (8204.791972101972 - (360628087*CST_PISQ)/4.35456e6)*EILN(8) + (-2614.863656011088 + (1765991975*CST_PISQ)/2.22953472e8)*EILN(10))))))))))))
        +(e*((2808679*EIL(1))/415800. + (2808679*EILN(1))/415800. + e*((35246719*EIL(2))/1.6632e6 + (35246719*EILN(2))/1.6632e6 + e*((-56035127*EIL(1))/1.6632e6 + (12697933*EIL(3))/184800. - (56035127*EILN(1))/1.6632e6 + (12697933*EILN(3))/184800. + e*((-257175379*EIL(2))/4.9896e6 + (10412803*EIL(4))/56700. - (257175379*EILN(2))/4.9896e6 + (10412803*EILN(4))/56700. + e*((-8162674601*EIL(1))/7.98336e7 - (834407369*EIL(3))/5.9136e6 + (2693240413*EIL(5))/6.386688e6 - (8162674601*EILN(1))/7.98336e7 - (834407369*EILN(3))/5.9136e6 + (2693240413*EILN(5))/6.386688e6 + e*((-1355880343*EIL(2))/1.99584e7 - (670725839*EIL(4))/1.782e6 + (3228255017*EIL(6))/3.696e6 - (1355880343*EILN(2))/1.99584e7 - (670725839*EILN(4))/1.782e6 + (3228255017*EILN(6))/3.696e6 + e*((-202770349471*EIL(1))/9.580032e8 - (1489582841*EIL(3))/5.9136e7 - (173988549563*EIL(5))/1.9160064e8 + (104759435123*EIL(7))/6.2208e7 - (202770349471*EILN(1))/9.580032e8 - (1489582841*EILN(3))/5.9136e7 - (173988549563*EILN(5))/1.9160064e8 + (104759435123*EILN(7))/6.2208e7 + e*((-187709086793*EIL(2))/1.197504e9 + (7812898921*EIL(4))/7.4844e7 - (69443564979*EIL(6))/3.4496e7 + (201873521957*EIL(8))/6.54885e7 - (187709086793*EILN(2))/1.197504e9 + (7812898921*EILN(4))/7.4844e7 - (69443564979*EILN(6))/3.4496e7 + (201873521957*EILN(8))/6.54885e7 + e*((-9926539925371*EIL(1))/2.7869184e10 - (90030904567*EIL(3))/6.30784e8 + (4104421372873*EIL(5))/8.583708672e9 - (2563121810208331*EIL(7))/6.13122048e11 + (95844230290033*EIL(9))/1.7661952e10 - (9926539925371*EILN(1))/2.7869184e10 - (90030904567*EILN(3))/6.30784e8 + (4104421372873*EILN(5))/8.583708672e9 - (2563121810208331*EILN(7))/6.13122048e11 + (95844230290033*EILN(9))/1.7661952e10 + e*((-548201681663*EIL(2))/2.052864e9 - (76634463569*EIL(4))/4.49064e8 + (787326257491*EIL(6))/5.51936e8 - (311227493320943*EIL(8))/3.7721376e10 + (1162298638603*EIL(10))/1.25411328e8 - (548201681663*EILN(2))/2.052864e9 - (76634463569*EILN(4))/4.49064e8 + (787326257491*EILN(6))/5.51936e8 - (311227493320943*EILN(8))/3.7721376e10 + (1162298638603*EILN(10))/1.25411328e8)))))))))));

    REAL8 lgxOx0p = log(v2/CONST_X0P);
    COMPLEX16 PN3_tail = (CST_2PI*(e*(EIL(1)/8. + EILN(1)/8. + e*(EIL(2)/4. + EILN(2)/4. + e*(-0.015625*EIL(1) + (27*EIL(3))/64. - EILN(1)/64. + (27*EILN(3))/64. + e*(-0.08333333333333333*EIL(2) + (2*EIL(4))/3. - EILN(2)/12. + (2*EILN(4))/3. + e*(EIL(1)/1536. - (243*EIL(3))/1024. + (3125*EIL(5))/3072. + EILN(1)/1536. - (243*EILN(3))/1024. + (3125*EILN(5))/3072. + e*(EIL(2)/96. - (8*EIL(4))/15. + (243*EIL(6))/160. + EILN(2)/96. - (8*EILN(4))/15. + (243*EILN(6))/160. + e*(-0.000013563368055555555*EIL(1) + (2187*EIL(3))/40960. - (78125*EIL(5))/73728. + (823543*EIL(7))/368640. - EILN(1)/73728. + (2187*EILN(3))/40960. - (78125*EILN(5))/73728. + (823543*EILN(7))/368640. + e*(-0.0006944444444444445*EIL(2) + (8*EIL(4))/45. - (2187*EIL(6))/1120. + (1024*EIL(8))/315. - EILN(2)/1440. + (8*EILN(4))/45. - (2187*EILN(6))/1120. + (1024*EILN(8))/315. + e*(EIL(1)/5.89824e6 - (2187*EIL(3))/327680. + (1953125*EIL(5))/4.128768e6 - (40353607*EIL(7))/1.179648e7 + (43046721*EIL(9))/9.17504e6 + EILN(1)/5.89824e6 - (2187*EILN(3))/327680. + (1953125*EILN(5))/4.128768e6 - (40353607*EILN(7))/1.179648e7 + (43046721*EILN(9))/9.17504e6 + e*(EIL(2)/34560. - (32*EIL(4))/945. + (19683*EIL(6))/17920. - (16384*EIL(8))/2835. + (1953125*EIL(10))/290304. + EILN(2)/34560. - (32*EILN(4))/945. + (19683*EILN(6))/17920. - (16384*EILN(8))/2835. + (1953125*EILN(10))/290304.)))))))))))
        +I*(e*(lgxOx0p*((-3*EIL(1))/4. + (3*EILN(1))/4.) + (EIL(1)*CST_LN2)/2. - (EILN(1)*CST_LN2)/2. + e*(lgxOx0p*((-3*EIL(2))/2. + (3*EILN(2))/2.) + e*(lgxOx0p*((3*EIL(1))/32. - (81*EIL(3))/32. - (3*EILN(1))/32. + (81*EILN(3))/32.) - (EIL(1)*CST_LN2)/16. + (EILN(1)*CST_LN2)/16. + EIL(3)*((27*CST_LN2)/16. - (27*CST_LN3)/16.) + EILN(3)*((-27*CST_LN2)/16. + (27*CST_LN3)/16.) + e*(lgxOx0p*(EIL(2)/2. - 4*EIL(4) - EILN(2)/2. + 4*EILN(4)) - (8*EIL(4)*CST_LN2)/3. + (8*EILN(4)*CST_LN2)/3. + e*(lgxOx0p*(-0.00390625*EIL(1) + (729*EIL(3))/512. - (3125*EIL(5))/512. + EILN(1)/256. - (729*EILN(3))/512. + (3125*EILN(5))/512.) + (EIL(1)*CST_LN2)/384. - (EILN(1)*CST_LN2)/384. + EILN(3)*((243*CST_LN2)/256. - (243*CST_LN3)/256.) + EIL(3)*((-243*CST_LN2)/256. + (243*CST_LN3)/256.) + EIL(5)*((3125*CST_LN2)/768. - (3125*CST_LN5)/768.) + EILN(5)*((-3125*CST_LN2)/768. + (3125*CST_LN5)/768.) + e*(lgxOx0p*(-0.0625*EIL(2) + (16*EIL(4))/5. - (729*EIL(6))/80. + EILN(2)/16. - (16*EILN(4))/5. + (729*EILN(6))/80.) + (32*EIL(4)*CST_LN2)/15. - (32*EILN(4)*CST_LN2)/15. - (243*EIL(6)*CST_LN3)/40. + (243*EILN(6)*CST_LN3)/40. + e*(lgxOx0p*(EIL(1)/12288. - (6561*EIL(3))/20480. + (78125*EIL(5))/12288. - (823543*EIL(7))/61440. - EILN(1)/12288. + (6561*EILN(3))/20480. - (78125*EILN(5))/12288. + (823543*EILN(7))/61440.) - (EIL(1)*CST_LN2)/18432. + (EILN(1)*CST_LN2)/18432. + EIL(3)*((2187*CST_LN2)/10240. - (2187*CST_LN3)/10240.) + EILN(3)*((-2187*CST_LN2)/10240. + (2187*CST_LN3)/10240.) + EILN(5)*((78125*CST_LN2)/18432. - (78125*CST_LN5)/18432.) + EIL(5)*((-78125*CST_LN2)/18432. + (78125*CST_LN5)/18432.) + EIL(7)*((823543*CST_LN2)/92160. - (823543*CST_LN7)/92160.) + EILN(7)*((-823543*CST_LN2)/92160. + (823543*CST_LN7)/92160.) + e*(lgxOx0p*(EIL(2)/240. - (16*EIL(4))/15. + (6561*EIL(6))/560. - (2048*EIL(8))/105. - EILN(2)/240. + (16*EILN(4))/15. - (6561*EILN(6))/560. + (2048*EILN(8))/105.) - (32*EIL(4)*CST_LN2)/45. - (8192*EIL(8)*CST_LN2)/315. + (32*EILN(4)*CST_LN2)/45. + (8192*EILN(8)*CST_LN2)/315. + (2187*EIL(6)*CST_LN3)/280. - (2187*EILN(6)*CST_LN3)/280. + e*(lgxOx0p*(-1.0172526041666667e-6*EIL(1) + (6561*EIL(3))/163840. - (1953125*EIL(5))/688128. + (40353607*EIL(7))/1.96608e6 - (129140163*EIL(9))/4.58752e6 + EILN(1)/983040. - (6561*EILN(3))/163840. + (1953125*EILN(5))/688128. - (40353607*EILN(7))/1.96608e6 + (129140163*EILN(9))/4.58752e6) + (EIL(1)*CST_LN2)/1.47456e6 - (EILN(1)*CST_LN2)/1.47456e6 + EIL(9)*((43046721*CST_LN2)/2.29376e6 - (43046721*CST_LN3)/1.14688e6) + EILN(3)*((2187*CST_LN2)/81920. - (2187*CST_LN3)/81920.) + EIL(3)*((-2187*CST_LN2)/81920. + (2187*CST_LN3)/81920.) + EILN(9)*((-43046721*CST_LN2)/2.29376e6 + (43046721*CST_LN3)/1.14688e6) + EIL(5)*((1953125*CST_LN2)/1.032192e6 - (1953125*CST_LN5)/1.032192e6) + EILN(5)*((-1953125*CST_LN2)/1.032192e6 + (1953125*CST_LN5)/1.032192e6) + e*(lgxOx0p*(-0.00017361111111111112*EIL(2) + (64*EIL(4))/315. - (59049*EIL(6))/8960. + (32768*EIL(8))/945. - (1953125*EIL(10))/48384. + EILN(2)/5760. - (64*EILN(4))/315. + (59049*EILN(6))/8960. - (32768*EILN(8))/945. + (1953125*EILN(10))/48384.) + (128*EIL(4)*CST_LN2)/945. + (131072*EIL(8)*CST_LN2)/2835. - (128*EILN(4)*CST_LN2)/945. - (131072*EILN(8)*CST_LN2)/2835. - (19683*EIL(6)*CST_LN3)/4480. + (19683*EILN(6)*CST_LN3)/4480. - (1953125*EIL(10)*CST_LN5)/72576. + (1953125*EILN(10)*CST_LN5)/72576.) + EILN(7)*((40353607*CST_LN2)/2.94912e6 - (40353607*CST_LN7)/2.94912e6) + EIL(7)*((-40353607*CST_LN2)/2.94912e6 + (40353607*CST_LN7)/2.94912e6))))))))))));

                    
    COMPLEX16 PN5_tail = eta*(CST_2PI*(e*((-2*EIL(1))/21. - (2*EILN(1))/21. + e*((-13*EIL(2))/84. - (13*EILN(2))/84. + e*((5*EIL(1))/672. - (45*EIL(3))/224. + (5*EILN(1))/672. - (45*EILN(3))/224. + e*(EIL(2)/36. - (2*EIL(4))/9. + EILN(2)/36. - (2*EILN(4))/9. + e*(-0.0001240079365079365*EIL(1) + (81*EIL(3))/1792. - (3125*EIL(5))/16128. - EILN(1)/8064. + (81*EILN(3))/1792. - (3125*EILN(5))/16128. + e*(-0.000496031746031746*EIL(2) + (8*EIL(4))/315. - (81*EIL(6))/1120. - EILN(2)/2016. + (8*EILN(4))/315. - (81*EILN(6))/1120. + e*(-1.2917493386243386e-6*EIL(1) + (729*EIL(3))/143360. - (78125*EIL(5))/774144. + (117649*EIL(7))/552960. - EILN(1)/774144. + (729*EILN(3))/143360. - (78125*EILN(5))/774144. + (117649*EILN(7))/552960. + e*(-0.00016534391534391533*EIL(2) + (8*EIL(4))/189. - (729*EIL(6))/1568. + (1024*EIL(8))/1323. - EILN(2)/6048. + (8*EILN(4))/189. - (729*EILN(6))/1568. + (1024*EILN(8))/1323. + e*(EIL(1)/1.548288e7 - (729*EIL(3))/286720. + (1953125*EIL(5))/1.0838016e7 - (5764801*EIL(7))/4.42368e6 + (14348907*EIL(9))/8.02816e6 + EILN(1)/1.548288e7 - (729*EILN(3))/286720. + (1953125*EILN(5))/1.0838016e7 - (5764801*EILN(7))/4.42368e6 + (14348907*EILN(9))/8.02816e6 + e*((11*EIL(2))/725760. - (352*EIL(4))/19845. + (72171*EIL(6))/125440. - (180224*EIL(8))/59535. + (21484375*EIL(10))/6.096384e6 + (11*EILN(2))/725760. - (352*EILN(4))/19845. + (72171*EILN(6))/125440. - (180224*EILN(8))/59535. + (21484375*EILN(10))/6.096384e6)))))))))))
        +I*(e*(lgxOx0p*((4*EIL(1))/7. - (4*EILN(1))/7.) - (8*EIL(1)*CST_LN2)/21. + (8*EILN(1)*CST_LN2)/21. + e*(lgxOx0p*((13*EIL(2))/14. - (13*EILN(2))/14.) + e*(lgxOx0p*((-5*EIL(1))/112. + (135*EIL(3))/112. + (5*EILN(1))/112. - (135*EILN(3))/112.) + (5*EIL(1)*CST_LN2)/168. - (5*EILN(1)*CST_LN2)/168. + EILN(3)*((45*CST_LN2)/56. - (45*CST_LN3)/56.) + EIL(3)*((-45*CST_LN2)/56. + (45*CST_LN3)/56.) + e*(lgxOx0p*(-0.16666666666666666*EIL(2) + (4*EIL(4))/3. + EILN(2)/6. - (4*EILN(4))/3.) + (8*EIL(4)*CST_LN2)/9. - (8*EILN(4)*CST_LN2)/9. + e*(lgxOx0p*(EIL(1)/1344. - (243*EIL(3))/896. + (3125*EIL(5))/2688. - EILN(1)/1344. + (243*EILN(3))/896. - (3125*EILN(5))/2688.) - (EIL(1)*CST_LN2)/2016. + (EILN(1)*CST_LN2)/2016. + EIL(3)*((81*CST_LN2)/448. - (81*CST_LN3)/448.) + EILN(3)*((-81*CST_LN2)/448. + (81*CST_LN3)/448.) + EILN(5)*((3125*CST_LN2)/4032. - (3125*CST_LN5)/4032.) + EIL(5)*((-3125*CST_LN2)/4032. + (3125*CST_LN5)/4032.) + e*(lgxOx0p*(EIL(2)/336. - (16*EIL(4))/105. + (243*EIL(6))/560. - EILN(2)/336. + (16*EILN(4))/105. - (243*EILN(6))/560.) - (32*EIL(4)*CST_LN2)/315. + (32*EILN(4)*CST_LN2)/315. + (81*EIL(6)*CST_LN3)/280. - (81*EILN(6)*CST_LN3)/280. + e*(lgxOx0p*(EIL(1)/129024. - (2187*EIL(3))/71680. + (78125*EIL(5))/129024. - (117649*EIL(7))/92160. - EILN(1)/129024. + (2187*EILN(3))/71680. - (78125*EILN(5))/129024. + (117649*EILN(7))/92160.) - (EIL(1)*CST_LN2)/193536. + (EILN(1)*CST_LN2)/193536. + EIL(3)*((729*CST_LN2)/35840. - (729*CST_LN3)/35840.) + EILN(3)*((-729*CST_LN2)/35840. + (729*CST_LN3)/35840.) + EILN(5)*((78125*CST_LN2)/193536. - (78125*CST_LN5)/193536.) + EIL(5)*((-78125*CST_LN2)/193536. + (78125*CST_LN5)/193536.) + EIL(7)*((117649*CST_LN2)/138240. - (117649*CST_LN7)/138240.) + EILN(7)*((-117649*CST_LN2)/138240. + (117649*CST_LN7)/138240.) + e*(lgxOx0p*(EIL(2)/1008. - (16*EIL(4))/63. + (2187*EIL(6))/784. - (2048*EIL(8))/441. - EILN(2)/1008. + (16*EILN(4))/63. - (2187*EILN(6))/784. + (2048*EILN(8))/441.) - (32*EIL(4)*CST_LN2)/189. - (8192*EIL(8)*CST_LN2)/1323. + (32*EILN(4)*CST_LN2)/189. + (8192*EILN(8)*CST_LN2)/1323. + (729*EIL(6)*CST_LN3)/392. - (729*EILN(6)*CST_LN3)/392. + e*(lgxOx0p*(-3.875248015873016e-7*EIL(1) + (2187*EIL(3))/143360. - (1953125*EIL(5))/1.806336e6 + (5764801*EIL(7))/737280. - (43046721*EIL(9))/4.01408e6 + EILN(1)/2.58048e6 - (2187*EILN(3))/143360. + (1953125*EILN(5))/1.806336e6 - (5764801*EILN(7))/737280. + (43046721*EILN(9))/4.01408e6) + (EIL(1)*CST_LN2)/3.87072e6 - (EILN(1)*CST_LN2)/3.87072e6 + EIL(9)*((14348907*CST_LN2)/2.00704e6 - (14348907*CST_LN3)/1.00352e6) + EILN(3)*((729*CST_LN2)/71680. - (729*CST_LN3)/71680.) + EIL(3)*((-729*CST_LN2)/71680. + (729*CST_LN3)/71680.) + EILN(9)*((-14348907*CST_LN2)/2.00704e6 + (14348907*CST_LN3)/1.00352e6) + EIL(5)*((1953125*CST_LN2)/2.709504e6 - (1953125*CST_LN5)/2.709504e6) + EILN(5)*((-1953125*CST_LN2)/2.709504e6 + (1953125*CST_LN5)/2.709504e6) + e*(lgxOx0p*((-11*EIL(2))/120960. + (704*EIL(4))/6615. - (216513*EIL(6))/62720. + (360448*EIL(8))/19845. - (21484375*EIL(10))/1.016064e6 + (11*EILN(2))/120960. - (704*EILN(4))/6615. + (216513*EILN(6))/62720. - (360448*EILN(8))/19845. + (21484375*EILN(10))/1.016064e6) + (1408*EIL(4)*CST_LN2)/19845. + (1441792*EIL(8)*CST_LN2)/59535. - (1408*EILN(4)*CST_LN2)/19845. - (1441792*EILN(8)*CST_LN2)/59535. - (72171*EIL(6)*CST_LN3)/31360. + (72171*EILN(6)*CST_LN3)/31360. - (21484375*EIL(10)*CST_LN5)/1.524096e6 + (21484375*EILN(10)*CST_LN5)/1.524096e6) + EILN(7)*((5764801*CST_LN2)/1.10592e6 - (5764801*CST_LN7)/1.10592e6) + EIL(7)*((-5764801*CST_LN2)/1.10592e6 + (5764801*CST_LN7)/1.10592e6))))))))))))
        +(CST_2PI*(e*((-103*EIL(1))/112. - (103*EILN(1))/112. + e*((-155*EIL(2))/56. - (155*EILN(2))/56. + e*((-353*EIL(1))/896. - (5589*EIL(3))/896. - (353*EILN(1))/896. - (5589*EILN(3))/896. + e*((7*EIL(2))/24. - (37*EIL(4))/3. + (7*EILN(2))/24. - (37*EILN(4))/3. + e*((-12071*EIL(1))/21504. + (45333*EIL(3))/14336. - (971875*EIL(5))/43008. - (12071*EILN(1))/21504. + (45333*EILN(3))/14336. - (971875*EILN(5))/43008. + e*((-1483*EIL(2))/1344. + (1102*EIL(4))/105. - (88209*EIL(6))/2240. - (1483*EILN(2))/1344. + (1102*EILN(4))/105. - (88209*EILN(6))/2240. + e*((-567425*EIL(1))/1.032192e6 - (287361*EIL(3))/114688. + (27171875*EIL(5))/1.032192e6 - (9764867*EIL(7))/147456. - (567425*EILN(1))/1.032192e6 - (287361*EILN(3))/114688. + (27171875*EILN(5))/1.032192e6 - (9764867*EILN(7))/147456. + e*((-17383*EIL(2))/20160. - (2078*EIL(4))/315. + (902259*EIL(6))/15680. - (239104*EIL(8))/2205. - (17383*EILN(2))/20160. - (2078*EILN(4))/315. + (902259*EILN(6))/15680. - (239104*EILN(8))/2205. + e*((-45422119*EIL(1))/8.257536e7 - (4323267*EIL(3))/4.58752e6 - (1001421875*EIL(5))/5.7802752e7 + (2728397959*EIL(7))/2.359296e7 - (22341248199*EIL(9))/1.2845056e8 - (45422119*EILN(1))/8.257536e7 - (4323267*EILN(3))/4.58752e6 - (1001421875*EILN(5))/5.7802752e7 + (2728397959*EILN(7))/2.359296e7 - (22341248199*EILN(9))/1.2845056e8 + e*((-427291*EIL(2))/483840. - (1154*EIL(4))/6615. - (10694673*EIL(6))/250880. + (4355072*EIL(8))/19845. - (1115234375*EIL(10))/4.064256e6 - (427291*EILN(2))/483840. - (1154*EILN(4))/6615. - (10694673*EILN(6))/250880. + (4355072*EILN(8))/19845. - (1115234375*EILN(10))/4.064256e6)))))))))))
        +I*(e*((-309*lgxOx0p*EILN(1))/56. + EIL(1)*(1.5 + (309*lgxOx0p)/56. - (103*CST_LN2)/28.) + EILN(1)*(-1.5 + (103*CST_LN2)/28.) + e*(3*EIL(2) + lgxOx0p*((465*EIL(2))/28. - (465*EILN(2))/28.) - 3*EILN(2) + e*(lgxOx0p*((16767*EIL(3))/448. - (1059*EILN(1))/448. - (16767*EILN(3))/448.) + EIL(1)*(1.3125 + (1059*lgxOx0p)/448. - (353*CST_LN2)/224.) + EILN(1)*(-1.3125 + (353*CST_LN2)/224.) + EILN(3)*(-5.0625 + (5589*CST_LN2)/224. - (5589*CST_LN3)/224.) + EIL(3)*(5.0625 - (5589*CST_LN2)/224. + (5589*CST_LN3)/224.) + e*(2*EIL(2) - 2*EILN(2) + lgxOx0p*((-7*EIL(2))/4. + 74*EIL(4) + (7*EILN(2))/4. - 74*EILN(4)) + EILN(4)*(-8 - (148*CST_LN2)/3.) + EIL(4)*(8 + (148*CST_LN2)/3.) + e*(lgxOx0p*((-135999*EIL(3))/7168. + (971875*EIL(5))/7168. - (12071*EILN(1))/3584. + (135999*EILN(3))/7168. - (971875*EILN(5))/7168.) + EIL(1)*(1.3203125 + (12071*lgxOx0p)/3584. - (12071*CST_LN2)/5376.) + EILN(1)*(-1.3203125 + (12071*CST_LN2)/5376.) + EIL(3)*(2.21484375 + (45333*CST_LN2)/3584. - (45333*CST_LN3)/3584.) + EILN(3)*(-2.21484375 - (45333*CST_LN2)/3584. + (45333*CST_LN3)/3584.) + EILN(5)*(-12.20703125 + (971875*CST_LN2)/10752. - (971875*CST_LN5)/10752.) + EIL(5)*(12.20703125 - (971875*CST_LN2)/10752. + (971875*CST_LN5)/10752.) + e*((17*EIL(2))/8. - (17*EILN(2))/8. + lgxOx0p*((1483*EIL(2))/224. - (2204*EIL(4))/35. + (264627*EIL(6))/1120. - (1483*EILN(2))/224. + (2204*EILN(4))/35. - (264627*EILN(6))/1120.) + EIL(4)*(1.6 - (4408*CST_LN2)/105.) + EILN(4)*(-1.6 + (4408*CST_LN2)/105.) + EILN(6)*(-18.225 - (88209*CST_LN3)/560.) + EIL(6)*(18.225 + (88209*CST_LN3)/560.) + e*(lgxOx0p*((862083*EIL(3))/57344. - (27171875*EIL(5))/172032. + (9764867*EIL(7))/24576. - (567425*EILN(1))/172032. - (862083*EILN(3))/57344. + (27171875*EILN(5))/172032. - (9764867*EILN(7))/24576.) + EIL(1)*(1.3201497395833333 + (567425*lgxOx0p)/172032. - (567425*CST_LN2)/258048.) + EILN(1)*(-1.3201497395833333 + (567425*CST_LN2)/258048.) + EILN(3)*(-2.85556640625 + (287361*CST_LN2)/28672. - (287361*CST_LN3)/28672.) + EIL(3)*(2.85556640625 - (287361*CST_LN2)/28672. + (287361*CST_LN3)/28672.) + EIL(5)*(-0.5086263020833334 + (27171875*CST_LN2)/258048. - (27171875*CST_LN5)/258048.) + EILN(5)*(0.5086263020833334 - (27171875*CST_LN2)/258048. + (27171875*CST_LN5)/258048.) + EILN(7)*(-26.808040364583334 + (9764867*CST_LN2)/36864. - (9764867*CST_LN7)/36864.) + EIL(7)*(26.808040364583334 - (9764867*CST_LN2)/36864. + (9764867*CST_LN7)/36864.) + e*((127*EIL(2))/60. - (127*EILN(2))/60. + lgxOx0p*((17383*EIL(2))/3360. + (4156*EIL(4))/105. - (2706777*EIL(6))/7840. + (478208*EIL(8))/735. - (17383*EILN(2))/3360. - (4156*EILN(4))/105. + (2706777*EILN(6))/7840. - (478208*EILN(8))/735.) + EILN(8)*(-39.00952380952381 - (1912832*CST_LN2)/2205.) + EILN(4)*(-3.7333333333333334 - (8312*CST_LN2)/315.) + EIL(4)*(3.7333333333333334 + (8312*CST_LN2)/315.) + EIL(8)*(39.00952380952381 + (1912832*CST_LN2)/2205.) + EIL(6)*(-5.207142857142857 - (902259*CST_LN3)/3920.) + EILN(6)*(5.207142857142857 + (902259*CST_LN3)/3920.) + e*(lgxOx0p*((12969801*EIL(3))/2.29376e6 + (1001421875*EIL(5))/9.633792e6 - (2728397959*EIL(7))/3.93216e6 + (67023744597*EIL(9))/6.422528e7 - (45422119*EILN(1))/1.376256e7 - (12969801*EILN(3))/2.29376e6 - (1001421875*EILN(5))/9.633792e6 + (2728397959*EILN(7))/3.93216e6 - (67023744597*EILN(9))/6.422528e7) + EIL(1)*(1.3201517740885416 + (45422119*lgxOx0p)/1.376256e7 - (45422119*CST_LN2)/2.064384e7) + EILN(1)*(-1.3201517740885416 + (45422119*CST_LN2)/2.064384e7) + EILN(9)*(-56.300643048967636 + (22341248199*CST_LN2)/3.211264e7 - (22341248199*CST_LN3)/1.605632e7) + EILN(3)*(-2.77547607421875 + (4323267*CST_LN2)/1.14688e6 - (4323267*CST_LN3)/1.14688e6) + EIL(3)*(2.77547607421875 - (4323267*CST_LN2)/1.14688e6 + (4323267*CST_LN3)/1.14688e6) + EIL(9)*(56.300643048967636 - (22341248199*CST_LN2)/3.211264e7 + (22341248199*CST_LN3)/1.605632e7) + EILN(5)*(-5.168006533668155 + (1001421875*CST_LN2)/1.4450688e7 - (1001421875*CST_LN5)/1.4450688e7) + EIL(5)*(5.168006533668155 - (1001421875*CST_LN2)/1.4450688e7 + (1001421875*CST_LN5)/1.4450688e7) + e*((6097*EIL(2))/2880. - (6097*EILN(2))/2880. + lgxOx0p*((427291*EIL(2))/80640. + (2308*EIL(4))/2205. + (32084019*EIL(6))/125440. - (8710144*EIL(8))/6615. + (1115234375*EIL(10))/677376. - (427291*EILN(2))/80640. - (2308*EILN(4))/2205. - (32084019*EILN(6))/125440. + (8710144*EILN(8))/6615. - (1115234375*EILN(10))/677376.) + EIL(8)*(-30.340740740740742 - (34840576*CST_LN2)/19845.) + EILN(4)*(-3.326984126984127 - (4616*CST_LN2)/6615.) + EIL(4)*(3.326984126984127 + (4616*CST_LN2)/6615.) + EILN(8)*(30.340740740740742 + (34840576*CST_LN2)/19845.) + EILN(6)*(-7.9734375 - (10694673*CST_LN3)/62720.) + EIL(6)*(7.9734375 + (10694673*CST_LN3)/62720.) + EILN(10)*(-80.73433366402116 - (1115234375*CST_LN5)/1.016064e6) + EIL(10)*(80.73433366402116 + (1115234375*CST_LN5)/1.016064e6)) + EIL(7)*(-14.241771443684895 + (2728397959*CST_LN2)/5.89824e6 - (2728397959*CST_LN7)/5.89824e6) + EILN(7)*(14.241771443684895 - (2728397959*CST_LN2)/5.89824e6 + (2728397959*CST_LN7)/5.89824e6))))))))))));

    COMPLEX16 PN6_tail = ((e*(EILN(1)*(-1.4599291383219954 + CST_PISQ/24. + (0.5095238095238095 - CST_LN2/2.)*CST_LN2) + lgxOx0p*((-9*lgxOx0p*EILN(1))/8. + EILN(1)*(-0.7642857142857142 + (3*CST_LN2)/2.)) + EIL(1)*(-1.4599291383219954 + CST_PISQ/24. + (0.5095238095238095 - CST_LN2/2.)*CST_LN2 + lgxOx0p*(-0.7642857142857142 - (9*lgxOx0p)/8. + (3*CST_LN2)/2.)) + e*((-5.839716553287982 + CST_PISQ/6.)*EIL(2) + lgxOx0p*((-107*EIL(2))/35. + lgxOx0p*((-9*EIL(2))/2. - (9*EILN(2))/2.) - (107*EILN(2))/35.) + (-5.839716553287982 + CST_PISQ/6.)*EILN(2) + e*(EILN(1)*(0.18249114229024943 - CST_PISQ/192. + (-0.06369047619047619 + CST_LN2/16.)*CST_LN2) + EIL(1)*(0.18249114229024943 - CST_PISQ/192. + lgxOx0p*(0.09553571428571428 + (9*lgxOx0p)/64. - (3*CST_LN2)/16.) + (-0.06369047619047619 + CST_LN2/16.)*CST_LN2) + lgxOx0p*(lgxOx0p*((-729*EIL(3))/64. + (9*EILN(1))/64. - (729*EILN(3))/64.) + EILN(1)*(0.09553571428571428 - (3*CST_LN2)/16.) + EIL(3)*(-7.738392857142857 + (243*CST_LN2)/16. - (243*CST_LN3)/16.) + EILN(3)*(-7.738392857142857 + (243*CST_LN2)/16. - (243*CST_LN3)/16.)) + EIL(3)*(-14.781782525510204 + (27*CST_PISQ)/64. + (-5.1589285714285715 - (81*CST_LN3)/16.)*CST_LN3 + CST_LN2*(5.1589285714285715 - (81*CST_LN2)/16. + (81*CST_LN3)/8.)) + EILN(3)*(-14.781782525510204 + (27*CST_PISQ)/64. + (-5.1589285714285715 - (81*CST_LN3)/16.)*CST_LN3 + CST_LN2*(5.1589285714285715 - (81*CST_LN2)/16. + (81*CST_LN3)/8.)) + e*((1.9465721844293273 - CST_PISQ/18.)*EIL(2) + (1.9465721844293273 - CST_PISQ/18.)*EILN(2) + lgxOx0p*((107*EIL(2))/105. + (107*EILN(2))/105. + lgxOx0p*((3*EIL(2))/2. - 24*EIL(4) + (3*EILN(2))/2. - 24*EILN(4)) + EIL(4)*(-16.304761904761904 - 32*CST_LN2) + EILN(4)*(-16.304761904761904 - 32*CST_LN2)) + EIL(4)*(-31.145154950869237 + (8*CST_PISQ)/9. + (-10.86984126984127 - (32*CST_LN2)/3.)*CST_LN2) + EILN(4)*(-31.145154950869237 + (8*CST_PISQ)/9. + (-10.86984126984127 - (32*CST_LN2)/3.)*CST_LN2) + e*(EILN(1)*(-0.00760379759542706 + CST_PISQ/4608. + (0.0026537698412698414 - CST_LN2/384.)*CST_LN2) + EIL(1)*(-0.00760379759542706 + CST_PISQ/4608. + lgxOx0p*(-0.003980654761904762 - (3*lgxOx0p)/512. + CST_LN2/128.) + (0.0026537698412698414 - CST_LN2/384.)*CST_LN2) + EIL(3)*(8.31475267059949 - (243*CST_PISQ)/1024. + CST_LN2*(-2.9018973214285713 + (729*CST_LN2)/256. - (729*CST_LN3)/128.) + CST_LN3*(2.9018973214285713 + (729*CST_LN3)/256.)) + EILN(3)*(8.31475267059949 - (243*CST_PISQ)/1024. + CST_LN2*(-2.9018973214285713 + (729*CST_LN2)/256. - (729*CST_LN3)/128.) + CST_LN3*(2.9018973214285713 + (729*CST_LN3)/256.)) + lgxOx0p*(lgxOx0p*((6561*EIL(3))/1024. - (46875*EIL(5))/1024. - (3*EILN(1))/512. + (6561*EILN(3))/1024. - (46875*EILN(5))/1024.) + EILN(1)*(-0.003980654761904762 + CST_LN2/128.) + EIL(3)*(4.352845982142857 - (2187*CST_LN2)/256. + (2187*CST_LN3)/256.) + EILN(3)*(4.352845982142857 - (2187*CST_LN2)/256. + (2187*CST_LN3)/256.) + EIL(5)*(-31.098865327380953 + (15625*CST_LN2)/256. - (15625*CST_LN5)/256.) + EILN(5)*(-31.098865327380953 + (15625*CST_LN2)/256. - (15625*CST_LN5)/256.)) + EIL(5)*(-59.4046687142739 + (15625*CST_PISQ)/9216. + (-20.732576884920636 - (15625*CST_LN5)/768.)*CST_LN5 + CST_LN2*(20.732576884920636 - (15625*CST_LN2)/768. + (15625*CST_LN5)/384.)) + EILN(5)*(-59.4046687142739 + (15625*CST_PISQ)/9216. + (-20.732576884920636 - (15625*CST_LN5)/768.)*CST_LN5 + CST_LN2*(20.732576884920636 - (15625*CST_LN2)/768. + (15625*CST_LN5)/384.)) + e*((-0.24332152305366592 + CST_PISQ/144.)*EIL(2) + (-0.24332152305366592 + CST_PISQ/144.)*EILN(2) + EIL(4)*(24.91612396069539 - (32*CST_PISQ)/45. + CST_LN2*(8.695873015873016 + (128*CST_LN2)/15.)) + EILN(4)*(24.91612396069539 - (32*CST_PISQ)/45. + CST_LN2*(8.695873015873016 + (128*CST_LN2)/15.)) + lgxOx0p*((-107*EIL(2))/840. - (107*EILN(2))/840. + lgxOx0p*((-3*EIL(2))/16. + (96*EIL(4))/5. - (6561*EIL(6))/80. - (3*EILN(2))/16. + (96*EILN(4))/5. - (6561*EILN(6))/80.) + EIL(4)*(13.043809523809523 + (128*CST_LN2)/5.) + EILN(4)*(13.043809523809523 + (128*CST_LN2)/5.) + EIL(6)*(-55.71642857142857 - (2187*CST_LN3)/20.) + EILN(6)*(-55.71642857142857 - (2187*CST_LN3)/20.)) + EIL(6)*(-106.42883418367347 + (243*CST_PISQ)/80. + (-37.144285714285715 - (729*CST_LN3)/20.)*CST_LN3) + EILN(6)*(-106.42883418367347 + (243*CST_PISQ)/80. + (-37.144285714285715 - (729*CST_LN3)/20.)*CST_LN3) + e*(EILN(1)*(0.00015841244990473041 - CST_PISQ/221184. + (-0.000055286871693121695 + CST_LN2/18432.)*CST_LN2) + EIL(1)*(0.00015841244990473041 - CST_PISQ/221184. + lgxOx0p*(0.00008293030753968254 + lgxOx0p/8192. - CST_LN2/6144.) + (-0.000055286871693121695 + CST_LN2/18432.)*CST_LN2) + EIL(3)*(-1.8708193508848852 + (2187*CST_PISQ)/40960. + (-0.6529268973214286 - (6561*CST_LN3)/10240.)*CST_LN3 + CST_LN2*(0.6529268973214286 - (6561*CST_LN2)/10240. + (6561*CST_LN3)/5120.)) + EILN(3)*(-1.8708193508848852 + (2187*CST_PISQ)/40960. + (-0.6529268973214286 - (6561*CST_LN3)/10240.)*CST_LN3 + CST_LN2*(0.6529268973214286 - (6561*CST_LN2)/10240. + (6561*CST_LN3)/5120.)) + EIL(5)*(61.87986324403531 - (390625*CST_PISQ)/221184. + CST_LN2*(-21.59643425512566 + (390625*CST_LN2)/18432. - (390625*CST_LN5)/9216.) + CST_LN5*(21.59643425512566 + (390625*CST_LN5)/18432.)) + EILN(5)*(61.87986324403531 - (390625*CST_PISQ)/221184. + CST_LN2*(-21.59643425512566 + (390625*CST_LN2)/18432. - (390625*CST_LN5)/9216.) + CST_LN5*(21.59643425512566 + (390625*CST_LN5)/18432.)) + lgxOx0p*(lgxOx0p*((-59049*EIL(3))/40960. + (390625*EIL(5))/8192. - (5764801*EIL(7))/40960. + EILN(1)/8192. - (59049*EILN(3))/40960. + (390625*EILN(5))/8192. - (5764801*EILN(7))/40960.) + EILN(1)*(0.00008293030753968254 - CST_LN2/6144.) + EIL(3)*(-0.9793903459821428 + (19683*CST_LN2)/10240. - (19683*CST_LN3)/10240.) + EILN(3)*(-0.9793903459821428 + (19683*CST_LN2)/10240. - (19683*CST_LN3)/10240.) + EIL(5)*(32.394651382688494 - (390625*CST_LN2)/6144. + (390625*CST_LN5)/6144.) + EILN(5)*(32.394651382688494 - (390625*CST_LN2)/6144. + (390625*CST_LN5)/6144.) + EIL(7)*(-95.61534396701389 + (5764801*CST_LN2)/30720. - (5764801*CST_LN7)/30720.) + EILN(7)*(-95.61534396701389 + (5764801*CST_LN2)/30720. - (5764801*CST_LN7)/30720.)) + EIL(7)*(-182.64324992464796 + (5764801*CST_PISQ)/1.10592e6 + (-63.743562644675926 - (5764801*CST_LN7)/92160.)*CST_LN7 + CST_LN2*(63.743562644675926 - (5764801*CST_LN2)/92160. + (5764801*CST_LN7)/46080.)) + EILN(7)*(-182.64324992464796 + (5764801*CST_PISQ)/1.10592e6 + (-63.743562644675926 - (5764801*CST_LN7)/92160.)*CST_LN7 + CST_LN2*(63.743562644675926 - (5764801*CST_LN2)/92160. + (5764801*CST_LN7)/46080.)) + e*((0.016221434870244394 - CST_PISQ/2160.)*EIL(2) + (0.016221434870244394 - CST_PISQ/2160.)*EILN(2) + EIL(8)*(-303.73941590181045 + (8192*CST_PISQ)/945. + (-212.01366591080875 - (131072*CST_LN2)/315.)*CST_LN2) + EILN(8)*(-303.73941590181045 + (8192*CST_PISQ)/945. + (-212.01366591080875 - (131072*CST_LN2)/315.)*CST_LN2) + EIL(4)*(-8.30537465356513 + (32*CST_PISQ)/135. + (-2.8986243386243387 - (128*CST_LN2)/45.)*CST_LN2) + EILN(4)*(-8.30537465356513 + (32*CST_PISQ)/135. + (-2.8986243386243387 - (128*CST_LN2)/45.)*CST_LN2) + EIL(6)*(136.83707252186588 - (2187*CST_PISQ)/560. + CST_LN3*(47.75693877551021 + (6561*CST_LN3)/140.)) + EILN(6)*(136.83707252186588 - (2187*CST_PISQ)/560. + CST_LN3*(47.75693877551021 + (6561*CST_LN3)/140.)) + lgxOx0p*((107*EIL(2))/12600. + (107*EILN(2))/12600. + lgxOx0p*(EIL(2)/80. - (32*EIL(4))/5. + (59049*EIL(6))/560. - (8192*EIL(8))/35. + EILN(2)/80. - (32*EILN(4))/5. + (59049*EILN(6))/560. - (8192*EILN(8))/35.) + EIL(8)*(-159.01024943310657 - (65536*CST_LN2)/105.) + EILN(8)*(-159.01024943310657 - (65536*CST_LN2)/105.) + EIL(4)*(-4.347936507936508 - (128*CST_LN2)/15.) + EILN(4)*(-4.347936507936508 - (128*CST_LN2)/15.) + EIL(6)*(71.63540816326531 + (19683*CST_LN3)/140.) + EILN(6)*(71.63540816326531 + (19683*CST_LN3)/140.)) + e*(EILN(1)*(-1.98015562380913e-6 + CST_PISQ/1.769472e7 + (6.910858961640212e-7 - CST_LN2/1.47456e6)*CST_LN2) + EIL(1)*(-1.98015562380913e-6 + CST_PISQ/1.769472e7 + lgxOx0p*(-1.0366288442460317e-6 - lgxOx0p/655360. + CST_LN2/491520.) + (6.910858961640212e-7 - CST_LN2/1.47456e6)*CST_LN2) + EIL(3)*(0.23385241886061064 - (2187*CST_PISQ)/327680. + CST_LN2*(-0.08161586216517858 + (6561*CST_LN2)/81920. - (6561*CST_LN3)/40960.) + (0.08161586216517858 + (6561*CST_LN3)/81920.)*CST_LN3) + EILN(3)*(0.23385241886061064 - (2187*CST_PISQ)/327680. + CST_LN2*(-0.08161586216517858 + (6561*CST_LN2)/81920. - (6561*CST_LN3)/40960.) + (0.08161586216517858 + (6561*CST_LN3)/81920.)*CST_LN3) + EIL(9)*(-493.16969576072137 + (129140163*CST_PISQ)/9.17504e6 + (-344.2382174994021 - (387420489*CST_LN3)/573440.)*CST_LN3 + CST_LN2*(172.11910874970104 - (387420489*CST_LN2)/2.29376e6 + (387420489*CST_LN3)/573440.)) + EILN(9)*(-493.16969576072137 + (129140163*CST_PISQ)/9.17504e6 + (-344.2382174994021 - (387420489*CST_LN3)/573440.)*CST_LN3 + CST_LN2*(172.11910874970104 - (387420489*CST_LN2)/2.29376e6 + (387420489*CST_LN3)/573440.)) + EIL(5)*(-27.624938948230053 + (9765625*CST_PISQ)/1.2386304e7 + (-9.641265292466812 - (9765625*CST_LN5)/1.032192e6)*CST_LN5 + CST_LN2*(9.641265292466812 - (9765625*CST_LN2)/1.032192e6 + (9765625*CST_LN5)/516096.)) + EILN(5)*(-27.624938948230053 + (9765625*CST_PISQ)/1.2386304e7 + (-9.641265292466812 - (9765625*CST_LN5)/1.032192e6)*CST_LN5 + CST_LN2*(9.641265292466812 - (9765625*CST_LN2)/1.032192e6 + (9765625*CST_LN5)/516096.)) + e*((-0.0006758931195935164 + CST_PISQ/51840.)*EIL(2) + (-0.0006758931195935164 + CST_PISQ/51840.)*EILN(2) + EIL(4)*(1.5819761244885961 - (128*CST_PISQ)/2835. + (0.5521189216427311 + (512*CST_LN2)/945.)*CST_LN2) + EILN(4)*(1.5819761244885961 - (128*CST_PISQ)/2835. + (0.5521189216427311 + (512*CST_LN2)/945.)*CST_LN2) + EIL(8)*(539.9811838254408 - (131072*CST_PISQ)/8505. + CST_LN2*(376.9131838414378 + (2097152*CST_LN2)/2835.)) + EILN(8)*(539.9811838254408 - (131072*CST_PISQ)/8505. + CST_LN2*(376.9131838414378 + (2097152*CST_LN2)/2835.)) + EIL(6)*(-76.97085329354957 + (19683*CST_PISQ)/8960. + (-26.86327806122449 - (59049*CST_LN3)/2240.)*CST_LN3) + EILN(6)*(-76.97085329354957 + (19683*CST_PISQ)/8960. + (-26.86327806122449 - (59049*CST_LN3)/2240.)*CST_LN3) + lgxOx0p*((-107*EIL(2))/302400. - (107*EILN(2))/302400. + lgxOx0p*(-0.0005208333333333333*EIL(2) + (128*EIL(4))/105. - (531441*EIL(6))/8960. + (131072*EIL(8))/315. - (9765625*EIL(10))/16128. - EILN(2)/1920. + (128*EILN(4))/105. - (531441*EILN(6))/8960. + (131072*EILN(8))/315. - (9765625*EILN(10))/16128.) + EIL(4)*(0.8281783824640967 + (512*CST_LN2)/315.) + EILN(4)*(0.8281783824640967 + (512*CST_LN2)/315.) + EIL(8)*(282.68488788107834 + (1048576*CST_LN2)/945.) + EILN(8)*(282.68488788107834 + (1048576*CST_LN2)/945.) + EIL(6)*(-40.29491709183674 - (177147*CST_LN3)/2240.) + EILN(6)*(-40.29491709183674 - (177147*CST_LN3)/2240.) + EIL(10)*(-411.36065247858403 - (9765625*CST_LN5)/12096.) + EILN(10)*(-411.36065247858403 - (9765625*CST_LN5)/12096.)) + EIL(10)*(-785.7760411940992 + (9765625*CST_PISQ)/435456. + (-274.2404349857227 - (9765625*CST_LN5)/36288.)*CST_LN5) + EILN(10)*(-785.7760411940992 + (9765625*CST_PISQ)/435456. + (-274.2404349857227 - (9765625*CST_LN5)/36288.)*CST_LN5)) + EIL(7)*(279.67247644711716 - (282475249*CST_PISQ)/3.538944e7 + CST_LN2*(-97.60733029966 + (282475249*CST_LN2)/2.94912e6 - (282475249*CST_LN7)/1.47456e6) + CST_LN7*(97.60733029966 + (282475249*CST_LN7)/2.94912e6)) + EILN(7)*(279.67247644711716 - (282475249*CST_PISQ)/3.538944e7 + CST_LN2*(-97.60733029966 + (282475249*CST_LN2)/2.94912e6 - (282475249*CST_LN7)/1.47456e6) + CST_LN7*(97.60733029966 + (282475249*CST_LN7)/2.94912e6)) + lgxOx0p*(lgxOx0p*((59049*EIL(3))/327680. - (9765625*EIL(5))/458752. + (282475249*EIL(7))/1.31072e6 - (3486784401*EIL(9))/9.17504e6 - EILN(1)/655360. + (59049*EILN(3))/327680. - (9765625*EILN(5))/458752. + (282475249*EILN(7))/1.31072e6 - (3486784401*EILN(9))/9.17504e6) + EILN(1)*(-1.0366288442460317e-6 + CST_LN2/491520.) + EIL(9)*(-258.1786631245516 + (1162261467*CST_LN2)/2.29376e6 - (1162261467*CST_LN3)/1.14688e6) + EILN(9)*(-258.1786631245516 + (1162261467*CST_LN2)/2.29376e6 - (1162261467*CST_LN3)/1.14688e6) + EIL(3)*(0.12242379324776785 - (19683*CST_LN2)/81920. + (19683*CST_LN3)/81920.) + EILN(3)*(0.12242379324776785 - (19683*CST_LN2)/81920. + (19683*CST_LN3)/81920.) + EIL(5)*(-14.46189793870022 + (9765625*CST_LN2)/344064. - (9765625*CST_LN5)/344064.) + EILN(5)*(-14.46189793870022 + (9765625*CST_LN2)/344064. - (9765625*CST_LN5)/344064.) + EIL(7)*(146.41099544949003 - (282475249*CST_LN2)/983040. + (282475249*CST_LN7)/983040.) + EILN(7)*(146.41099544949003 - (282475249*CST_LN2)/983040. + (282475249*CST_LN7)/983040.))))))))))))
        +I*(e*((3*lgxOx0p*CST_PI*EILN(1))/4. + EILN(1)*((107*CST_PI)/420. - (CST_PI*CST_LN2)/2.) + EIL(1)*((-107*CST_PI)/420. - (3*lgxOx0p*CST_PI)/4. + (CST_PI*CST_LN2)/2.) + e*((-107*CST_PI*EIL(2))/105. + (107*CST_PI*EILN(2))/105. + lgxOx0p*(-3*CST_PI*EIL(2) + 3*CST_PI*EILN(2)) + e*(lgxOx0p*((-243*CST_PI*EIL(3))/32. - (3*CST_PI*EILN(1))/32. + (243*CST_PI*EILN(3))/32.) + EIL(1)*((107*CST_PI)/3360. + (3*lgxOx0p*CST_PI)/32. - (CST_PI*CST_LN2)/16.) + EILN(1)*((-107*CST_PI)/3360. + (CST_PI*CST_LN2)/16.) + EIL(3)*((-2889*CST_PI)/1120. + (81*CST_PI*CST_LN2)/16. - (81*CST_PI*CST_LN3)/16.) + EILN(3)*((2889*CST_PI)/1120. - (81*CST_PI*CST_LN2)/16. + (81*CST_PI*CST_LN3)/16.) + e*((107*CST_PI*EIL(2))/315. - (107*CST_PI*EILN(2))/315. + lgxOx0p*(CST_PI*EIL(2) - 16*CST_PI*EIL(4) - CST_PI*EILN(2) + 16*CST_PI*EILN(4)) + EIL(4)*((-1712*CST_PI)/315. - (32*CST_PI*CST_LN2)/3.) + EILN(4)*((1712*CST_PI)/315. + (32*CST_PI*CST_LN2)/3.) + e*(lgxOx0p*((2187*CST_PI*EIL(3))/512. - (15625*CST_PI*EIL(5))/512. + (CST_PI*EILN(1))/256. - (2187*CST_PI*EILN(3))/512. + (15625*CST_PI*EILN(5))/512.) + EILN(1)*((107*CST_PI)/80640. - (CST_PI*CST_LN2)/384.) + EIL(1)*((-107*CST_PI)/80640. - (lgxOx0p*CST_PI)/256. + (CST_PI*CST_LN2)/384.) + EILN(3)*((-26001*CST_PI)/17920. + (729*CST_PI*CST_LN2)/256. - (729*CST_PI*CST_LN3)/256.) + EIL(3)*((26001*CST_PI)/17920. - (729*CST_PI*CST_LN2)/256. + (729*CST_PI*CST_LN3)/256.) + EIL(5)*((-334375*CST_PI)/32256. + (15625*CST_PI*CST_LN2)/768. - (15625*CST_PI*CST_LN5)/768.) + EILN(5)*((334375*CST_PI)/32256. - (15625*CST_PI*CST_LN2)/768. + (15625*CST_PI*CST_LN5)/768.) + e*((-107*CST_PI*EIL(2))/2520. + (107*CST_PI*EILN(2))/2520. + lgxOx0p*(-0.125*(CST_PI*EIL(2)) + (64*CST_PI*EIL(4))/5. - (2187*CST_PI*EIL(6))/40. + (CST_PI*EILN(2))/8. - (64*CST_PI*EILN(4))/5. + (2187*CST_PI*EILN(6))/40.) + EILN(4)*((-6848*CST_PI)/1575. - (128*CST_PI*CST_LN2)/15.) + EIL(4)*((6848*CST_PI)/1575. + (128*CST_PI*CST_LN2)/15.) + EIL(6)*((-26001*CST_PI)/1400. - (729*CST_PI*CST_LN3)/20.) + EILN(6)*((26001*CST_PI)/1400. + (729*CST_PI*CST_LN3)/20.) + e*(lgxOx0p*((-19683*CST_PI*EIL(3))/20480. + (390625*CST_PI*EIL(5))/12288. - (5764801*CST_PI*EIL(7))/61440. - (CST_PI*EILN(1))/12288. + (19683*CST_PI*EILN(3))/20480. - (390625*CST_PI*EILN(5))/12288. + (5764801*CST_PI*EILN(7))/61440.) + EIL(1)*((107*CST_PI)/3.87072e6 + (lgxOx0p*CST_PI)/12288. - (CST_PI*CST_LN2)/18432.) + EILN(1)*((-107*CST_PI)/3.87072e6 + (CST_PI*CST_LN2)/18432.) + EIL(3)*((-234009*CST_PI)/716800. + (6561*CST_PI*CST_LN2)/10240. - (6561*CST_PI*CST_LN3)/10240.) + EILN(3)*((234009*CST_PI)/716800. - (6561*CST_PI*CST_LN2)/10240. + (6561*CST_PI*CST_LN3)/10240.) + EILN(5)*((-8359375*CST_PI)/774144. + (390625*CST_PI*CST_LN2)/18432. - (390625*CST_PI*CST_LN5)/18432.) + EIL(5)*((8359375*CST_PI)/774144. - (390625*CST_PI*CST_LN2)/18432. + (390625*CST_PI*CST_LN5)/18432.) + EIL(7)*((-88119101*CST_PI)/2.7648e6 + (5764801*CST_PI*CST_LN2)/92160. - (5764801*CST_PI*CST_LN7)/92160.) + EILN(7)*((88119101*CST_PI)/2.7648e6 - (5764801*CST_PI*CST_LN2)/92160. + (5764801*CST_PI*CST_LN7)/92160.) + e*((107*CST_PI*EIL(2))/37800. - (107*CST_PI*EILN(2))/37800. + lgxOx0p*((CST_PI*EIL(2))/120. - (64*CST_PI*EIL(4))/15. + (19683*CST_PI*EIL(6))/280. - (16384*CST_PI*EIL(8))/105. - (CST_PI*EILN(2))/120. + (64*CST_PI*EILN(4))/15. - (19683*CST_PI*EILN(6))/280. + (16384*CST_PI*EILN(8))/105.) + EIL(8)*((-1753088*CST_PI)/33075. - (65536*CST_PI*CST_LN2)/315.) + EIL(4)*((-6848*CST_PI)/4725. - (128*CST_PI*CST_LN2)/45.) + EILN(4)*((6848*CST_PI)/4725. + (128*CST_PI*CST_LN2)/45.) + EILN(8)*((1753088*CST_PI)/33075. + (65536*CST_PI*CST_LN2)/315.) + EILN(6)*((-234009*CST_PI)/9800. - (6561*CST_PI*CST_LN3)/140.) + EIL(6)*((234009*CST_PI)/9800. + (6561*CST_PI*CST_LN3)/140.) + e*(lgxOx0p*((19683*CST_PI*EIL(3))/163840. - (9765625*CST_PI*EIL(5))/688128. + (282475249*CST_PI*EIL(7))/1.96608e6 - (1162261467*CST_PI*EIL(9))/4.58752e6 + (CST_PI*EILN(1))/983040. - (19683*CST_PI*EILN(3))/163840. + (9765625*CST_PI*EILN(5))/688128. - (282475249*CST_PI*EILN(7))/1.96608e6 + (1162261467*CST_PI*EILN(9))/4.58752e6) + EILN(1)*((107*CST_PI)/3.096576e8 - (CST_PI*CST_LN2)/1.47456e6) + EIL(1)*((-107*CST_PI)/3.096576e8 - (lgxOx0p*CST_PI)/983040. + (CST_PI*CST_LN2)/1.47456e6) + EIL(9)*((-13817997441*CST_PI)/1.605632e8 + (387420489*CST_PI*CST_LN2)/2.29376e6 - (387420489*CST_PI*CST_LN3)/1.14688e6) + EILN(3)*((-234009*CST_PI)/5.7344e6 + (6561*CST_PI*CST_LN2)/81920. - (6561*CST_PI*CST_LN3)/81920.) + EIL(3)*((234009*CST_PI)/5.7344e6 - (6561*CST_PI*CST_LN2)/81920. + (6561*CST_PI*CST_LN3)/81920.) + EILN(9)*((13817997441*CST_PI)/1.605632e8 - (387420489*CST_PI*CST_LN2)/2.29376e6 + (387420489*CST_PI*CST_LN3)/1.14688e6) + EIL(5)*((-208984375*CST_PI)/4.3352064e7 + (9765625*CST_PI*CST_LN2)/1.032192e6 - (9765625*CST_PI*CST_LN5)/1.032192e6) + EILN(5)*((208984375*CST_PI)/4.3352064e7 - (9765625*CST_PI*CST_LN2)/1.032192e6 + (9765625*CST_PI*CST_LN5)/1.032192e6) + e*((-107*CST_PI*EIL(2))/907200. + (107*CST_PI*EILN(2))/907200. + lgxOx0p*(-0.00034722222222222224*(CST_PI*EIL(2)) + (256*CST_PI*EIL(4))/315. - (177147*CST_PI*EIL(6))/4480. + (262144*CST_PI*EIL(8))/945. - (9765625*CST_PI*EIL(10))/24192. + (CST_PI*EILN(2))/2880. - (256*CST_PI*EILN(4))/315. + (177147*CST_PI*EILN(6))/4480. - (262144*CST_PI*EILN(8))/945. + (9765625*CST_PI*EILN(10))/24192.) + EILN(8)*((-28049408*CST_PI)/297675. - (1048576*CST_PI*CST_LN2)/2835.) + EILN(4)*((-27392*CST_PI)/99225. - (512*CST_PI*CST_LN2)/945.) + EIL(4)*((27392*CST_PI)/99225. + (512*CST_PI*CST_LN2)/945.) + EIL(8)*((28049408*CST_PI)/297675. + (1048576*CST_PI*CST_LN2)/2835.) + EIL(6)*((-2106081*CST_PI)/156800. - (59049*CST_PI*CST_LN3)/2240.) + EILN(6)*((2106081*CST_PI)/156800. + (59049*CST_PI*CST_LN3)/2240.) + EIL(10)*((-208984375*CST_PI)/1.524096e6 - (9765625*CST_PI*CST_LN5)/36288.) + EILN(10)*((208984375*CST_PI)/1.524096e6 + (9765625*CST_PI*CST_LN5)/36288.)) + EILN(7)*((-4317835949*CST_PI)/8.84736e7 + (282475249*CST_PI*CST_LN2)/2.94912e6 - (282475249*CST_PI*CST_LN7)/2.94912e6) + EIL(7)*((4317835949*CST_PI)/8.84736e7 - (282475249*CST_PI*CST_LN2)/2.94912e6 + (282475249*CST_PI*CST_LN7)/2.94912e6))))))))))));


#undef EIL
#undef EILN

    return SQ2O3*(PN0 + v2*(PN2 + v*(PN3_tail + v*(PN4 + v*(I*ImPN5 + PN5_tail + v*(PN6 + PN6_tail))))));
}



DEFFUNC(2,2)
DEFFUNC(2,1)
DEFFUNC(2,0)
#undef DEFFUNC