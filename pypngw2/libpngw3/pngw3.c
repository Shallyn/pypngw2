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

static inline REAL8 evaluate_dchi(REAL8 chi, REAL8 beta)
{
    REAL8 x = 1. - beta*cos(chi);
    REAL8 y = sin(chi)*beta;
    return 2.0*atan2(y, x);
}

static int evaluate_PACorrections(REAL8 eta, REAL8 chi, DynParamsCache *dc,
    REAL8 *ret_eOsc, REAL8 *ret_vOsc, REAL8 *ret_nOsc,
    REAL8 *ret_clOsc, REAL8 *ret_clbdOsc)
{
    REAL8 *ePow = dc->ePow, *sq1mesqPow = dc->sq1mesqPow;
    REAL8 beta = eval_beta(ePow[1]);
    REAL8 dchi = evaluate_dchi(chi, beta);
    REAL8 schi = sin(chi), cchi = cos(chi);
    REAL8 up = 1./(1. - ePow[1]*cchi);
    REAL8 up2 = up*up;
    REAL8 up3 = up2*up;
    REAL8 up4 = up3*up;
    REAL8 up5 = up4*up;
    REAL8 up6 = up5*up;
    REAL8 up7 = up6*up;

    // e
    REAL8 eOsc_PN5, eOsc_PN7;
    eOsc_PN5 = -2.*schi*up3*(49. + 36.*ePow[2]) - (schi*up2*(134. + 291.*ePow[2]))/sq1mesqPow[2] - 
        210.*schi*up4*sq1mesqPow[2] - (schi*up*(134. + 1069.*ePow[2] + 72.*ePow[4]))/sq1mesqPow[4] - 
        (3.*ePow[1]*(dchi + schi*ePow[1])*(304. + 121.*ePow[2]))/sq1mesqPow[5];
    eOsc_PN5 /= 45.;
    eOsc_PN7 = 42*schi*up4*(28296 - 18700*eta + 5*(-6957 + 7712*eta)*ePow[2]) + 
        (2*schi*up3*(4*(9204 + 49525*eta) + (-820497 + 512960*eta)*ePow[2] + 360*(-573 + 728*eta)*ePow[4]))/
        sq1mesqPow[2] + 1512*schi*up5*(-939 + 1000*eta + 30*(19 + 2*eta)*ePow[2])*sq1mesqPow[2] + 
        (schi*up2*(484032 - 1960*eta + (-3739182 + 3600940*eta)*ePow[2] + 15*(-112977 + 84448*eta)*ePow[4]))/
        sq1mesqPow[4] - 2520*(137 + 140*eta)*schi*up6*sq1mesqPow[4] + 
        (schi*up*(484032 - 1960*eta + 14*(-393381 + 502390*eta)*ePow[2] + (-9640353 + 7434560*eta)*ePow[4] + 
        1440*(-129 + 91*eta)*ePow[6]))/sq1mesqPow[6] + 604800*schi*up7*sq1mesqPow[6] + 
        (15*ePow[1]*(dchi + schi*ePow[1])*(-85608 + 228704*eta + 12*(-64916 + 54271*eta)*ePow[2] + 
        (-125361 + 93184*eta)*ePow[4]))/sq1mesqPow[7];
    eOsc_PN7 /= 37800.;
    *ret_eOsc = eta * dc->vPow[5] * (eOsc_PN5 + dc->vPow[2]*eOsc_PN7);

    // v
    REAL8 vOsc_PN5, vOsc_PN7;
    vOsc_PN5 = 210*schi*up4*ePow[1] + (2*schi*up3*ePow[1]*(49 + 36*ePow[2]))/sq1mesqPow[2] + 
        (schi*up2*ePow[1]*(314 + 111*ePow[2]))/sq1mesqPow[4] + (schi*up*ePow[1]*(602 + 673*ePow[2]))/
        sq1mesqPow[6] + (3*(dchi + schi*ePow[1])*(96 + 292*ePow[2] + 37*ePow[4]))/sq1mesqPow[7];
    vOsc_PN5 /= 45.;

    vOsc_PN7 = -504*schi*up5*ePow[1]*(-339 + 1000*eta + 30*(-1 + 2*eta)*ePow[2]) - 
        (42*schi*up4*ePow[1]*(2322 - 3700*eta + 15*(-299 + 688*eta)*ePow[2]))/sq1mesqPow[2] + 
        840*(137 + 140*eta)*schi*up6*ePow[1]*sq1mesqPow[2] - 
        (2*schi*up3*ePow[1]*(2*(-63479 + 84910*eta) + (-192229 + 119280*eta)*ePow[2] + 
        720*(-15 + 49*eta)*ePow[4]))/sq1mesqPow[4] - 201600*schi*up7*ePow[1]*sq1mesqPow[4] - 
        (schi*up2*ePow[1]*(4*(-45869 + 94990*eta) + 4*(-322676 + 279405*eta)*ePow[2] + 
        15*(-11717 + 8288*eta)*ePow[4]))/sq1mesqPow[6] - 
        (schi*up*ePow[1]*(-247076 + 601720*eta + 4*(-789717 + 782915*eta)*ePow[2] + 
        (-1543861 + 1132320*eta)*ePow[4]))/sq1mesqPow[8] - 
        (15*(dchi + schi*ePow[1])*(16*(-265 + 924*eta) + 24*(-5699 + 6650*eta)*ePow[2] + 
        14*(-12661 + 10122*eta)*ePow[4] + (-11717 + 8288*eta)*ePow[6]))/sq1mesqPow[9];
    vOsc_PN7 /= 12600;
    *ret_vOsc = eta * dc->vPow[6] * (vOsc_PN5 + dc->vPow[2]*vOsc_PN7);

    // n
    REAL8 nOsc_PN5, nOsc_PN7;
    nOsc_PN5 = 210*schi*up4*ePow[1] + (2*schi*up3*ePow[1]*(49 + 36*ePow[2]))/sq1mesqPow[2] + 
        (schi*up2*ePow[1]*(314 + 111*ePow[2]))/sq1mesqPow[4] + (schi*up*ePow[1]*(602 + 673*ePow[2]))/
        sq1mesqPow[6] + (3*(dchi + schi*ePow[1])*(96 + 292*ePow[2] + 37*ePow[4]))/sq1mesqPow[7];
    nOsc_PN5 /= 15.;
    nOsc_PN7 = -504*schi*up5*ePow[1]*(-339 + 1000*eta + 30*(-1 + 2*eta)*ePow[2]) - 
        (42*schi*up4*ePow[1]*(6522 - 3700*eta + 15*(-299 + 688*eta)*ePow[2]))/sq1mesqPow[2] + 
        840*(137 + 140*eta)*schi*up6*ePow[1]*sq1mesqPow[2] + 
        (2*schi*up3*ePow[1]*(85798 - 169820*eta + (161989 - 119280*eta)*ePow[2] - 
        720*(-15 + 49*eta)*ePow[4]))/sq1mesqPow[4] - 201600*schi*up7*ePow[1]*sq1mesqPow[4] - 
        (schi*up2*ePow[1]*(4*(45271 + 94990*eta) + 4*(-324566 + 279405*eta)*ePow[2] + 
        15*(-11717 + 8288*eta)*ePow[4]))/sq1mesqPow[6] - 
        (schi*up*ePow[1]*(520684 + 601720*eta + 4*(-703827 + 782915*eta)*ePow[2] + 
        (-1584181 + 1132320*eta)*ePow[4]))/sq1mesqPow[8] - 
        (15*(dchi + schi*ePow[1])*(16*(1415 + 924*eta) + 24*(-3711 + 6650*eta)*ePow[2] + 
        14*(-12889 + 10122*eta)*ePow[4] + (-11717 + 8288*eta)*ePow[6]))/sq1mesqPow[9];
    nOsc_PN7 /= 4200.;
    *ret_nOsc = eta * dc->vPow[5]*dc->vPow[3] * (nOsc_PN5 + dc->vPow[2]*nOsc_PN7);

    // cl
    REAL8 clOsc_PN5, clOsc_PN7;
    clOsc_PN5 = -9 + 21*up - 35*up2 + (-27 + 39*up + 35*up2)*ePow[2];
    clOsc_PN5 *= 8.*up / (45.*ePow[2]);

    clOsc_PN7 = -180*(570 - 553*eta + 58*(3 + 7*eta)*ePow[2] + 3*(32 + 49*eta)*ePow[4]) + 
        15*up*(13935 - 14770*eta - 6*(-339 + 70*eta)*ePow[2] + (831 + 15190*eta)*ePow[4]) + 
        350*up2*(2*(-621 + 547*eta) + (1269 + 172*eta)*ePow[2] + 9*(-55 + 6*eta)*ePow[4])*sq1mesqPow[2] - 
        315*up3*(-533 + 925*eta + 3*(211 + 85*eta)*ePow[2])*sq1mesqPow[4] + 
        378*up4*(377 + 140*eta + 280*ePow[2])*sq1mesqPow[6] - 88200*up5*sq1mesqPow[8];
    clOsc_PN7 *= -up / (4725.*ePow[2]*sq1mesqPow[2]);
    *ret_clOsc = eta * dc->vPow[5]* (clOsc_PN5 + dc->vPow[2]*clOsc_PN7);

    // clbd
    REAL8 clbdOsc_PN5, clbdOsc_PN7;
    clbdOsc_PN5 = 3*(-12 + 35*up)*up2*ePow[4] + ePow[2]*(144 - 258*up + up2*(92 - 36*sq1mesqPow[1]) + 
        105*up3*(-2 + sq1mesqPow[1])) - up*(18 - 56*up + 105*up2)*(-1 + sq1mesqPow[1]);
    clbdOsc_PN5 *= -2.*up/(45.*ePow[2]);

    clbdOsc_PN7 = 315*up4*ePow[8]*(-6*(55 + 6*eta*(-1 + sq1mesqPow[1]) - 47*sq1mesqPow[1]) - 
        (377 + 140*eta)*up*(-1 + sq1mesqPow[1]) + 1200*up2*(-1 + sq1mesqPow[1])) + 
        up*(90*(444 - 553*eta) + 10*(-10407 + 14770*eta)*up - 525*(-495 + 547*eta)*up2 + 
        252*(-533 + 925*eta)*up3 - 315*(377 + 140*eta)*up4 + 75600*up5)*(-1 + sq1mesqPow[1]) - 
        75600*up6*ePow[10]*(-1 + sq1mesqPow[1]) - 
        3*up2*ePow[6]*(-420*(377 + 140*eta)*up3*(-1 + sq1mesqPow[1]) + 252000*up4*(-1 + sq1mesqPow[1]) + 
        42*up2*(-1409 + 1580*eta*(-1 + sq1mesqPow[1]) + 1049*sq1mesqPow[1]) - 
        175*up*(267 - 327*sq1mesqPow[1] + 17*eta*(-17 + 35*sq1mesqPow[1])) + 
        60*(7*eta*(-7 + 25*sq1mesqPow[1]) - 2*(16 + 61*sq1mesqPow[1]))) + 
        ePow[4]*(378*up4*(241 + 1760*eta*(-1 + sq1mesqPow[1]) - 361*sq1mesqPow[1]) - 
        1890*(377 + 140*eta)*up5*(-1 + sq1mesqPow[1]) + 756000*up6*(-1 + sq1mesqPow[1]) + 
        720*(-87 + 196*eta)*sq1mesqPow[1] + 10*up2*(-825 - 16114*eta - 15639*sq1mesqPow[1] + 
        67270*eta*sq1mesqPow[1]) - 1575*up3*(343 - 383*sq1mesqPow[1] + eta*(-375 + 579*sq1mesqPow[1])) - 
        30*up*(-12*(15 + 496*sq1mesqPow[1]) + 7*eta*(-219 + 2557*sq1mesqPow[1]))) + 
        ePow[2]*(-378*up4*(791 + 1820*eta*(-1 + sq1mesqPow[1]) - 831*sq1mesqPow[1]) + 
        1260*(377 + 140*eta)*up5*(-1 + sq1mesqPow[1]) - 378000*up6*(-1 + sq1mesqPow[1]) - 
        2880*(-27 + 49*eta)*sq1mesqPow[1] + 1575*up3*(419 - 439*sq1mesqPow[1] + 
        eta*(-461 + 563*sq1mesqPow[1])) + 60*up*(576 - 3432*sq1mesqPow[1] + 
        7*eta*(-228 + 1397*sq1mesqPow[1])) - 20*up2*(5079 - 11925*sq1mesqPow[1] + 
        7*eta*(-2143 + 5635*sq1mesqPow[1])));
    clbdOsc_PN7 *= up / (4725.*ePow[2]*sq1mesqPow[3]);
    REAL8 e = ePow[1], z = cos(chi);
    REAL8 z2 = z*z;
    REAL8 eccj = e*(4*z + e*((15*z2)/2. + e*(z*(0.6666666666666666 + (209*z2)/18.) + 
           e*(2.0260416666666665 + z2*(1.3333333333333333 + (389*z2)/24.) + 
              e*(z*(0.3 + z2*(2.175 + (2121*z2)/100.)) + 
                 e*(4.545833333333333 + z2*(0.6125 + z2*(3.16875 + (3187*z2)/120.)) + 
                    e*(z*(0.17857142857142858 + 
                          z2*(1.0172619047619047 + z2*(4.297142857142857 + (15783*z2)/490.))) + 
                       e*(7.310853794642857 + 
                          z2*(0.36830357142857145 + 
                             z2*(1.5053571428571428 + z2*(5.546428571428572 + (10677*z2)/280.))) + 
                          e*(z*(0.12152777777777778 + 
                                z2*(0.617352843915344 + 
                                   z2*(2.0694940476190475 + z2*(6.90547052154195 + (200933*z2)/4536.))
                                   )) + e*
                              (10.22533199249752 + 
                                z2*(0.2521701388888889 + 
                                   z2*(0.9211619543650794 + 
                                    z2*(2.70359623015873 + 
                                    z2*(8.365054563492064 + (1277239*z2)/25200.)))) + 
                                e*(z*(0.08948863636363637 + 
                                    z2*(0.4250710227272727 + 
                                    z2*(1.2758522727272728 + 
                                    z2*(3.402516233766234 + 
                                   z2*(9.917424242424243 + (138603*z2)/2420.))))) + 
                                   e*(13.241405306435825 + 
                                    z2*(0.1864346590909091 + 
                                    z2*(0.6375374842171717 + 
                                    z2*(1.6780791546416547 + 
                                    z2*(4.1618337391774896 + 
                                    z2*(11.555955988455988 + (10653061*z2)/166320.))))) + 
                                    e*(z*
                                    (0.0694110576923077 + 
                                    z2*(0.31545700393356646 + 
                                    z2*(0.8872143915112665 + 
                                    z2*(2.1249262889887888 + 
                                    z2*(4.977708980833981 + 
                                    z2*(13.274924948788586 + (47518529*z2)/669240.)))))) + 
                                    e*(16.33153368514043 + 
                                    z2*(0.14501953125 + 
                                    z2*(0.47480912642045453 + 
                                    z2*(1.1720229640151516 + 
                                    z2*(2.6138266030844157 + 
                                    z2*(5.846773538961039 + 
                                   z2*(15.06933170995671 + (673697*z2)/8624.)))))) + 
                                    e*(z*
                                    (0.055859375 + 
                                    z2*(0.24605118189102565 + 
                                    z2*(0.662933784965035 + 
                                    z2*(1.490113619019869 + 
                                    z2*(3.142501788335122 + 
                                    z2*(6.766047210365392 + 
                                    z2*(16.93477121169429 + (115387637*z2)/1.35135e6))))))) + 
                                    e*(19.478346392855933 + 
                                    z2*(0.11695556640625 + 
                                    z2*(0.3712975135216346 + 
                                    z2*(0.8784329927884615 + 
                                    z2*(1.8398291432959402 + 
                                    z2*(3.708915216727717 + 
                                    z2*(7.73287355006105 + 
                                    z2*(18.867331676260246 + (24320783*z2)/262080.))))))) + 
                                    e*(z*
                                    (0.04620720358455882 + 
                                    z2*(0.19884499942555148 + 
                                    z2*(0.5196653505903563 + 
                                    z2*(1.1200441246483297 + 
                                    z2*(2.219675478091149 + 
                                    z2*(4.3112347002930225 + 
                                    z2*(8.74486908029753 + 
                                    z2*(20.863515381677146 + (2322351459*z2)/2.314312e7)))))))) + 
                                    e*(22.67014560823999 + 
                                    z2*(0.09690677418428309 + 
                                    z2*(0.30065059287875306 + 
                                    z2*(0.690160601938592 + 
                                    z2*(1.3866206593986208 + 
                                    z2*(2.628298225020996 + 
                                    z2*(4.947802625019996 + 
                                    z2*(9.799882700020568 + 
                                    z2*(22.920175650023996 + (1134467353*z2)/1.050192e7)))))))) + 
                                    e*(z*
                                    (0.039046438116776314 + 
                                    z2*(0.16502411864982425 + 
                                    z2*(0.42157002577963754 + 
                                    z2*(0.8818749940749959 + 
                                    z2*(1.677116472245451 + 
                                    z2*(3.06446352274839 + 
                                    z2*(5.617111703317124 + 
                                    z2*(10.895962838116477 + 
                                    z2*(25.0344660840586 + (17076667765*z2)/1.47435288e8))))))))) + 
                                    e*(25.898655734039824 + 
                                    z2*(0.08199752004523027 + 
                                    z2*(0.24989889280714855 + 
                                    z2*(0.5608662797201529 + 
                                    z2*(1.0939751156513735 + 
                                    z2*(1.9905728850230169 + 
                                    z2*(3.5270423125213117 + 
                                    z2*(6.317785071449449 + 
                                    z2*(12.03133057502192 + 
                                    z2*(27.20379920002598 + (182896377*z2)/1.478048e6)))))))))))))))))
                                    ))))))))))));
    clbdOsc_PN7 += 48.*eccj/(5.*sq1mesqPow[2]);
    *ret_clbdOsc = eta * dc->vPow[5]* (clbdOsc_PN5 + dc->vPow[2]*clbdOsc_PN7);
    return X_SUCCESS;
}


static REAL8 evaluate_KeplerEQCoeff(INT p, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe)
{
    if (p==0) return 0.0;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 *ePow = dc->ePow, *sq1m2Pow = dc->sq1mesqPow;
    REAL8 PN0, PN4, PN6;
#define EVALEJ(pp, qq, aa, bb) evaluate_Jpqab((pp), (qq), (aa), (bb), pe)
#define SQ1M2(n) sq1m2Pow[n]
#define EPOW(n) ePow[n]

    PN0 = (2.*EVALEJ(p,p,0,0))/p;
    PN4 = 24*(-5 + 2*eta)*EPOW(1)*EVALEJ(-1 + p,p,0,0) - 
        12*(-5 + 2*eta)*(-1 + p)*EPOW(1)*EVALEJ(-1 + p,p,0,1) + (120 - 48*eta)*EVALEJ(p,p,0,0) + 
        (120 - 48*eta)*EVALEJ(p,p,0,1) + 12*(-5 + 2*eta)*(1 + p)*EPOW(1)*EVALEJ(1 + p,p,0,1) + 
        (60 - 39*eta + eta2)*EPOW(1)*EVALEJ(-1 + p,p,1,0)*SQ1M2(1) - 
        (60 - 39*eta + eta2)*EPOW(1)*EVALEJ(1 + p,p,1,0)*SQ1M2(1) - 
        6*(-5 + 2*eta)*p*EVALEJ(-2 + p,p,0,1)*(-1 + SQ1M2(2)) + 
        6*(-5 + 2*eta)*p*EVALEJ(2 + p,p,0,1)*(-1 + SQ1M2(2));
    PN4 *= 1. / (8.*sq1m2Pow[1]);
    PN6 = 70*(-1 + p)*EPOW(1)*EVALEJ(-1 + p,p,0,1)*
        (5760 + 2016*eta2 + eta*(-13664 + 123*CST_PISQ) - 96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        140*EVALEJ(p,p,0,1)*(5760 + 2016*eta2 + eta*(-13664 + 123*CST_PISQ) - 
        96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        35*p*EVALEJ(2 + p,p,0,1)*(1 - SQ1M2(2))*
        (5760 + 2016*eta2 + eta*(-13664 + 123*CST_PISQ) - 96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        140*EPOW(1)*EVALEJ(-1 + p,p,0,0)*(-5760 - 2016*eta2 + eta*(13664 - 123*CST_PISQ) + 
        96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        70*(1 + p)*EPOW(1)*EVALEJ(1 + p,p,0,1)*
        (-5760 - 2016*eta2 + eta*(13664 - 123*CST_PISQ) + 96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        35*p*EVALEJ(-2 + p,p,0,1)*(1 - SQ1M2(2))*
        (-5760 - 2016*eta2 + eta*(13664 - 123*CST_PISQ) + 96*(30 - 29*eta + 11*eta2)*SQ1M2(2)) + 
        4*EPOW(1)*EVALEJ(-1 + p,p,1,0)*SQ1M2(1)*
        (-33600 + 36960*eta2 + eta*(-147464 + 4305*CST_PISQ) + 
        70*(-720 + 688*eta - 347*eta2 + 13*eta3)*SQ1M2(2)) - 
        4*EPOW(1)*EVALEJ(1 + p,p,1,0)*SQ1M2(1)*
        (-33600 + 36960*eta2 + eta*(-147464 + 4305*CST_PISQ) + 
        70*(-720 + 688*eta - 347*eta2 + 13*eta3)*SQ1M2(2)) + 
        140*EVALEJ(p,p,0,0)*(5760 + 2016*eta2 + eta*(-13664 + 123*CST_PISQ) - 
        96*(30 - 29*eta + 11*eta2)*SQ1M2(2) + 8*eta*(-505 + 131*eta + eta2)*p*SQ1M2(3)) - 
        280*eta*(23 - 73*eta + 13*eta2)*p*EVALEJ(p,p,1,0)*SQ1M2(5);
    PN6 *= 1./(13440.*sq1m2Pow[3]);
    // print_debug("PN0 = %.16e, PN4 = %.16e, PN6 = %.16e\n", PN0, PN4, PN6);
#undef SQ1M2
#undef EPOW
#undef EVALEJ
    return PN0 + dc->vPow[4]*(PN4 + dc->vPow[2]*PN6);
}

REAL8 evaluate_chi_from_l(REAL8 l, REAL8 eta, DynParamsCache *dc, PNEllipticEvaluator *pe)
{
    REAL8 ret = 0.0;
    size_t min_iter = 5;
    size_t cons_iter = 3, ici = 0;
    REAL8 eps = 1e-16;
    REAL8 cum = 0.0;
    size_t p;
    for (p = 1; ; p++)
    {
        REAL8 coeff = evaluate_KeplerEQCoeff(p, eta, dc, pe) ;
        ret += coeff * sin(p*l);
        cum += fabs(coeff);
        if (p > min_iter && coeff < eps*cum) {
            if (ici >cons_iter)
                break;
            else
                ici++;
        } else 
            ici = 0;
    }
    return l + ret;
}

REAL8 evaluate_l_from_chi(REAL8 chi, REAL8 eta, DynParamsCache *dc)
{
    REAL8 *ePow = dc->ePow, *sq1mesqPow = dc->sq1mesqPow;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 schi = sin(chi), s2chi = sin(2.*chi), s3chi = sin(3.*chi);
    REAL8 dchi = evaluate_dchi(chi, eval_beta(ePow[1]));
    REAL8 up = 1./(1. - ePow[1]*cos(chi));
    REAL8 up2 = up*up;
    REAL8 up3 = up2*up;
    REAL8 PN0 = chi - ePow[1]*schi;
    REAL8 PN4 = -((-15 + eta)*eta*schi*up*ePow[1]) + (12*dchi*(5 - 2*eta))/sq1mesqPow[1];
    PN4 /= 8.0;
    REAL8 PN6 = 403200*dchi - 1523200*dchi*eta + 134400*dchi*eta2 + 17220*dchi*eta*CST_PISQ + 403200*dchi*ePow[2] - 
        389760*dchi*eta*ePow[2] + 147840*dchi*eta2*ePow[2] + 67200*s3chi*up3*ePow[3]*sq1mesqPow[1] - 
        51152*eta*s3chi*up3*ePow[3]*sq1mesqPow[1] - 11900*eta2*s3chi*up3*ePow[3]*sq1mesqPow[1] - 
        1820*eta3*s3chi*up3*ePow[3]*sq1mesqPow[1] - 4305*eta*CST_PISQ*s3chi*up3*ePow[3]*sq1mesqPow[1] + 
        schi*up*ePow[1]*(40320*(20 - 13*eta + 2*eta2) + 
        up2*(-4*(-67200 - 93468*eta + 58800*eta2 + 280*eta3 + 4305*eta*CST_PISQ) - 
        3*(-22400 + 215664*eta - 56700*eta2 + 2660*eta3 + 1435*eta*CST_PISQ)*ePow[2] + 
        140*eta*(84 - 375*eta + 65*eta2)*ePow[4]))*sq1mesqPow[1] - 
        1120*eta*s3chi*up3*ePow[5]*sq1mesqPow[1] - 11620*eta2*s3chi*up3*ePow[5]*sq1mesqPow[1] + 
        1820*eta3*s3chi*up3*ePow[5]*sq1mesqPow[1] - 4*s2chi*up3*ePow[2]*sq1mesqPow[1]*
        (67200 + 3360*eta2*(-9 + 2*ePow[2]) - 3*eta*(-6516 + 1435*CST_PISQ + 23940*ePow[2]) - 
        1960*eta3*sq1mesqPow[2]);
    PN6 /= 26880.*sq1mesqPow[3];
    // print_debug("PN0 = %.16e, PN4 = %.16e, PN6 = %.16e\n", PN0, PN4, PN6);
    return PN0 + dc->vPow[4]*(PN4 + dc->vPow[2]*PN6);
}

REAL8 evaluate_k(REAL8 eta, DynParamsCache *dc)
{
    REAL8 eta2 = eta*eta;
    REAL8 *ePow = dc->ePow, *sq1mesqPow = dc->sq1mesqPow;
    REAL8 PN2 = 3/sq1mesqPow[2];
    REAL8 PN4 = (54 - 28*eta + (51 - 26*eta)*ePow[2])/(4*sq1mesqPow[4]);
    REAL8 PN6 = (896*eta2 + ePow[2]*(5120*eta2 + eta*(-22848 + 123*CST_PISQ - 768/sq1mesqPow[1]) + 
        96*(191 + 20/sq1mesqPow[1])) + eta*(-20000 + 492*CST_PISQ - 768/sq1mesqPow[1]) + 
        16*ePow[4]*(156 + 65*eta2 + eta*(-110 + 96/sq1mesqPow[1]) - 240/sq1mesqPow[1]) + 
        960*(7 + 2/sq1mesqPow[1]))/(128*sq1mesqPow[6]);
    return dc->vPow[2] * (PN2 + dc->vPow[2]*(PN4 + dc->vPow[2]*PN6));
}


INT evaluate_PACorr_from_tdynamics(REAL8 eta, REAL8 chi0, 
    PNGW3Dynamics *dyn)
{
    if (dyn->array_extended) return X_FAILURE;
    size_t length = dyn->length;
    dyn->array_extended = CreateREAL8ArrayL(2, PNGW3_NDYNAMICS+1, length);
    dyn->chiVec = dyn->array_extended->data;
    dyn->vPAVec = dyn->chiVec + length;
    dyn->ePAVec = dyn->vPAVec + length;
    dyn->nPAVec = dyn->ePAVec + length;
    dyn->lPAVec = dyn->nPAVec + length;
    dyn->lbdPAVec = dyn->lPAVec + length;
    DynParamsCache dc;
    evaluate_DynParamsCache(dyn->vVec[0], dyn->eVec[0], &dc);

    REAL8 l0 = evaluate_l_from_chi(chi0, eta, &dc);
    REAL8 lbd0, k0;
    k0 = evaluate_k(eta, &dc);
    lbd0 = (1. + k0)*l0;
    REAL8 hlmr, hlmi;
    PNEllipticEvaluator *pe = NULL;
    for (size_t i = 0; i<length; i++) {
        evaluate_DynParamsCache(dyn->vVec[i], dyn->eVec[i], &dc);
        REAL8 l = dyn->lVec[i] + l0;
        dyn->lVec[i] = l;
        dyn->lbdVec[i] += lbd0;
        pe = CreatePNEllipticEvaluator(dyn->eVec[i], 1e-16, 1e-16);
        REAL8 chi = chi0 + evaluate_chi_from_l(l, eta, &dc, pe);
        STRUCTFREE(pe, PNEllipticEvaluator);
        dyn->chiVec[i] = chi;
        evaluate_PACorrections(eta, chi, &dc, 
            dyn->ePAVec+i, dyn->vPAVec+i, dyn->nPAVec+i,
            dyn->lPAVec+i, dyn->lbdPAVec+i);
    }
    return X_SUCCESS;
}

INT evaluate_hSphericalModes_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSphericalModes_TD **ret_hlms,
        REAL8 eps_h, REAL8 eps_pe)
{
    if (!dyn->array_extended) return X_FAILURE;
    size_t length = dyn->length;
    (*ret_hlms) = CreatePNGW3HSphericalModes_TD(length);
    REAL8 hlmr, hlmi;
    PNEllipticEvaluator *pe = NULL;
    DynParamsCache dc;
    // REAL8 eps_h = 1e-5, epspe = 1e-10;
    REAL8 hPref = 4.*eta*sqrt(CST_PI/5.);
    for (size_t i = 0; i<length; i++) {
        (*ret_hlms)->tVec[i] = dyn->tVec[i];
        REAL8 ei = dyn->eVec[i] + dyn->ePAVec[i];
        if (ei < 0) ei = 0.0;
        REAL8 vi = dyn->vVec[i] + dyn->vPAVec[i];
        REAL8 lbdi = dyn->lbdVec[i] + dyn->lbdPAVec[i];
        REAL8 li = dyn->lVec[i] + dyn->lPAVec[i];
        print_err("%.1f[%d/%d]:e = %.5e, v = %.5e\r", 
            100.*i / (REAL8)length, i, length,
            ei, vi);
        pe = CreatePNEllipticEvaluator(ei, eps_pe, eps_pe);
        evaluate_DynParamsCache(vi, ei, &dc);
        COMPLEX16 hlm;
#define EVALHLM(L,M) evaluate_h##L##M(eta, li, &dc, pe, &hlmr, &hlmi, eps_h);\
hlm = (hlmr + I*hlmi) * cexp((M##.)*I*(li - lbdi));\
(*ret_hlms)->h##L##M##rVec[i] = dc.vPow[2] * hPref * creal(hlm);\
(*ret_hlms)->h##L##M##iVec[i] = dc.vPow[2] * hPref * cimag(hlm);

        // L=2
        EVALHLM(2,2)
        EVALHLM(2,1)
        EVALHLM(2,0)

        // L=3
        EVALHLM(3,3)
        EVALHLM(3,2)
        EVALHLM(3,1)
        EVALHLM(3,0)

        // L=4
        EVALHLM(4,4)
        EVALHLM(4,3)
        EVALHLM(4,2)
        EVALHLM(4,1)
        EVALHLM(4,0)

        // L=5
        EVALHLM(5,5)
        EVALHLM(5,4)
        EVALHLM(5,3)
        EVALHLM(5,2)
        EVALHLM(5,1)
        EVALHLM(5,0)

        STRUCTFREE(pe, PNEllipticEvaluator);

    }
#undef EVALHLM
    return X_SUCCESS;
}

INT evaluate_hSphericalModes_pce10_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSphericalModes_TD **ret_hlms)
{
    if (!dyn->array_extended) return X_FAILURE;
    size_t length = dyn->length;
    (*ret_hlms) = CreatePNGW3HSphericalModes_TD(length);
    REAL8 hlmr, hlmi;
    DynParamsCache dc;
    // REAL8 eps_h = 1e-5, epspe = 1e-10;
    REAL8 hPref = 4.*eta*sqrt(CST_PI/5.);
    for (size_t i = 0; i<length; i++) {
        (*ret_hlms)->tVec[i] = dyn->tVec[i];
        REAL8 ei = dyn->eVec[i] + dyn->ePAVec[i];
        if (ei < 0) ei = 0.0;
        REAL8 vi = dyn->vVec[i] + dyn->vPAVec[i];
        REAL8 lbdi = dyn->lbdVec[i] + dyn->lbdPAVec[i];
        REAL8 li = dyn->lVec[i] + dyn->lPAVec[i];
        print_err("%.1f[%d/%d]:e = %.5e, v = %.5e\r", 
            100.*i / (REAL8)length, i, length,
            ei, vi);
        evaluate_DynParamsCache(vi, ei, &dc);
        COMPLEX16 hlm;
#define EVALHLM(L,M) hlm = calculate_h##L##M##_pc_e10(eta, li, &dc) * cexp((M##.)*I*(li - lbdi));\
(*ret_hlms)->h##L##M##rVec[i] = dc.vPow[2] * hPref * creal(hlm);\
(*ret_hlms)->h##L##M##iVec[i] = dc.vPow[2] * hPref * cimag(hlm);

        // L=2
        EVALHLM(2,2)
        EVALHLM(2,1)
        EVALHLM(2,0)

        // L=3
        EVALHLM(3,3)
        EVALHLM(3,2)
        EVALHLM(3,1)
        EVALHLM(3,0)

        // L=4
        EVALHLM(4,4)
        EVALHLM(4,3)
        EVALHLM(4,2)
        EVALHLM(4,1)
        EVALHLM(4,0)

        // L=5
        EVALHLM(5,5)
        EVALHLM(5,4)
        EVALHLM(5,3)
        EVALHLM(5,2)
        EVALHLM(5,1)
        EVALHLM(5,0)


    }
#undef EVALHLM
    return X_SUCCESS;
}


INT evaluate_h22_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSingleSphericalMode_TD **ret_h22, REAL8 eps_h, REAL8 eps_pe)
{
    if (!dyn->array_extended) return X_FAILURE;
    size_t length = dyn->length;
    (*ret_h22) = CreatePNGW3HSingleSphericalMode_TD(length);
    REAL8 hlmr, hlmi;
    PNEllipticEvaluator *pe = NULL;
    DynParamsCache dc;
    REAL8 hPref = 4.*eta*sqrt(CST_PI/5.);
    for (size_t i = 0; i<length; i++) {
        (*ret_h22)->tVec[i] = dyn->tVec[i];
        REAL8 ei = dyn->eVec[i] + dyn->ePAVec[i];
        if (ei < 0) ei = 0.0;
        REAL8 vi = dyn->vVec[i] + dyn->vPAVec[i];
        REAL8 lbdi = dyn->lbdVec[i] + dyn->lbdPAVec[i];
        REAL8 li = dyn->lVec[i] + dyn->lPAVec[i];
        print_err("%.1f[%d/%d]:e = %.5e, v = %.5e\r", 
            100.*i / (REAL8)length, i, length,
            ei, vi);
        pe = CreatePNEllipticEvaluator(ei, eps_pe, eps_pe);
        evaluate_DynParamsCache(vi, ei, &dc);
        COMPLEX16 hlm;

        evaluate_h22(eta, li, &dc, pe, &hlmr, &hlmi, eps_h);
        hlm = (hlmr + I*hlmi) * cexp((2.)*I*(li - lbdi));\
        (*ret_h22)->hrVec[i] = dc.vPow[2] * hPref * creal(hlm);\
        (*ret_h22)->hiVec[i] = dc.vPow[2] * hPref * cimag(hlm);

        STRUCTFREE(pe, PNEllipticEvaluator);

    }
    return X_SUCCESS;
}

INT evaluate_h22_pce10_from_PADynamics(REAL8 eta, PNGW3Dynamics *dyn, 
        PNGW3HSingleSphericalMode_TD **ret_h22, REAL8 eps_h, REAL8 eps_pe)
{
    if (!dyn->array_extended) return X_FAILURE;
    size_t length = dyn->length;
    (*ret_h22) = CreatePNGW3HSingleSphericalMode_TD(length);
    REAL8 hlmr, hlmi;
    PNEllipticEvaluator *pe = NULL;
    DynParamsCache dc;
    REAL8 hPref = 4.*eta*sqrt(CST_PI/5.);
    for (size_t i = 0; i<length; i++) {
        (*ret_h22)->tVec[i] = dyn->tVec[i];
        REAL8 ei = dyn->eVec[i] + dyn->ePAVec[i];
        if (ei < 0) ei = 0.0;
        REAL8 vi = dyn->vVec[i] + dyn->vPAVec[i];
        REAL8 lbdi = dyn->lbdVec[i] + dyn->lbdPAVec[i];
        REAL8 li = dyn->lVec[i] + dyn->lPAVec[i];
        print_err("%.1f[%d/%d]:e = %.5e, v = %.5e\r", 
            100.*i / (REAL8)length, i, length,
            ei, vi);
        pe = CreatePNEllipticEvaluator(ei, eps_pe, eps_pe);
        evaluate_DynParamsCache(vi, ei, &dc);
        COMPLEX16 hlm = calculate_h22_pc_e10(eta, li, &dc)* cexp((2.)*I*(li - lbdi));

        (*ret_h22)->hrVec[i] = dc.vPow[2] * hPref * creal(hlm);
        (*ret_h22)->hiVec[i] = dc.vPow[2] * hPref * cimag(hlm);
        
        STRUCTFREE(pe, PNEllipticEvaluator);

    }
    return X_SUCCESS;
}
