/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#include "pngw3_evolve.h"
#include "pngw3_ehfunc.h"
#include <gsl/gsl_errno.h>



void evaluate_DynParamsCache(REAL8 v, REAL8 e, DynParamsCache *cache)
{
    if (!cache || e<0.0) return;
    cache->ePow[0] = 1.0;
    cache->sq1mesqPow[0] = 1.0;
    cache->betaPow[0] = 1.0;
    cache->beta2p1Pow[0] = 1.0;
    cache->beta2m1Pow[0] = 1.0;
    cache->vPow[0] = 1.0;
    REAL8 sq1mesq = sqrt(1 - e*e);
    REAL8 beta = e < 1e-16 ? e : (1. - sq1mesq) / e;
    REAL8 beta2 = beta*beta;
    REAL8 beta2p1 = 1. + beta2;
    REAL8 beta2m1 = -1. + beta2;
    for (int i=1; i<MAXPOW_EPOW+1; i++ ) cache->ePow[i] = cache->ePow[i-1]*e;
    for (int i=1; i<MAXPOW_SQ1MESQ+1; i++ ) cache->sq1mesqPow[i] = cache->sq1mesqPow[i-1]*sq1mesq;
    for (int i=1; i<MAXPOW_BETA+1; i++ ) cache->betaPow[i] = cache->betaPow[i-1]*beta;
    for (int i=1; i<MAXPOW_BETA2P1+1; i++ ) cache->beta2p1Pow[i] = cache->beta2p1Pow[i-1]*beta2p1;
    for (int i=1; i<MAXPOW_BETA2M1+1; i++ ) cache->beta2m1Pow[i] = cache->beta2m1Pow[i-1]*beta2m1;
    for (int i=1; i<MAXPOW_V+1; i++ ) cache->vPow[i] = cache->vPow[i-1]*v;
    return;
}

static REAL8 evaluate_n_from_DynParamsCache(REAL8 eta, DynParamsCache *cache)
{
    REAL8 ret_n = 0.0;
    REAL8 PN0, PN2, PN4, PN6;
    REAL8 *vPow = cache->vPow;
    REAL8 *betaPow = cache->betaPow;
    REAL8 *beta2p1Pow = cache->beta2p1Pow;
    REAL8 *beta2m1Pow = cache->beta2m1Pow;
    REAL8 eta2 = eta*eta;
    REAL8 pisq123 = 123.*CST_PISQ;
    PN0 = 1.0;
    PN2 = -3.*beta2p1Pow[2] / beta2m1Pow[2];
    PN4 = 0.5 * beta2p1Pow[2] * ( 40*betaPow[2] * (-3. + 2.*eta) + (1. + betaPow[4])*(-9. + 14.*eta) ) / beta2m1Pow[4];

    PN6 = -0.3125*beta2p1Pow[2] * (432. + (-3848 + pisq123)*eta + 224.*eta2 + 
        betaPow[2]*(13152. + 5.*(-6880+pisq123)*eta + 6016.*eta2) + 
        betaPow[4]*8.*(3348. + (-8086.+pisq123)*eta + 1968.*eta2) + 
        betaPow[6]*(3552. + 5.*(-6112.+pisq123)*eta + 6016.*eta2) +
        betaPow[8]*(-528. + (-3464.+pisq123)*eta + 224.*eta2) ) / beta2m1Pow[6];
    // PN8 = ...
    // print_debug("PN0 = %.16e, PN2 = %.16e\n\tPN4 = %.16e, PN6 = %.16e\n", 
    //         PN0, PN2, PN4, PN6);
    return vPow[3]*(PN0 + vPow[2]*(PN2 + vPow[2]*(PN4 + vPow[2]*PN6)));
}

static REAL8 evaluate_dvdt_from_DynParamsCache(REAL8 eta, DynParamsCache *cache)
{
    REAL8 *vPow = cache->vPow;
    REAL8 *ePow = cache->ePow;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 *sq1m2Pow = cache->sq1mesqPow;
    REAL8 dm2 = 1. - 4.*eta;
    REAL8 feh_varphi, feh_varphiT;
    REAL8 feh_varphi01, feh_varphi10;
    REAL8 feh_beta, feh_gamma;
    REAL8 feh_F, feh_chi;
    REAL8 PN0, PN2, PN4, PN6;
    REAL8 PN3Tail, PN5Tail, PN6Tail;

    feh_varphi = evaluate_ehfunc_varphi(ePow[2]);
    feh_varphiT = evaluate_ehfunc_varphiTilde(ePow[2]);
    feh_beta = evaluate_ehfunc_beta(ePow[2]);
    feh_gamma = evaluate_ehfunc_gamma(ePow[2]);
    feh_varphi01 = evaluate_ehfunc_varphi01(eta, ePow[2]);
    feh_varphi10 = evaluate_ehfunc_varphi10(ePow[2]);
    feh_F = evaluate_ehfunc_F(ePow[2]);
    feh_chi = evaluate_ehfunc_chi(ePow[2]);

    REAL8 pnv2 = vPow[2] / sq1m2Pow[2];
    // REAL8 pnv2 = pnv*pnv;

    PN0 = (96.+292.*ePow[2]+37.*ePow[4])/(15.);
    PN2 = -1./30.*(eta*(528. + 5700.*ePow[2] + 5061.*ePow[4] + 296.*ePow[6])) + 
        (-11888. + 87720.*ePow[2] + 171038.*ePow[4] + 11717.*ePow[6])/840.;
    PN4 = (sq1m2Pow[1]*(1451520 + 64532160*ePow[2] + 66316320*ePow[4] + 2646000*ePow[6]))/90720 + 
        eta*((sq1m2Pow[1]*(-580608 - 25812864*ePow[2] - 26526528*ePow[4] - 1058400*ePow[6]))/90720 + 
        (4514976 + 15464736*ePow[2] - 207204264*ePow[4] - 123108426*ePow[6] - 3259980*ePow[8])/90720) + 
        (eta2*(1903104 + 61282032*ePow[2] + 166506060*ePow[4] + 64828848*ePow[6] + 1964256*ePow[8]))/
        90720 + (-360224 - 92846560*ePow[2] + 783768*ePow[4] + 83424402*ePow[6] + 3523113*ePow[8])/90720;
    REAL8 lgterm = log(0.5*vPow[2]*(1. + sq1m2Pow[1])/(CONST_X0 *(1-ePow[2]) ));
    PN6 = 614389219./445500. + (19898670811.*ePow[2])/2079000. + (8036811073.*ePow[4])/24948000. + 
        (985878037.*ePow[6])/16632000. + (2814019181.*ePow[8])/1056000. + 
        lgterm*(54784./525. + (465664.*ePow[2])/315. + (4426376.*ePow[4])/1575. + (1498856.*ePow[6])/1575. + 
        (10593.*ePow[8])/350.) + eta*(-57265081./34020. + (123.*CST_PISQ)/2. + (2678401319.*ePow[2])/340200. + 
        (3239.*CST_PISQ*ePow[2])/48. + (43741211273.*ePow[4])/1360800. - (197087.*CST_PISQ*ePow[4])/960. + 
        (54136669.*ePow[6])/43200. - (261211.*CST_PISQ*ePow[6])/1920. - (4342403.*ePow[8])/1008. - 
        (4059.*CST_PISQ*ePow[8])/640. - (1874543.*ePow[10])/30240.) + 
        eta3*(-1121./81. - (1287385.*ePow[2])/972. - (33769597.*ePow[4])/3888. - (3200965.*ePow[6])/324. - 
        (982645.*ePow[8])/486. - (8288.*ePow[10])/243.) + (33332681.*ePow[10])/591360. + 
        eta2*(-16073./420. - (9657701.*ePow[2])/2520. + (1306589.*ePow[4])/2016. + (62368205.*ePow[6])/4032. + 
        (3542389.*ePow[8])/840. + (109733.*ePow[10])/2520.) + 
        (-1425319./3375. + (933454.*ePow[2])/1125. + (840635951.*ePow[4])/63000. + (702667207.*ePow[6])/94500. + 
        eta*(9874./315. - (41.*CST_PISQ)/30. - (2257181.*ePow[2])/189. + (45961.*CST_PISQ*ePow[2])/720. - 
        (4927789.*ePow[4])/180. + (6191.*CST_PISQ*ePow[4])/96. - (6830419.*ePow[6])/756. + 
        (287.*CST_PISQ*ePow[6])/2880. - (427733.*ePow[8])/2520.) + (18801.*ePow[8])/112. + 
        eta2*(632./15. + (125278.*ePow[2])/45. + (317273.*ePow[4])/45. + (232177.*ePow[6])/90. + 
        (4739.*ePow[8])/90.))*sq1m2Pow[1];

    PN3Tail = 128.*CST_PI*feh_varphi/5.;
    PN5Tail = (16403.*dm2*feh_beta + 112.*(dm2*feh_gamma + 12.*(-129.-2*eta+ePow[2]*(15.+2.*eta))*feh_varphi / sq1m2Pow[2] + 
        72.*(feh_varphi01 + 21.*feh_varphi10/sq1m2Pow[2]  + 2.*feh_varphiT/sq1m2Pow[3] ) ) ) * CST_PI / 315.;
    PN6Tail = 32.*(-515063./55125. + 16*CST_PISQ/15. - 856./175. * log(vPow[2]/CONST_X0P))*feh_F - 54784.*feh_chi/525;
    // print_debug("PN0 = %.16e\n", PN0);
    // print_debug("PN2 = %.16e\n", PN2);
    // print_debug("PN4 = %.16e\n", PN4);
    // print_debug("PN6 = %.16e\n", PN6);
    // print_debug("PN3Tail = %.16e\n", PN3Tail);
    // print_debug("PN5Tail = %.16e\n", PN5Tail);

    // print_debug("feh_varphi = %.16e\n", feh_varphi);
    REAL8 inst = vPow[3]*vPow[6]*eta * (PN0 + pnv2*(PN2 + pnv2*(PN4 + pnv2*PN6))) / sq1m2Pow[7];
    REAL8 tail = vPow[6]*vPow[6]*eta * (PN3Tail + vPow[2]*PN5Tail + vPow[3]*PN6Tail);
    // print_debug("tail = %.16e\n", tail);
    // print_debug("inst = %.16e\n", inst);
    return inst + tail;
}

static REAL8 evaluate_dndt_from_DynParamsCache(REAL8 eta, DynParamsCache *cache)
{
    REAL8 *vPow = cache->vPow;
    REAL8 *ePow = cache->ePow;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 *sq1m2Pow = cache->sq1mesqPow;
    REAL8 dm2 = 1. - 4.*eta;
    REAL8 feh_varphi, feh_varphiT;
    REAL8 feh_varphi01, feh_varphi10;
    REAL8 feh_beta, feh_gamma;
    REAL8 feh_F, feh_chi;
    REAL8 PN0, PN2, PN4, PN6;
    REAL8 PN3Tail, PN5Tail, PN6Tail;

    feh_varphi = evaluate_ehfunc_varphi(ePow[2]);
    feh_varphiT = evaluate_ehfunc_varphiTilde(ePow[2]);
    feh_beta = evaluate_ehfunc_beta(ePow[2]);
    feh_gamma = evaluate_ehfunc_gamma(ePow[2]);
    feh_varphi01 = evaluate_ehfunc_varphi01(eta, ePow[2]);
    feh_varphi10 = evaluate_ehfunc_varphi10(ePow[2]);
    feh_F = evaluate_ehfunc_F(ePow[2]);
    feh_chi = evaluate_ehfunc_chi(ePow[2]);

    REAL8 pnv2 = vPow[2] / sq1m2Pow[2];
    // REAL8 pnv2 = pnv*pnv;

    PN0 = (96.+292.*ePow[2]+37.*ePow[4])/(5.);
    PN2 = (-16.*(2423. + 924.*eta) - 24.*(-1667. + 6650.*eta)*ePow[2] + (174230. - 141708.*eta)*ePow[4] + 
        (11717. - 8288.*eta)*ePow[6])/280.;
    PN4 = (-1159. + 686925.*eta + 59472.*eta2)/945. - (4819994.*ePow[2])/945. + 3231.*eta*ePow[2] + 
        (182387.*eta2*ePow[2])/90. - (203957.*ePow[4])/90. - (764357.*eta*ePow[4])/140. + (396443.*eta2*ePow[4])/72. + 
        (4760347.*ePow[6])/1680. - (993011.*eta*ePow[6])/240. + (192943.*eta2*ePow[6])/90. + 
        (391457./3360. - (6037.*eta)/56. + (2923.*eta2)/45.)*ePow[8] + 
        (48. - (96.*eta)/5 - (2134.*(-5. + 2.*eta)*ePow[2])/5. - (2193.*(-5. + 2.*eta)*ePow[4])/5. + (175.*ePow[6])/2. - 
        35.*eta*ePow[6])*sq1m2Pow[1];
    REAL8 lgterm = log(0.5*vPow[2]*(1. + sq1m2Pow[1])/(CONST_X0 *(1-ePow[2]) ));
    PN6 = lgterm*((107.*(3072. + 43520.*ePow[2] + 82736.*ePow[4] + 28016.*ePow[6] + 891.*ePow[8]))/1050.) +
        (4915859933./1039500. + (1463719.*eta)/2268. - (711931.*eta2)/420. - (1121.*eta3)/27. - (369.*eta*CST_PISQ)/10. + 
        (10928916619.*ePow[2])/297000. + (4697941919.*eta*ePow[2])/113400. - (23667137.*eta2*ePow[2])/840. - (1287385.*eta3*ePow[2])/324. - 
        (1599.*eta*CST_PISQ*ePow[2])/80. - (127363208627.*ePow[4])/8316000. + (9286298159.*eta*ePow[4])/64800. - (5331901.*eta2*ePow[4])/224. - 
        (33769597.*eta3*ePow[4])/1296. - (94423.*eta*CST_PISQ*ePow[4])/160. - (82502370763.*ePow[6])/5544000 + (2015302783.*eta*ePow[6])/100800. + 
        (18129215.*eta2*ePow[6])/448. - (3200965.*eta3*ePow[6])/108. - (12751.*eta*CST_PISQ*ePow[6])/32. + (59641969601.*ePow[8])/7392000. - 
        (2779943.*eta*ePow[8])/210. + (5399701.*eta2*ePow[8])/420. - (982645.*eta3*ePow[8])/162. - (12177.*eta*CST_PISQ*ePow[8])/640. + 
        (33332681.*ePow[10])/197120. - (1874543.*eta*ePow[10])/10080. + (109733.*eta2*ePow[10])/840. - (8288.*eta3*ePow[10])/81. + 
        (-2667319./1125. + (56242.*eta)/105. + (632.*eta2)/5. - (41.*eta*CST_PISQ)/10. - (2673296.*ePow[2])/375. - (10074037.*eta*ePow[2])/315. + 
        (125278.*eta2*ePow[2])/15. + (45961.*eta*CST_PISQ*ePow[2])/240. + (700397951.*ePow[4])/21000. - (4767517.*eta*ePow[4])/60. + 
        (317273.*eta2*ePow[4])/15. + (6191.*eta*CST_PISQ*ePow[4])/32. + (708573457.*ePow[6])/31500. - (6849319.*eta*ePow[6])/252. + 
        (232177.*eta2*ePow[6])/30. + (287.*eta*CST_PISQ*ePow[6])/960. + (56403.*ePow[8])/112. - (427733.*eta*ePow[8])/840. + (4739.*eta2*ePow[8])/30.)*sq1m2Pow[1]);

    PN3Tail = 384.*CST_PI*feh_varphi/5.;
    PN5Tail = (16403.*dm2*feh_beta + 112.*(dm2*feh_gamma + (12./sq1m2Pow[2])*((-147.-2*eta+ePow[2]*(15.+2.*eta))*feh_varphi + 6.*((1.-ePow[2])*feh_varphi01 + 21.*feh_varphi10)  ) ) ) * CST_PI / 105.;
    PN6Tail = 32.*((-515063./18375. + 16*CST_PISQ/5. - 2568./175. * log(vPow[2]/CONST_X0P))*feh_F - 1712.*feh_chi/175.);
    // print_debug("PN0 = %.16e\n", PN0);
    // print_debug("PN2 = %.16e\n", PN2);
    // print_debug("PN4 = %.16e\n", PN4);
    // print_debug("PN6 = %.16e\n", PN6);
    // print_debug("PN3Tail = %.16e\n", PN3Tail);
    // print_debug("PN5Tail = %.16e\n", PN5Tail);

    // print_debug("feh_varphi = %.16e\n", feh_varphi);
    REAL8 inst = vPow[5]*vPow[6]*eta * (PN0 + pnv2*(PN2 + pnv2*(PN4 + pnv2*PN6))) / sq1m2Pow[7];
    REAL8 tail = vPow[5]*vPow[6]*vPow[3]*eta * (PN3Tail + vPow[2]*PN5Tail + vPow[3]*PN6Tail);
    // print_debug("tail = %.16e\n", tail);
    // print_debug("inst = %.16e\n", inst);
    return inst + tail;
}


static REAL8 evaluate_dedt_from_DynParamsCache(REAL8 eta, DynParamsCache *cache)
{
    REAL8 *vPow = cache->vPow;
    REAL8 *ePow = cache->ePow;
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 *sq1m2Pow = cache->sq1mesqPow;
    REAL8 dm2 = 1. - 4.*eta;
    REAL8 feh_varphi, feh_varphiT;
    REAL8 feh_varphi01, feh_varphi10;
    REAL8 feh_beta, feh_gamma;
    REAL8 feh_F, feh_chi;
    REAL8 PN0, PN2, PN4, PN6;
    REAL8 PN3Tail, PN5Tail, PN6Tail;

    PN0 = (-304. - 121.*ePow[2])/15.;
    PN2 = (8.*(8451. + 28588.*eta) + 12.*(-59834. + 54271.*eta)*ePow[2] + (-125361. + 93184.*eta)*ePow[4])/2520.;
    PN4 = 949877./1890. + (3082783.*ePow[2])/2520. - (23289859.*ePow[4])/15120. + 
        sq1m2Pow[1]*(-1336./3. - (2321.*ePow[2])/2. - (565.*ePow[4])/6. + 
        eta*(2672./15. + (2321.*ePow[2])/5. + (113.*ePow[4])/3.)) + 
        eta2*(-752./5. - (64433.*ePow[2])/40. - (127411.*ePow[4])/90. - (821.*ePow[6])/9.) - 
        (420727.*ePow[6])/3360. + eta*(-18763./42. + (988423.*ePow[2])/840. + (13018711.*ePow[4])/5040. + 
        (362071.*ePow[6])/2520.);
    REAL8 lgterm = log(0.5*vPow[2]*(1. + sq1m2Pow[1])/(CONST_X0 *(1-ePow[2]) ));
    PN6 = -(107./3150.)*(24608.+89024.*ePow[2]+42884.*ePow[4]+1719.*ePow[6])*lgterm - 730168./(23625.*(1. + sq1m2Pow[1])) + 
        (-54208557619./6237000. - (46226320013.*ePow[2])/6237000. + (116987170177.*ePow[4])/16632000. - 
        (5891934893.*ePow[6])/1232000. + eta2*(4088921./2520. + (21283907.*ePow[2])/3024. - (61093675.*ePow[4])/4032. - 
        (24757667.*ePow[6])/2520. - (41179.*ePow[8])/216.) - (302322169.*ePow[8])/1774080. + 
        eta3*(61001./486. + (86910509.*ePow[2])/19440. + (2223241.*ePow[4])/180. + (11792069.*ePow[6])/2430. + 
        (193396.*ePow[8])/1215.) + eta*(-50099023./113400. - (779.*CST_PISQ)/10. - (28141879.*ePow[2])/900. + 
        (139031.*CST_PISQ*ePow[2])/960. - (11499615139.*ePow[4])/907200. + (271871.*CST_PISQ*ePow[4])/1920. + 
        (5028323.*ePow[6])/560. + (6519.*CST_PISQ*ePow[6])/640. + (1921387.*ePow[8])/10080.) + 
        (22713049./15750. - (89395687.*ePow[2])/7875. - (5321445613.*ePow[4])/378000. + 
        eta2*(-54332./45. - (681989.*ePow[2])/90. - (225106.*ePow[4])/45. - (3197.*ePow[6])/18.) - (186961.*ePow[6])/336. + 
        eta*(5526991./945. - (8323.*CST_PISQ)/180. + (38295557.*ePow[2])/1260. - (94177.*CST_PISQ*ePow[2])/960. + 
        (26478311.*ePow[4])/1512. - (2501.*CST_PISQ*ePow[4])/2880. + (289691.*ePow[6])/504.))*sq1m2Pow[1]);
    
    PN3Tail = -(394./3.)*sq1m2Pow[1]*CST_PI*evaluate_ehfunc_varphiEPN3(ePow[2]);
    PN5Tail = CST_PI*((128./5.)*(sq1m2Pow[1]*evaluate_ehfunc_varphiE01PN5(eta, ePow[2]) + evaluate_ehfunc_varphiE10PN5(eta, ePow[2])) +
        (197./(117.*sq1m2Pow[1]))*(1785.-104.*eta + ePow[2]*(-261.+62.*eta)) * evaluate_ehfunc_varphiEPN3(ePow[2]) - 
        dm2*sq1m2Pow[1]*((-220781./630.)*evaluate_ehfunc_betaEPN5(ePow[2]) - (248./45.)*evaluate_ehfunc_gammaEPN5(ePow[2]) ));
    PN6Tail = sq1m2Pow[1]*( (396083447./165375. - 12304.*CST_PISQ/45. + 658264.*log(vPow[2]/CONST_X0P)/525.)*evaluate_ehfunc_FEPN6(ePow[2]) - 54784.*evaluate_ehfunc_chiEPN6(ePow[2])/525. );
    REAL8 pnv2 = vPow[2] / sq1m2Pow[2];
    REAL8 inst = vPow[3]*vPow[5]* ePow[1] *eta  * (PN0 + pnv2*(PN2 + pnv2*(PN4 + pnv2*PN6))) / sq1m2Pow[5];
    REAL8 tail = vPow[5]*vPow[6]* ePow[1] *eta  * (PN3Tail + vPow[2]*(PN5Tail + vPow[1]*PN6Tail));
    // print_debug("PN0 = %.16e\n", PN0);
    // print_debug("PN2 = %.16e\n", PN2);
    // print_debug("PN4 = %.16e\n", PN4);
    // print_debug("PN6 = %.16e\n", PN6);

    // print_debug("PN3Tail = %.16e\n", PN3Tail);
    // print_debug("PN5Tail = %.16e\n", PN5Tail);
    // print_debug("PN6Tail = %.16e\n", PN6Tail);

    // print_debug("inst = %.16e\n", inst);
    // print_debug("tail = %.16e\n", tail);
    return inst + tail;
}

static INT evaluate_init_conditions(REAL8 eta, REAL8 v0, REAL8 e0, REAL8 *init)
{
    DynParamsCache dyncache;
    evaluate_DynParamsCache(v0, e0, &dyncache);
    REAL8 n0 = evaluate_n_from_DynParamsCache(eta, &dyncache);
    REAL8 l0 = 0.0;
    REAL8 lbd0 = 0.0; // TBD
    //  v, e, n, l, lbd
    init[0] = v0;
    init[1] = e0;
    init[2] = n0;
    init[3] = l0;
    init[4] = lbd0;
    return X_SUCCESS;
}

static int GSL_PNGW3_vEvolve_RK4_func_dydt(
    double t,    
    const REAL8     values[],   /**<< The dynamics v, e, n, l, lbd */
    REAL8           dvalues[],  /**<< The derivatives dr/dt, dphi/dt. dpr/dt and dpphi/dt */
    void            *funcParams /**<< Structure containing all the necessary parameters */
)
{
    DynParamsCache dyncache;
    PNGW3Core *params = (PNGW3Core*)((RK4Params*)funcParams)->params;
    REAL8 eta = params->eta;
    REAL8 v = values[0];
    REAL8 e = values[1];
    REAL8 n = values[2];
    evaluate_DynParamsCache(v, e, &dyncache);
    REAL8 dvdt, dedt, dndt;

    dvdt = evaluate_dvdt_from_DynParamsCache(eta, &dyncache);
    dedt = evaluate_dedt_from_DynParamsCache(eta, &dyncache);
    dndt = evaluate_dndt_from_DynParamsCache(eta, &dyncache);

    dvalues[0] = dvdt;
    dvalues[1] = dedt;
    dvalues[2] = dndt;
    dvalues[3] = n;
    dvalues[4] = dyncache.vPow[3];
    // print_debug("t, v, e, n, l, lbd = %.16e, %.16e, %.16e, %.16e, %.16e, %.16e\n",
    //     t, v, e, n, values[3], values[4]);

    return GSL_SUCCESS;
}

static int GSL_PNGW3_vEvolve_RK4_stop(
    double t,    
    const REAL8     values[],   /**<< The dynamics r, phi, pr, pphi */
    REAL8           dvalues[],  /**<< The derivatives dr/dt, dphi/dt. dpr/dt and dpphi/dt */
    void            *funcParams /**<< Structure containing all the necessary parameters */
)
{
    PNGW3Core *params = (PNGW3Core*)((RK4Params*)funcParams)->params;
    if (values[0] >= params->hParams->vmax) {
        return 1;
    }
    return GSL_SUCCESS;
}



INT evolve_binary_orbit_t(
    PNGW3Core       *pnparams,
    PNGW3Dynamics   **outdyn
)
{
    int is_failed = 0;
    RK4Integrator *integrator = CreateRK4Integrator(PNGW3_NDYNAMICS, 
            GSL_PNGW3_vEvolve_RK4_func_dydt, 
            GSL_PNGW3_vEvolve_RK4_stop,
            pnparams->hParams->EPS_ABS,
            pnparams->hParams->EPS_REL);

    RK4Params rk4Params;
    rk4Params.params = (void *)pnparams;
    rk4Params.cache = NULL;

    REAL8Array *dynarr = NULL;
    // v, e, n, l, lbd
    
    REAL8 init_conds[PNGW3_NDYNAMICS];
    evaluate_init_conditions(pnparams->eta, 
        pnparams->initial_v, pnparams->initial_e,
        init_conds);
    REAL8 dt = 1.0;
    if (run_Adaptive_RK4Integrator(&dynarr, integrator, &rk4Params, 
            init_conds, 0, 100., 
            dt, dt/20., 
            FALSE) != X_SUCCESS) 
        {is_failed = 1; goto EXIT;}
    UINT lendyn = dynarr->dimLength->data[1];
    *outdyn = CreatePNGW3Dynamics(dynarr);
EXIT:
    STRUCTFREE(integrator, RK4Integrator);
    if (is_failed) {
        X_PRINT_INFO("ERROR! Integration Failed!");
        return X_FAILURE;
    }
    return X_SUCCESS;
}

void test_func()
{
    // DynParamsCache cache;
    // REAL8 eta = 1./5.;
    // REAL8 v = 0.2, e = 0.5;
    // evaluate_DynParamsCache(v, e, &cache);
    // print_debug("dvdt = %.16e\n", evaluate_dvdt_from_DynParamsCache(eta, &cache));
    // print_debug("dedt = %.16e\n", evaluate_dedt_from_DynParamsCache(eta, &cache));
    // print_debug("dndt = %.16e\n", evaluate_dndt_from_DynParamsCache(eta, &cache));
    PNGW3HyperParams hParams;
    REAL8 eta = 0.25;
    REAL8 e0 = 0.2;
    REAL8 v0 = 0.2;
    REAL8 chi0 = 0.0;
    REAL8 phic = 0.0;
    REAL8 tc = 0.0;
    REAL8 iota = 0.0;
    hParams.EPS_ABS = 1e-13;
    hParams.EPS_REL = 1e-13;
    hParams.PNOrder = 6;
    hParams.vmax = VISCO;
    PNGW3Core *pnparams = CreatePNGW3Core(eta, e0, v0, chi0, phic, tc, iota, &hParams);
    PNGW3Dynamics *dyn = NULL;
    evolve_binary_orbit_t(pnparams, &dyn);
    for (size_t i=0; i<dyn->length; i++) {
        // t, v, e, n, l, lbd
        print_err("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", 
            dyn->tVec[i], dyn->vVec[i], dyn->eVec[i], 
            dyn->nVec[i], dyn->lVec[i], dyn->lbdVec[i]);
    }
    return;
}
