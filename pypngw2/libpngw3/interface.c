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

typedef struct {
    INT         length;
    REAL8Vector *tVec;
    REAL8Vector *vVec;
    REAL8Vector *eVec;
    REAL8Vector *nVec;
    REAL8Vector *lVec;
    REAL8Vector *lbdVec;

    REAL8Vector *chiVec;
    REAL8Vector *vPAVec;
    REAL8Vector *ePAVec;
    REAL8Vector *nPAVec;
    REAL8Vector *lPAVec;
    REAL8Vector *lbdPAVec;

    // h22
    REAL8Vector *h22r;
    REAL8Vector *h22i;

    // h21
    REAL8Vector *h21r;
    REAL8Vector *h21i;

    // h20
    REAL8Vector *h20r;
    REAL8Vector *h20i;

    // h33
    REAL8Vector *h33r;
    REAL8Vector *h33i;

    // h32
    REAL8Vector *h32r;
    REAL8Vector *h32i;

    // h31
    REAL8Vector *h31r;
    REAL8Vector *h31i;

    // h30
    REAL8Vector *h30r;
    REAL8Vector *h30i;

    // h44
    REAL8Vector *h44r;
    REAL8Vector *h44i;

    // h43
    REAL8Vector *h43r;
    REAL8Vector *h43i;

    // h42
    REAL8Vector *h42r;
    REAL8Vector *h42i;

    // h41
    REAL8Vector *h41r;
    REAL8Vector *h41i;

    // h40
    REAL8Vector *h40r;
    REAL8Vector *h40i;

    // h55
    REAL8Vector *h55r;
    REAL8Vector *h55i;

    // h54
    REAL8Vector *h54r;
    REAL8Vector *h54i;

    // h53
    REAL8Vector *h53r;
    REAL8Vector *h53i;

    // h52
    REAL8Vector *h52r;
    REAL8Vector *h52i;

    // h51
    REAL8Vector *h51r;
    REAL8Vector *h51i;

    // h50
    REAL8Vector *h50r;
    REAL8Vector *h50i;

}ITF_PNGW_TEvolve;



INT itf_tevolve_orbit_3PN(PNGW3Core *pnparams,
                          ITF_PNGW_TEvolve **outdyn)
{
    // print_debug("EPS_ABS = %.e\n", pnparams->hParams->EPS_ABS);
    // print_debug("EPS_REL = %.e\n", pnparams->hParams->EPS_REL);
    // print_debug("vmax = %.e\n", pnparams->hParams->vmax);
    // print_debug("PNOrder = %d\n", pnparams->hParams->PNOrder);

    // return X_FAILURE;
    PNGW3Dynamics *dyn = NULL;
    if (evolve_binary_orbit_t(pnparams, &dyn) != X_SUCCESS)
        return X_FAILURE;
    if (evaluate_PACorr_from_tdynamics(pnparams->eta, pnparams->initial_chi, dyn) != X_SUCCESS)
        return X_FAILURE;

    // print_debug("len = %d\n", dyn->length);
    (*outdyn) = (ITF_PNGW_TEvolve*)MYMalloc(sizeof(ITF_PNGW_TEvolve));
    (*outdyn)->tVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->tVec->data, dyn->tVec, dyn->length*sizeof(REAL8));
    (*outdyn)->vVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vVec->data, dyn->vVec, dyn->length*sizeof(REAL8));
    (*outdyn)->eVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->eVec->data, dyn->eVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nVec->data, dyn->nVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lVec->data, dyn->lVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdVec->data, dyn->lbdVec, dyn->length*sizeof(REAL8));

    (*outdyn)->chiVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->chiVec->data, dyn->chiVec, dyn->length*sizeof(REAL8));

    (*outdyn)->vPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vPAVec->data, dyn->vPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->ePAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->ePAVec->data, dyn->ePAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nPAVec->data, dyn->nPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lPAVec->data, dyn->lPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdPAVec->data, dyn->lbdPAVec, dyn->length*sizeof(REAL8));

    PNGW3HSphericalModes_TD *hlms = NULL;
    if (evaluate_hSphericalModes_from_PADynamics(pnparams->eta, dyn, &hlms,
        1e-4, 1e-5) != X_SUCCESS)
        return X_FAILURE;
#define CPHLMS(L,M) (*outdyn)->h##L##M##r = CreateREAL8Vector(dyn->length);\
    memcpy((*outdyn)->h##L##M##r->data, hlms->h##L##M##rVec, dyn->length*sizeof(REAL8));\
    (*outdyn)->h##L##M##i = CreateREAL8Vector(dyn->length);\
    memcpy((*outdyn)->h##L##M##i->data, hlms->h##L##M##iVec, dyn->length*sizeof(REAL8));

    CPHLMS(2,2)
    CPHLMS(2,1)
    CPHLMS(2,0)

    CPHLMS(3,3)
    CPHLMS(3,2)
    CPHLMS(3,1)
    CPHLMS(3,0)

    CPHLMS(4,4)
    CPHLMS(4,3)
    CPHLMS(4,2)
    CPHLMS(4,1)
    CPHLMS(4,0)

    CPHLMS(5,5)
    CPHLMS(5,4)
    CPHLMS(5,3)
    CPHLMS(5,2)
    CPHLMS(5,1)
    CPHLMS(5,0)

#undef CPHLMS
    // print_debug("h42 = %.e , %.e, H42 = %e, %e\n", 
    //     creal(hlm), cimag(hlm), hlmr, hlmi);

    STRUCTFREE(dyn, PNGW3Dynamics);
    STRUCTFREE(hlms, PNGW3HSphericalModes_TD);
    return X_SUCCESS;
}

INT itf_tevolve_orbit_pce10_3PN(PNGW3Core *pnparams,
                          ITF_PNGW_TEvolve **outdyn)
{
    // print_debug("EPS_ABS = %.e\n", pnparams->hParams->EPS_ABS);
    // print_debug("EPS_REL = %.e\n", pnparams->hParams->EPS_REL);
    // print_debug("vmax = %.e\n", pnparams->hParams->vmax);
    // print_debug("PNOrder = %d\n", pnparams->hParams->PNOrder);

    // return X_FAILURE;
    PNGW3Dynamics *dyn = NULL;
    if (evolve_binary_orbit_t(pnparams, &dyn) != X_SUCCESS)
        return X_FAILURE;
    if (evaluate_PACorr_from_tdynamics(pnparams->eta, pnparams->initial_chi, dyn) != X_SUCCESS)
        return X_FAILURE;

    // print_debug("len = %d\n", dyn->length);
    (*outdyn) = (ITF_PNGW_TEvolve*)MYMalloc(sizeof(ITF_PNGW_TEvolve));
    (*outdyn)->tVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->tVec->data, dyn->tVec, dyn->length*sizeof(REAL8));
    (*outdyn)->vVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vVec->data, dyn->vVec, dyn->length*sizeof(REAL8));
    (*outdyn)->eVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->eVec->data, dyn->eVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nVec->data, dyn->nVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lVec->data, dyn->lVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdVec->data, dyn->lbdVec, dyn->length*sizeof(REAL8));

    (*outdyn)->chiVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->chiVec->data, dyn->chiVec, dyn->length*sizeof(REAL8));

    (*outdyn)->vPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vPAVec->data, dyn->vPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->ePAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->ePAVec->data, dyn->ePAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nPAVec->data, dyn->nPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lPAVec->data, dyn->lPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdPAVec->data, dyn->lbdPAVec, dyn->length*sizeof(REAL8));

    PNGW3HSphericalModes_TD *hlms = NULL;
    if (evaluate_hSphericalModes_pce10_from_PADynamics(pnparams->eta, dyn, &hlms) != X_SUCCESS)
        return X_FAILURE;
#define CPHLMS(L,M) (*outdyn)->h##L##M##r = CreateREAL8Vector(dyn->length);\
    memcpy((*outdyn)->h##L##M##r->data, hlms->h##L##M##rVec, dyn->length*sizeof(REAL8));\
    (*outdyn)->h##L##M##i = CreateREAL8Vector(dyn->length);\
    memcpy((*outdyn)->h##L##M##i->data, hlms->h##L##M##iVec, dyn->length*sizeof(REAL8));

    CPHLMS(2,2)
    CPHLMS(2,1)
    CPHLMS(2,0)

    CPHLMS(3,3)
    CPHLMS(3,2)
    CPHLMS(3,1)
    CPHLMS(3,0)

    CPHLMS(4,4)
    CPHLMS(4,3)
    CPHLMS(4,2)
    CPHLMS(4,1)
    CPHLMS(4,0)

    CPHLMS(5,5)
    CPHLMS(5,4)
    CPHLMS(5,3)
    CPHLMS(5,2)
    CPHLMS(5,1)
    CPHLMS(5,0)

#undef CPHLMS
    // print_debug("h42 = %.e , %.e, H42 = %e, %e\n", 
    //     creal(hlm), cimag(hlm), hlmr, hlmi);

    STRUCTFREE(dyn, PNGW3Dynamics);
    STRUCTFREE(hlms, PNGW3HSphericalModes_TD);
    return X_SUCCESS;
}

typedef struct {
    INT         length;
    REAL8Vector *tVec;
    REAL8Vector *vVec;
    REAL8Vector *eVec;
    REAL8Vector *nVec;
    REAL8Vector *lVec;
    REAL8Vector *lbdVec;

    REAL8Vector *chiVec;
    REAL8Vector *vPAVec;
    REAL8Vector *ePAVec;
    REAL8Vector *nPAVec;
    REAL8Vector *lPAVec;
    REAL8Vector *lbdPAVec;

    // h22
    REAL8Vector *h22r_full;
    REAL8Vector *h22i_full;

    REAL8Vector *h22r_pce10;
    REAL8Vector *h22i_pce10;

}ITF_PNGW_TEvolve_h22;


INT itf_tevolve_orbit_h22_3PN(PNGW3Core *pnparams,
                          ITF_PNGW_TEvolve_h22 **outdyn)
{
    // print_debug("EPS_ABS = %.e\n", pnparams->hParams->EPS_ABS);
    // print_debug("EPS_REL = %.e\n", pnparams->hParams->EPS_REL);
    // print_debug("vmax = %.e\n", pnparams->hParams->vmax);
    // print_debug("PNOrder = %d\n", pnparams->hParams->PNOrder);

    // return X_FAILURE;
    PNGW3Dynamics *dyn = NULL;
    if (evolve_binary_orbit_t(pnparams, &dyn) != X_SUCCESS)
        return X_FAILURE;
    if (evaluate_PACorr_from_tdynamics(pnparams->eta, pnparams->initial_chi, dyn) != X_SUCCESS)
        return X_FAILURE;

    // print_debug("len = %d\n", dyn->length);
    (*outdyn) = (ITF_PNGW_TEvolve_h22*)MYMalloc(sizeof(ITF_PNGW_TEvolve_h22));
    (*outdyn)->tVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->tVec->data, dyn->tVec, dyn->length*sizeof(REAL8));
    (*outdyn)->vVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vVec->data, dyn->vVec, dyn->length*sizeof(REAL8));
    (*outdyn)->eVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->eVec->data, dyn->eVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nVec->data, dyn->nVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lVec->data, dyn->lVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdVec->data, dyn->lbdVec, dyn->length*sizeof(REAL8));

    (*outdyn)->chiVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->chiVec->data, dyn->chiVec, dyn->length*sizeof(REAL8));

    (*outdyn)->vPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->vPAVec->data, dyn->vPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->ePAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->ePAVec->data, dyn->ePAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->nPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->nPAVec->data, dyn->nPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lPAVec->data, dyn->lPAVec, dyn->length*sizeof(REAL8));
    (*outdyn)->lbdPAVec = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->lbdPAVec->data, dyn->lbdPAVec, dyn->length*sizeof(REAL8));

    PNGW3HSingleSphericalMode_TD *h22 = NULL;
    PNGW3HSingleSphericalMode_TD *h22_pce10 = NULL;

    if (evaluate_h22_from_PADynamics(pnparams->eta, dyn, &h22,
        1e-5, 1e-8) != X_SUCCESS)
        return X_FAILURE;

    if (evaluate_h22_pce10_from_PADynamics(pnparams->eta, dyn, &h22_pce10,
        1e-5, 1e-8) != X_SUCCESS)
        return X_FAILURE;

    // L=2
    (*outdyn)->h22r_full = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->h22r_full->data, h22->hrVec, dyn->length*sizeof(REAL8));
    (*outdyn)->h22i_full = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->h22i_full->data, h22->hiVec, dyn->length*sizeof(REAL8));

    (*outdyn)->h22r_pce10 = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->h22r_pce10->data, h22_pce10->hrVec, dyn->length*sizeof(REAL8));
    (*outdyn)->h22i_pce10 = CreateREAL8Vector(dyn->length);
    memcpy((*outdyn)->h22i_pce10->data, h22_pce10->hiVec, dyn->length*sizeof(REAL8));

    STRUCTFREE(dyn, PNGW3Dynamics);
    STRUCTFREE(h22, PNGW3HSingleSphericalMode_TD);
    STRUCTFREE(h22_pce10, PNGW3HSingleSphericalMode_TD);

    return X_SUCCESS;
}


INT debug_evaluate_hlmpFourierMode(INT p, REAL8 eta, REAL8 v, REAL8 e, REAL8 *out_r, REAL8 *out_i)
{
    PNEllipticEvaluator *pe = CreatePNEllipticEvaluator(e, 1e-5, 1e-8);
    DynParamsCache dc;
    evaluate_DynParamsCache(v, e, &dc);
    COMPLEX16 hlm = calculate_h22FourierModep(p, eta, &dc, pe);
    // calculate_h22FourierModep(p, e, pe, modes);
    // memset(modes, 0, sizeof(hlmpFourierMode));
    DestroyPNEllipticEvaluator(pe);
    *out_r = creal(hlm);
    *out_i = cimag(hlm);
    return X_SUCCESS;
}
