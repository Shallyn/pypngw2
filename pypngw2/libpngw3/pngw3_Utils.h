/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW3_UTILS__
#define __INCLUDE_PNGW3_UTILS__
#include "pngw3_Datatypes.h"
#include "../libpnfourier/pnFourier.h"

#define INTERNAL_DEBUG_LEVEL_PNGW3 get_debug_level_pngw3()
#define SET_INTERNAL_DEBUG_LEVEL_PNGW3(flag) set_debug_level_pngw3(flag)

void set_debug_level_pngw3(INT flag);
INT get_debug_level_pngw3();

#define MAXPOW_EPOW 10
#define MAXPOW_SQ1MESQ 20
#define MAXPOW_BETA 15
#define MAXPOW_BETA2P1 10
#define MAXPOW_BETA2M1 10
#define MAXPOW_V 6

#define CONST_X0 0.824360635350064073424325393907081785826888050
#define CONST_X0P 0.437507968365647931205489897855240998438501594

typedef struct {
    // e
    REAL8 ePow[MAXPOW_EPOW+1];

    // sqrt[1-e^2]
    REAL8 sq1mesqPow[MAXPOW_SQ1MESQ+1];

    // beta = (1-sqrt[1-e^2])/e
    REAL8 betaPow[MAXPOW_BETA+1];

    // beta2p1
    REAL8 beta2p1Pow[MAXPOW_BETA2P1+1];

        // beta2p1
    REAL8 beta2m1Pow[MAXPOW_BETA2M1+1];

    // v
    REAL8 vPow[MAXPOW_V+1];
}DynParamsCache;

PNGW3Core* CreatePNGW3Core(
    REAL8           eta,
    REAL8           initial_e,
    REAL8           initial_v,
    REAL8           initial_chi,
    REAL8           phic,
    REAL8           tc,
    REAL8           iota,
    PNGW3HyperParams *hParams
);
void DestroyPNGW3Core(PNGW3Core *core);

PNGW3Dynamics* CreatePNGW3Dynamics(REAL8Array *arr);
void DestroyPNGW3Dynamics(PNGW3Dynamics* dyn);

PNGW3HSphericalModes_TD* CreatePNGW3HSphericalModes_TD(size_t length);
void DestroyPNGW3HSphericalModes_TD(PNGW3HSphericalModes_TD* hlms);

PNGW3HSingleSphericalMode_TD* CreatePNGW3HSingleSphericalMode_TD(size_t length);
void DestroyPNGW3HSingleSphericalMode_TD(PNGW3HSingleSphericalMode_TD* hlms);

RK4Integrator *CreateRK4Integrator(int dim,
    int (*dydt) (double t, const double y[], double dydt[], void *params),   /* These are XLAL functions! */
    int (*stop) (double t, const double y[], double dydt[], void *params), 
    double eps_abs, double eps_rel);

void DestroyRK4Integrator(RK4Integrator *integrator);

int run_Adaptive_RK4Integrator(REAL8Array **output,
    RK4Integrator *integrator,
    RK4Params * params, REAL8 * yinit, REAL8 tinit, REAL8 tend,
    REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
    BOOLEAN do_interpolate
);

COMPLEX16 WignerDMatrix(
                        int l,        /**< mode number l */
                        int mp,        /**< mode number m' */
                        int m,        /**< mode number m */
                        double alpha,  /**< euler angle (rad) */
                        double beta, /**< euler angle (rad) */
                        double gam  /**< euler angle (rad) */
);
INT SpinWeightedSphericalHarmonic(REAL8 theta,
                                  REAL8 phi,
                                  INT s,
                                  INT l,
                                  INT m,
                                  COMPLEX16 *ret);


#endif
