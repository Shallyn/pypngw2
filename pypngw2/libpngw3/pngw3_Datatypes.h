/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PNGW3_DATATYPES__
#define __INCLUDE_PNGW3_DATATYPES__
#include "../libcore/pngw_Datatypes.h"
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#define VISCO 0.40824829046386301636621401245098

typedef struct {
    REAL8   EPS_ABS;
    REAL8   EPS_REL;
    UINT    PNOrder;
    REAL8   vmax;
}PNGW3HyperParams;

typedef struct {
    // Core Parameters
    REAL8 eta;              // mass ratio
    REAL8 initial_e;        // initial orbital eccentricity
    REAL8 initial_v;        // initial orbital (averaged) frequency
    REAL8 initial_chi;      // initial orbital relativistic anomaly
    REAL8 phic;             // coalescing phase
    REAL8 tc;               // coalescing time
    REAL8 iota;             // inclination angle

    // Aux Parameters
    PNGW3HyperParams *hParams;
}PNGW3Core;


#define PNGW3_NDYNAMICS 5
typedef struct {
    size_t      length;
    REAL8Array  *array;

    // evolve
    REAL8 *tVec;

    // dynamics
    REAL8 *vVec;
    REAL8 *eVec;
    REAL8 *nVec;
    REAL8 *lVec;
    REAL8 *lbdVec;

    // extended dynamics
    REAL8Array *array_extended;
    REAL8 *chiVec;
    REAL8 *vPAVec;
    REAL8 *ePAVec;
    REAL8 *nPAVec;
    REAL8 *lPAVec;
    REAL8 *lbdPAVec;
}PNGW3Dynamics;

typedef struct {
    size_t length;
    REAL8Array *array;
    REAL8* tVec;
    REAL8* hrVec;
    REAL8* hiVec;
}PNGW3HSingleSphericalMode_TD;

#define PNGW3_NSPHERICALMODES 36 // 2*(3+4+5+6)
typedef struct {
    size_t              length;
    REAL8Array          *array; 

    REAL8       *tVec;
    // L=2
    REAL8       *h22rVec;
    REAL8       *h22iVec;

    REAL8       *h21rVec;
    REAL8       *h21iVec;

    REAL8       *h20rVec;
    REAL8       *h20iVec;

    // L=3
    REAL8       *h33rVec;
    REAL8       *h33iVec;

    REAL8       *h32rVec;
    REAL8       *h32iVec;

    REAL8       *h31rVec;
    REAL8       *h31iVec;

    REAL8       *h30rVec;
    REAL8       *h30iVec;

    // L=4
    REAL8       *h44rVec;
    REAL8       *h44iVec;

    REAL8       *h43rVec;
    REAL8       *h43iVec;

    REAL8       *h42rVec;
    REAL8       *h42iVec;

    REAL8       *h41rVec;
    REAL8       *h41iVec;

    REAL8       *h40rVec;
    REAL8       *h40iVec;

    // L=5
    REAL8       *h55rVec;
    REAL8       *h55iVec;

    REAL8       *h54rVec;
    REAL8       *h54iVec;

    REAL8       *h53rVec;
    REAL8       *h53iVec;

    REAL8       *h52rVec;
    REAL8       *h52iVec;

    REAL8       *h51rVec;
    REAL8       *h51iVec;

    REAL8       *h50rVec;
    REAL8       *h50iVec;
}PNGW3HSphericalModes_TD;

typedef struct {
    REAL8Vector *cache;
    void *params;
}RK4Params;

typedef struct
{
    gsl_odeiv_step    *step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve  *evolve;

    gsl_odeiv_system  *sys;

    int (* dydt) (double t, const double y[], double dydt[], void * params);
    int (* stop) (double t, const double y[], double dydt[], void * params);

    int retries;		/* retries with smaller step when derivatives encounter singularity */
    int stopontestonly;	/* stop only on test, use tend to size buffers only */
    
    int returncode;
} RK4Integrator;



#endif
