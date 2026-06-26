/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel:
* 20xx.xx.xx, LOC
**/

#include "pngw3_Utils.h"
#include "pngw3_Datatypes.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

INT INTERNAL_DEBUG_FLAG_PNGW3 = 0;
void set_debug_level_pngw3(INT flag) {
    INTERNAL_DEBUG_FLAG_PNGW3 = flag;
}
INT get_debug_level_pngw3() {
    return INTERNAL_DEBUG_FLAG_PNGW3;
}

PNGW3Core* CreatePNGW3Core(
    REAL8            eta,
    REAL8            initial_e,
    REAL8            initial_v,
    REAL8            initial_chi,
    REAL8            phic,
    REAL8            tc,
    REAL8            iota,
    PNGW3HyperParams *hParams
)
{
    PNGW3Core *ret = (PNGW3Core*)MYMalloc(sizeof(PNGW3Core));
    ret->eta = eta;
    ret->initial_e = initial_e;
    ret->initial_v = initial_v;
    ret->initial_chi = initial_chi;
    ret->phic = phic;
    ret->tc = tc;
    ret->iota = iota;
    ret->hParams = (PNGW3HyperParams*)MYMalloc(sizeof(PNGW3HyperParams));
    memcpy(ret->hParams, hParams, sizeof(PNGW3HyperParams));
    return ret;
}

void DestroyPNGW3Core(PNGW3Core *core)
{
    if (!core) return;
    if (core->hParams) MYFree(core->hParams);
    MYFree(core);
    core = NULL;
    return;
}

PNGW3Dynamics* CreatePNGW3Dynamics(REAL8Array *arr)
{
    if (!arr) return NULL;
    if (arr->dimLength->length != 2) return NULL;
    if (arr->dimLength->data[0] != PNGW3_NDYNAMICS+1) return NULL;
    PNGW3Dynamics* dyn = (PNGW3Dynamics*) MYMalloc(sizeof(PNGW3Dynamics));
    if (!dyn) {
        X_ERROR_NULL(X_ENOMEM);
    }
    dyn->array = arr;
    dyn->length = arr->dimLength->data[1];
    dyn->tVec = arr->data;
    dyn->vVec = dyn->tVec + dyn->length;
    dyn->eVec = dyn->vVec + dyn->length;
    dyn->nVec = dyn->eVec + dyn->length;
    dyn->lVec = dyn->nVec + dyn->length;
    dyn->lbdVec = dyn->lVec + dyn->length;

    dyn->array_extended = NULL;
    dyn->chiVec = NULL;
    dyn->vPAVec = NULL;
    dyn->ePAVec = NULL;
    dyn->nPAVec = NULL;
    dyn->lPAVec = NULL;
    dyn->lbdPAVec = NULL;
    return dyn;
}

void DestroyPNGW3Dynamics(PNGW3Dynamics* dyn)
{
    if (!dyn) return;
    STRUCTFREE(dyn->array, REAL8Array);
    STRUCTFREE(dyn->array_extended, REAL8Array);
    MYFree(dyn);    
    return;
}

PNGW3HSphericalModes_TD* CreatePNGW3HSphericalModes_TD(size_t length)
{
    PNGW3HSphericalModes_TD *ret = (PNGW3HSphericalModes_TD*)MYMalloc(sizeof(PNGW3HSphericalModes_TD));
    ret->array = CreateREAL8ArrayL(2, PNGW3_NSPHERICALMODES+1, length);
    ret->tVec = ret->array->data;
    // L=2
    ret->h22rVec = ret->tVec + length;
    ret->h22iVec = ret->h22rVec + length;

    ret->h21rVec = ret->h22iVec + length;
    ret->h21iVec = ret->h21rVec + length;

    ret->h20rVec = ret->h21iVec + length;
    ret->h20iVec = ret->h20rVec + length;

    // L=3
    ret->h33rVec = ret->h20iVec + length;
    ret->h33iVec = ret->h33rVec + length;

    ret->h32rVec = ret->h33iVec + length;
    ret->h32iVec = ret->h32rVec + length;

    ret->h31rVec = ret->h32iVec + length;
    ret->h31iVec = ret->h31rVec + length;

    ret->h30rVec = ret->h31iVec + length;
    ret->h30iVec = ret->h30rVec + length;

    // L=4
    ret->h44rVec = ret->h30iVec + length;
    ret->h44iVec = ret->h44rVec + length;

    ret->h43rVec = ret->h44iVec + length;
    ret->h43iVec = ret->h43rVec + length;

    ret->h42rVec = ret->h43iVec + length;
    ret->h42iVec = ret->h42rVec + length;

    ret->h41rVec = ret->h42iVec + length;
    ret->h41iVec = ret->h41rVec + length;

    ret->h40rVec = ret->h41iVec + length;
    ret->h40iVec = ret->h40rVec + length;

    // L=5
    ret->h55rVec = ret->h40iVec + length;
    ret->h55iVec = ret->h55rVec + length;

    ret->h54rVec = ret->h55iVec + length;
    ret->h54iVec = ret->h54rVec + length;

    ret->h53rVec = ret->h54iVec + length;
    ret->h53iVec = ret->h53rVec + length;

    ret->h52rVec = ret->h53iVec + length;
    ret->h52iVec = ret->h52rVec + length;

    ret->h51rVec = ret->h52iVec + length;
    ret->h51iVec = ret->h51rVec + length;

    ret->h50rVec = ret->h51iVec + length;
    ret->h50iVec = ret->h50rVec + length;

    return ret;
}

void DestroyPNGW3HSphericalModes_TD(PNGW3HSphericalModes_TD* hlms)
{
    if (!hlms) return;
    STRUCTFREE(hlms->array, REAL8Array);
    MYFree(hlms);    
    return;
}


PNGW3HSingleSphericalMode_TD* CreatePNGW3HSingleSphericalMode_TD(size_t length)
{
    PNGW3HSingleSphericalMode_TD *ret = (PNGW3HSingleSphericalMode_TD*)MYMalloc(sizeof(PNGW3HSingleSphericalMode_TD));
    ret->array = CreateREAL8ArrayL(2, 3, length);
    ret->tVec = ret->array->data;
    // L=2
    ret->hrVec = ret->tVec + length;
    ret->hiVec = ret->hrVec + length;

    return ret;
}

void DestroyPNGW3HSingleSphericalMode_TD(PNGW3HSingleSphericalMode_TD* hlms)
{
    if (!hlms) return;
    STRUCTFREE(hlms->array, REAL8Array);
    MYFree(hlms);    
    return;
}

#define BEGIN_GSL_ERRHANDLER \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define END_GSL_ERRHANDLER \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
        }

RK4Integrator *CreateRK4Integrator(
    int dim,
    int (*dydt) (double t, const double y[], double dydt[], void *params),   /* These are XLAL functions! */
    int (*stop) (double t, const double y[], double dydt[], void *params), 
    double eps_abs, double eps_rel
)
{
    RK4Integrator *integrator;

    /* allocate our custom integrator structure */
    if (!(integrator = (RK4Integrator *) MYCalloc(1, sizeof(RK4Integrator)))) 
    {
        return NULL;
    }

    /* allocate the GSL ODE components */
    integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, dim);
    //XLAL_CALLGSL(integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, dim));
    integrator->control = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    integrator->evolve = gsl_odeiv_evolve_alloc(dim);

    /* allocate the GSL system (functions, etc.) */
    integrator->sys = (gsl_odeiv_system *) MYCalloc(1, sizeof(gsl_odeiv_system));

    /* if something failed to be allocated, bail out */
    if (!(integrator->step) || !(integrator->control) || !(integrator->evolve) || !(integrator->sys)) 
    {
        STRUCTFREE(integrator, RK4Integrator);
        return NULL;
    }

    integrator->dydt = dydt;
    integrator->stop = stop;

    integrator->sys->function = dydt;
    integrator->sys->jacobian = NULL;
    integrator->sys->dimension = dim;
    integrator->sys->params = NULL;

    integrator->retries = 6;
    integrator->stopontestonly = 0;
    return integrator;
}

void DestroyRK4Integrator(RK4Integrator *integrator)
{
    if (!integrator)
        return;
    if (integrator->evolve)
        gsl_odeiv_evolve_free(integrator->evolve);
    if (integrator->control)
        gsl_odeiv_control_free(integrator->control);
    if (integrator->step)
        gsl_odeiv_step_free(integrator->step);
    MYFree(integrator->sys);
    MYFree(integrator);
    return;
}


int run_Adaptive_RK4Integrator(REAL8Array **output,
    RK4Integrator *integrator,
    RK4Params * params, REAL8 * yinit, REAL8 tinit, REAL8 tend,
    REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
    BOOLEAN do_interpolate
)
{
    int errnum = X_SUCCESS;
    int status; /* used throughout */
    register unsigned int i, j;

    /* needed for the integration */
    size_t dim, outputlength=0, bufferlength, retries;
    REAL8 t, tnew, h0, h0old;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *out_array = NULL;
    REAL8 *times, *vector;      /* aliases */

    BEGIN_GSL_ERRHANDLER;

    /* allocate the buffers!
        * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
        * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat_or_h0) + 2;   /* allow for the initial value and possibly a final semi-step */
    UINT dimn;/* Variable for different loop indices below */
    UINT ncache = params->cache ? params->cache->length : 0;
    dimn = dim + 1 + ncache;

    buffers = CreateREAL8ArrayL(2, dimn, bufferlength); 
    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp ) 
    {
        errnum = X_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = (void*) params;
    
    integrator->returncode = 0;

    t = tinit;
    h0 = deltat_or_h0;

    h0old = h0; /* initialized so that it will not trigger the check h0<h0old at the first step */
    memcpy(y, yinit, dim * sizeof(REAL8));
    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = X_FAILURE;
        goto bail_out;
    }

    if (ncache > 0)
        for (i=1; i<=dim; i++)
            buffers->data[(i+dim)*bufferlength] = params->cache->data[i];

    X_PRINT_INFO("Run RK4 integrator ...");
    while (1) 
    {
        if (integrator->stopontestonly && t >= tend)
            break;
        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) 
            {
                integrator->returncode = status;
                break;
            }
        }
        /* ready to try stepping! */

try_step:
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));
        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
            * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
            * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* Enforce a minimal allowed time step */
        /* To ignore this type of constraint, set min_deltat_or_h0 = 0 */
        if (h0 < min_deltat_or_h0) h0 = min_deltat_or_h0;
        /* did the error-checker reduce the stepsize?
            * note: previously, was using status == GSL_ODEIV_HADJ_DEC
            * other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
            * GSL_ODEIV_HADJ_NIL if it was unchanged
            * since we introduced a minimal step size we simply compare to the saved value of h0 */
        if (h0 < h0old) 
        {
            h0old = h0;
            // print_debug("[%zu/%zu]h_old = %.16e, h_new = %.16e\n", outputlength, bufferlength, h0old, h0);
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives, save the time step */
        t = tnew;
        h0old = h0;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        outputlength++;

        if (outputlength >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8ArrayL(2, dimn, 2 * bufferlength))) 
            {
                errnum = X_ENOMEM;   /* ouch, that hurt */
                goto bail_out;
            } else {
                for ( i = 0; i < dimn /* dim */; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], outputlength * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into output buffers */
        buffers->data[outputlength] = t;
        for ( i = 1; i <= dim; i++)
            buffers->data[i * bufferlength + outputlength] = y[i - 1];   /* y does not have time */
        
        if (ncache > 0)
            for (i=1; i<=dim; i++)
                buffers->data[(i+dim)*bufferlength + outputlength] = params->cache->data[i];

    } // End Loop
    memcpy(yinit, y, dim * sizeof(REAL8));
    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (outputlength == 0)
        goto bail_out;
    if (do_interpolate) {
        interp = gsl_spline_alloc(gsl_interp_cspline, outputlength + 1);
        accel = gsl_interp_accel_alloc();
        outputlen = (int)(t / deltat_or_h0) + 1;
        out_array = CreateREAL8ArrayL(2, dimn, outputlen);
        if (!interp || !accel || !out_array) {
            errnum = X_ENOMEM;   /* ouch again, ran out of memory */
            STRUCTFREE(out_array, REAL8Array);
            outputlen = 0;
            goto bail_out;
        }
        /* make an array of times */
        times = out_array->data;
        for ( j = 0; j < outputlen; j++)
            times[j] = deltat_or_h0 * j;
        /* interpolate! */
        for ( i = 1; i <= dimn-1; i++) {
            gsl_spline_init(interp, &buffers->data[0], &buffers->data[bufferlength * i], outputlength + 1);

            vector = out_array->data + outputlen * i;
            for ( j = 0; j < outputlen; j++) {
                gsl_spline_eval_e(interp, times[j], accel, &(vector[j]));
            }
        }
    } else {
        out_array = CreateREAL8ArrayL(2, dimn, outputlength);
        if (!out_array)
        {
            errnum = X_ENOMEM;   /* ouch again, ran out of memory */
            STRUCTFREE(out_array, REAL8Array);
            outputlength = 0;
            goto bail_out;
        }
        for( j=0;j<outputlength;j++) 
        {
            // (*t_and_y_and_dydt_out)->data[j] = buffers->data[j];
            for( i=0;i<dimn;i++) 
                out_array->data[i*outputlength + j] = buffers->data[i*bufferlength + j];

        }

    }
bail_out:
    END_GSL_ERRHANDLER;
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);
    if (interp)
        gsl_spline_free(interp);
    if (accel)
        gsl_interp_accel_free(accel);
    if (errnum)
        return errnum;
    *output = out_array;
    return X_SUCCESS;
}


static COMPLEX16 CX16polar(REAL8 r, REAL8 phi)
{
    REAL8 re,im;
    re = r*cos(phi);
    im = r*sin(phi);
    
    
    COMPLEX16 z=re+I*im;
    
    return z;
}

INT SpinWeightedSphericalHarmonic(REAL8 theta,
                                  REAL8 phi,
                                  INT s,
                                  INT l,
                                  INT m,
                                  COMPLEX16 *ret)
{
    REAL8 fac;
    COMPLEX16 ans;
    
    /* sanity checks ... */
    if ( l < abs(s) )
    {
        XPrintError("Error -%s : Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", __func__, s, l, m);
        return X_FAILURE;
    }
    if ( l < abs(m) )
    {
        XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n" ,__func__, s, l, m );
        return X_FAILURE;
    }
    
    if ( s == -2 )
    {
        if ( l == 2 )
        {
            switch ( m )
            {
                case -2:
                    fac = sqrt( 5.0 / ( 64.0 * CST_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
                    break;
                case -1:
                    fac = sqrt( 5.0 / ( 16.0 * CST_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
                    break;
                    
                case 0:
                    fac = sqrt( 15.0 / ( 32.0 * CST_PI ) ) * sin( theta )*sin( theta );
                    break;
                    
                case 1:
                    fac = sqrt( 5.0 / ( 16.0 * CST_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
                    break;
                    
                case 2:
                    fac = sqrt( 5.0 / ( 64.0 * CST_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
                    break;
                default:
                    XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            } /*  switch (m) */
        }  /* l==2*/
        else if ( l == 3 )
        {
            switch ( m )
            {
                case -3:
                    fac = sqrt(21.0/(2.0*CST_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
                    break;
                case -2:
                    fac = sqrt(7.0/(4.0*CST_PI))*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
                    break;
                case -1:
                    fac = sqrt(35.0/(2.0*CST_PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
                    break;
                case 0:
                    fac = (sqrt(105.0/(2.0*CST_PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
                    break;
                case 1:
                    fac = -sqrt(35.0/(2.0*CST_PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
                    break;
                    
                case 2:
                    fac = sqrt(7.0/CST_PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
                    break;
                    
                case 3:
                    fac = -sqrt(21.0/(2.0*CST_PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
                    break;
                    
                default:
                    XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        }   /* l==3 */
        else if ( l == 4 )
        {
            switch ( m )
            {
                case -4:
                    fac = 3.0*sqrt(7.0/CST_PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
                    break;
                case -3:
                    fac = 3.0*sqrt(7.0/(2.0*CST_PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
                    break;
                    
                case -2:
                    fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(CST_PI));
                    break;
                case -1:
                    fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*CST_PI));
                    break;
                case 0:
                    fac = (3.0*sqrt(5.0/(2.0*CST_PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
                    break;
                case 1:
                    fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*CST_PI));
                    break;
                case 2:
                    fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(CST_PI));
                    break;
                case 3:
                    fac = -3.0*sqrt(7.0/(2.0*CST_PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
                    break;
                case 4:
                    fac = 3.0*sqrt(7.0/CST_PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
                    break;
                default:
                    XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        }    /* l==4 */
        else if ( l == 5 )
        {
            switch ( m )
            {
                case -5:
                    fac = sqrt(330.0/CST_PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
                    break;
                case -4:
                    fac = sqrt(33.0/CST_PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
                    break;
                case -3:
                    fac = (sqrt(33.0/(2.0*CST_PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
                    break;
                case -2:
                    fac = (sqrt(11.0/CST_PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
                    break;
                case -1:
                    fac = (sqrt(77.0/CST_PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
                    break;
                case 0:
                    fac = (sqrt(1155.0/(2.0*CST_PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
                    break;
                case 1:
                    fac = sqrt(77.0/CST_PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
                    break;
                case 2:
                    fac = sqrt(11.0/CST_PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
                    break;
                case 3:
                    fac = -sqrt(33.0/(2.0*CST_PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
                    break;
                case 4:
                    fac = sqrt(33.0/CST_PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
                    break;
                case 5:
                    fac = -sqrt(330.0/CST_PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
                    break;
                default:
                    XPrintError("Error -%s: nvalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        }  /* l==5 */
        else if ( l == 6 )
        {
            switch ( m )
            {
                case -6:
                    fac = (3.*sqrt(715./CST_PI)*pow(cos(theta/2.0),4)*pow(sin(theta/2.0),8))/2.0;
                    break;
                case -5:
                    fac = (sqrt(2145./CST_PI)*pow(cos(theta/2.0),3)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),7))/2.0;
                    break;
                case -4:
                    fac = (sqrt(195./(2.0*CST_PI))*pow(cos(theta/2.0),2)*(35. + 44.*cos(theta)
                                                                                      + 33.*cos(2.*theta))*pow(sin(theta/2.0),6))/8.0;
                    break;
                case -3:
                    fac = (3.*sqrt(13./CST_PI)*cos(theta/2.0)*(98. + 185.*cos(theta) + 110.*cos(2*theta)
                                                                       + 55.*cos(3.*theta))*pow(sin(theta/2.0),5))/32.0;
                    break;
                case -2:
                    fac = (sqrt(13./CST_PI)*(1709. + 3096.*cos(theta) + 2340.*cos(2.*theta) + 1320.*cos(3.*theta)
                                                 + 495.*cos(4.*theta))*pow(sin(theta/2.0),4))/256.0;
                    break;
                case -1:
                    fac = (sqrt(65./(2.0*CST_PI))*cos(theta/2.0)*(161. + 252.*cos(theta) + 252.*cos(2.*theta)
                                                                          + 132.*cos(3.*theta) + 99.*cos(4.*theta))*pow(sin(theta/2.0),3))/64.0;
                    break;
                case 0:
                    fac = (sqrt(1365./CST_PI)*(35. + 60.*cos(2.*theta) + 33.*cos(4.*theta))*pow(sin(theta),2))/512.0;
                    break;
                case 1:
                    fac = (sqrt(65./(2.0*CST_PI))*pow(cos(theta/2.0),3)*(161. - 252.*cos(theta) + 252.*cos(2.*theta)
                                                                                     - 132.*cos(3.*theta) + 99.*cos(4.*theta))*sin(theta/2.0))/64.0;
                    break;
                case 2:
                    fac = (sqrt(13./CST_PI)*pow(cos(theta/2.0),4)*(1709. - 3096.*cos(theta) + 2340.*cos(2.*theta)
                                                                               - 1320*cos(3*theta) + 495*cos(4*theta)))/256.0;
                    break;
                case 3:
                    fac = (-3.*sqrt(13./CST_PI)*pow(cos(theta/2.0),5)*(-98. + 185.*cos(theta) - 110.*cos(2*theta)
                                                                                   + 55.*cos(3.*theta))*sin(theta/2.0))/32.0;
                    break;
                case 4:
                    fac = (sqrt(195./(2.0*CST_PI))*pow(cos(theta/2.0),6)*(35. - 44.*cos(theta)
                                                                                      + 33.*cos(2*theta))*pow(sin(theta/2.0),2))/8.0;
                    break;
                case 5:
                    fac = -(sqrt(2145./CST_PI)*pow(cos(theta/2.0),7)*(-1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
                    break;
                case 6:
                    fac = (3.*sqrt(715./CST_PI)*pow(cos(theta/2.0),8)*pow(sin(theta/2.0),4))/2.0;
                    break;
                default:
                    XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        } /* l==6 */
        else if ( l == 7 )
        {
            switch ( m )
            {
                case -7:
                    fac = sqrt(15015./(2.0*CST_PI))*pow(cos(theta/2.0),5)*pow(sin(theta/2.0),9);
                    break;
                case -6:
                    fac = (sqrt(2145./CST_PI)*pow(cos(theta/2.0),4)*(2. + 7.*cos(theta))*pow(sin(theta/2.0),8))/2.0;
                    break;
                case -5:
                    fac = (sqrt(165./(2.0*CST_PI))*pow(cos(theta/2.0),3)*(93. + 104.*cos(theta)
                                                                                      + 91.*cos(2.*theta))*pow(sin(theta/2.0),7))/8.0;
                    break;
                case -4:
                    fac = (sqrt(165./(2.0*CST_PI))*pow(cos(theta/2.0),2)*(140. + 285.*cos(theta)
                                                                                      + 156.*cos(2.*theta) + 91.*cos(3.*theta))*pow(sin(theta/2.0),6))/16.0;
                    break;
                case -3:
                    fac = (sqrt(15./(2.0*CST_PI))*cos(theta/2.0)*(3115. + 5456.*cos(theta) + 4268.*cos(2.*theta)
                                                                          + 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*pow(sin(theta/2.0),5))/128.0;
                    break;
                case -2:
                    fac = (sqrt(15./CST_PI)*(5220. + 9810.*cos(theta) + 7920.*cos(2.*theta) + 5445.*cos(3.*theta)
                                                 + 2860.*cos(4.*theta) + 1001.*cos(5.*theta))*pow(sin(theta/2.0),4))/512.0;
                    break;
                case -1:
                    fac = (3.*sqrt(5./(2.0*CST_PI))*cos(theta/2.0)*(1890. + 4130.*cos(theta) + 3080.*cos(2.*theta)
                                                                            + 2805.*cos(3.*theta) + 1430.*cos(4.*theta) + 1001.*cos(5*theta))*pow(sin(theta/2.0),3))/512.0;
                    break;
                case 0:
                    fac = (3.*sqrt(35./CST_PI)*cos(theta)*(109. + 132.*cos(2.*theta)
                                                                   + 143.*cos(4.*theta))*pow(sin(theta),2))/512.0;
                    break;
                case 1:
                    fac = (3.*sqrt(5./(2.0*CST_PI))*pow(cos(theta/2.0),3)*(-1890. + 4130.*cos(theta) - 3080.*cos(2.*theta)
                                                                                       + 2805.*cos(3.*theta) - 1430.*cos(4.*theta) + 1001.*cos(5.*theta))*sin(theta/2.0))/512.0;
                    break;
                case 2:
                    fac = (sqrt(15./CST_PI)*pow(cos(theta/2.0),4)*(-5220. + 9810.*cos(theta) - 7920.*cos(2.*theta)
                                                                               + 5445.*cos(3.*theta) - 2860.*cos(4.*theta) + 1001.*cos(5.*theta)))/512.0;
                    break;
                case 3:
                    fac = -(sqrt(15./(2.0*CST_PI))*pow(cos(theta/2.0),5)*(3115. - 5456.*cos(theta) + 4268.*cos(2.*theta)
                                                                                      - 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*sin(theta/2.0))/128.0;
                    break;
                case 4:
                    fac = (sqrt(165./(2.0*CST_PI))*pow(cos(theta/2.0),6)*(-140. + 285.*cos(theta) - 156.*cos(2*theta)
                                                                                      + 91.*cos(3.*theta))*pow(sin(theta/2.0),2))/16.0;
                    break;
                case 5:
                    fac = -(sqrt(165./(2.0*CST_PI))*pow(cos(theta/2.0),7)*(93. - 104.*cos(theta)
                                                                                       + 91.*cos(2.*theta))*pow(sin(theta/2.0),3))/8.0;
                    break;
                case 6:
                    fac = (sqrt(2145./CST_PI)*pow(cos(theta/2.0),8)*(-2. + 7.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
                    break;
                case 7:
                    fac = -(sqrt(15015./(2.0*CST_PI))*pow(cos(theta/2.0),9)*pow(sin(theta/2.0),5));
                    break;
                default:
                    XPrintError("Error %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        } /* l==7 */
        else if ( l == 8 )
        {
            switch ( m )
            {
                case -8:
                    fac = sqrt(34034./CST_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),10);
                    break;
                case -7:
                    fac = sqrt(17017./(2.0*CST_PI))*pow(cos(theta/2.0),5)*(1. + 4.*cos(theta))*pow(sin(theta/2.0),9);
                    break;
                case -6:
                    fac = sqrt(255255./CST_PI)*pow(cos(theta/2.0),4)*(1. + 2.*cos(theta))
                    *sin(CST_PI/4.0 - theta/2.0)*sin(CST_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),8);
                    break;
                case -5:
                    fac = (sqrt(12155./(2.0*CST_PI))*pow(cos(theta/2.0),3)*(19. + 42.*cos(theta)
                                                                                        + 21.*cos(2.*theta) + 14.*cos(3.*theta))*pow(sin(theta/2.0),7))/8.0;
                    break;
                case -4:
                    fac = (sqrt(935./(2.0*CST_PI))*pow(cos(theta/2.0),2)*(265. + 442.*cos(theta) + 364.*cos(2.*theta)
                                                                                      + 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),6))/32.0;
                    break;
                case -3:
                    fac = (sqrt(561./(2.0*CST_PI))*cos(theta/2.0)*(869. + 1660.*cos(theta) + 1300.*cos(2.*theta)
                                                                           + 910.*cos(3.*theta) + 455.*cos(4.*theta) + 182.*cos(5.*theta))*pow(sin(theta/2.0),5))/128.0;
                    break;
                case -2:
                    fac = (sqrt(17./CST_PI)*(7626. + 14454.*cos(theta) + 12375.*cos(2.*theta) + 9295.*cos(3.*theta)
                                                 + 6006.*cos(4.*theta) + 3003.*cos(5.*theta) + 1001.*cos(6.*theta))*pow(sin(theta/2.0),4))/512.0;
                    break;
                case -1:
                    fac = (sqrt(595./(2.0*CST_PI))*cos(theta/2.0)*(798. + 1386.*cos(theta) + 1386.*cos(2.*theta)
                                                                           + 1001.*cos(3.*theta) + 858.*cos(4.*theta) + 429.*cos(5.*theta) + 286.*cos(6.*theta))*pow(sin(theta/2.0),3))/512.0;
                    break;
                case 0:
                    fac = (3.*sqrt(595./CST_PI)*(210. + 385.*cos(2.*theta) + 286.*cos(4.*theta)
                                                     + 143.*cos(6.*theta))*pow(sin(theta),2))/4096.0;
                    break;
                case 1:
                    fac = (sqrt(595./(2.0*CST_PI))*pow(cos(theta/2.0),3)*(798. - 1386.*cos(theta) + 1386.*cos(2.*theta)
                                                                                      - 1001.*cos(3.*theta) + 858.*cos(4.*theta) - 429.*cos(5.*theta) + 286.*cos(6.*theta))*sin(theta/2.0))/512.0;
                    break;
                case 2:
                    fac = (sqrt(17./CST_PI)*pow(cos(theta/2.0),4)*(7626. - 14454.*cos(theta) + 12375.*cos(2.*theta)
                                                                               - 9295.*cos(3.*theta) + 6006.*cos(4.*theta) - 3003.*cos(5.*theta) + 1001.*cos(6.*theta)))/512.0;
                    break;
                case 3:
                    fac = -(sqrt(561./(2.0*CST_PI))*pow(cos(theta/2.0),5)*(-869. + 1660.*cos(theta) - 1300.*cos(2.*theta)
                                                                                       + 910.*cos(3.*theta) - 455.*cos(4.*theta) + 182.*cos(5.*theta))*sin(theta/2.0))/128.0;
                    break;
                case 4:
                    fac = (sqrt(935./(2.0*CST_PI))*pow(cos(theta/2.0),6)*(265. - 442.*cos(theta) + 364.*cos(2.*theta)
                                                                                      - 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),2))/32.0;
                    break;
                case 5:
                    fac = -(sqrt(12155./(2.0*CST_PI))*pow(cos(theta/2.0),7)*(-19. + 42.*cos(theta) - 21.*cos(2.*theta)
                                                                                         + 14.*cos(3.*theta))*pow(sin(theta/2.0),3))/8.0;
                    break;
                case 6:
                    fac = sqrt(255255./CST_PI)*pow(cos(theta/2.0),8)*(-1. + 2.*cos(theta))*sin(CST_PI/4.0 - theta/2.0)
                    *sin(CST_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),4);
                    break;
                case 7:
                    fac = -(sqrt(17017./(2.0*CST_PI))*pow(cos(theta/2.0),9)*(-1. + 4.*cos(theta))*pow(sin(theta/2.0),5));
                    break;
                case 8:
                    fac = sqrt(34034./CST_PI)*pow(cos(theta/2.0),10)*pow(sin(theta/2.0),6);
                    break;
                default:
                    XPrintError("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return X_FAILURE;
                    break;
            }
        } /* l==8 */
        else
        {
            XPrintError("Error -%s: Unsupported mode l=%d (only l in [2,8] implemented)\n", __func__, l);
            return X_FAILURE;
        }
    }
    else
    {
        XPrintError("Error -%s: Unsupported mode s=%d (only s=-2 implemented)\n", __func__, s);
        return X_FAILURE;
    }
    if (m)
        ans = CX16polar(1.0, m*phi) * fac;
    else
        ans = fac;
    *ret = ans;
    return X_SUCCESS;
}

/**
 * Computes the n-th Jacobi polynomial for polynomial weights alpha and beta.
 * The implementation here is only valid for real x -- enforced by the argument
 * type. An extension to complex values would require evaluation of several
 * gamma functions.
 *
 * See http://en.wikipedia.org/wiki/Jacobi_polynomials
 */
double XLALJacobiPolynomial(int n, int alpha, int beta, double x)
{
    double f1 = (x-1)/2.0, f2 = (x+1)/2.0;
    int s=0;
    double sum=0, val=0;
    if( n == 0 ) return 1.0;
    for( s=0; n-s >= 0; s++ )
    {
        val=1.0;
        val *= gsl_sf_choose( n+alpha, s );
        val *= gsl_sf_choose( n+beta, n-s );
        if( n-s != 0 ) val *= pow( f1, n-s );
        if( s != 0 ) val*= pow( f2, s );

        sum += val;
    }
    return sum;
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/**
 * These functions compute the amplitude and phase of a Wigner coefficient
 * Dlmmp, given Euler angles of an active rotation.
 */
double WignerdMatrix(
                    int l,        /**< mode number l */
                    int mp,        /**< mode number m' */
                    int m,        /**< mode number m */
                    double beta  /**< euler angle (rad) */
)
{

	int k = MIN( l+m, MIN( l-m, MIN( l+mp, l-mp )));
	double a=0, lam=0;
	if(k == l+m){
		a = mp-m;
		lam = mp-m;
	} else if(k == l-m) {
		a = m-mp;
		lam = 0;
	} else if(k == l+mp) {
		a = m-mp;
		lam = 0;
	} else if(k == l-mp) {
		a = mp-m;
		lam = mp-m;
	}

	int b = 2*l-2*k-a;
	double pref = pow(-1, lam) * sqrt(gsl_sf_choose( 2*l-k, k+a )) / sqrt(gsl_sf_choose( k+b, b ));

	return pref * pow(sin(beta/2.0), a) * pow( cos(beta/2.0), b) * XLALJacobiPolynomial(k, a, b, cos(beta));
}

/**
 * Computes the full Wigner D matrix for the Euler angle alpha, beta, and gamma
 * with major index 'l' and minor index transition from m to mp.
 *
 * Uses a slightly unconventional method since the intuitive version by Wigner
 * is less suitable to algorthmic development.
 *
 * See http://en.wikipedia.org/wiki/Wigner_D-matrix
 *
 * Currently only supports the modes which are implemented for the spin
 * weighted spherical harmonics.
 */
COMPLEX16 WignerDMatrix(
                        int l,        /**< mode number l */
                        int mp,        /**< mode number m' */
                        int m,        /**< mode number m */
                        double alpha,  /**< euler angle (rad) */
                        double beta, /**< euler angle (rad) */
                        double gam  /**< euler angle (rad) */
)
{
	 return cexp( -I*mp*alpha ) *
			WignerdMatrix( l, mp, m, beta ) * 
			cexp( -I*m*gam );
}
#undef MIN