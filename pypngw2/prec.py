#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed, 07 Aug 2024 02:13:40 +0000

@author: Shallyn
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from optparse import OptionParser
from pathlib import Path
from scipy.integrate import solve_ivp, RK45

#pwd = Path(__file__).absolute().parent
#pwd = Path(sys.path[0])
pwd = Path(os.getcwd())

class pyPrecDyn(object):
    def __init__(self, tVec:np.ndarray, DynArray:np.ndarray):
        raise NotImplementedError('pyPrecDyn is not Implemented')

class pyPrecDynVec(object):
    def __init__(self, eta:float, y:np.ndarray):
        if len(y) != 12:
            raise ValueError(f'length of y should be 12, rather than {len(y)}')
        self.rVec = y[:3]
        self.vVec = y[3:6]
        self.chi1Vec = y[6:9]
        self.chi2Vec = y[9:]
        self.r2 = r2 = self.rVec.dot(self.rVec)
        self.r = np.sqrt(r2)
        v2 = self.vVec.dot(self.vVec)
        self.v = np.sqrt(v2)
        self.vr = self.rVec.dot(self.vVec)/self.r
        self.LNVec = np.cross(self.rVec, self.vVec)
        self.LN = np.linalg.norm(self.LNVec)


MSUN_TSCALE = 4.927200467232548e-6
MSUN_RSCALE = 1.47662504e3
MPARSEC = 3.08567758e22

def make_unit_vec(theta, phi):
    return np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])

DEFAULT_V_ISCO = 1/np.sqrt(6)
DEFAULT_m1 = 30
DEFAULT_m2 = 30
DEFAULT_chi1 = 0.0
DEFAULT_chi1theta = 0.0
DEFAULT_chi1phi = 0.0
DEFAULT_chi2 = 0.0
DEFAULT_chi2theta = 0.0
DEFAULT_chi2phi = 0.0
DEFAULT_e0 = 0.0
DEFAULT_zeta0 = 0.0
DEFAULT_fmin = 20
DEFAULT_iota = 0.0
DEFAULT_dL = 100.0
DEFAULT_phic = 0.0
DEFAULT_kappa1 = 1.0
DEFAULT_kappa2 = 1.0
DEFAULT_vommax = 1./np.sqrt(6) # 0.277
class pypnprec(object):
    def __init__(self, **kwargs):
        self.__m1 = kwargs.get('m1', DEFAULT_m1)
        self.__m2 = kwargs.get('m2', DEFAULT_m2)
        self.__chi1 = kwargs.get('chi1', DEFAULT_chi1)
        self.__chi1theta = kwargs.get('chi1theta', DEFAULT_chi1theta)
        self.__chi1phi = kwargs.get('chi1phi', DEFAULT_chi1phi)
        self.__chi2 = kwargs.get('chi2', DEFAULT_chi2)
        self.__chi2theta = kwargs.get('chi2theta', DEFAULT_chi2theta)
        self.__chi2phi = kwargs.get('chi2phi', DEFAULT_chi2phi)
        self.__e0 = kwargs.get('e0', DEFAULT_e0)
        self.__fmin = kwargs.get('fmin', DEFAULT_fmin)
        self.__zeta0 = kwargs.get('zeta0', DEFAULT_zeta0)
        self.__Theta = kwargs.get('iota', DEFAULT_iota)
        self.__Phi = kwargs.get('phic', DEFAULT_phic)
        self.__dL = kwargs.get('distance', DEFAULT_dL)
    @property
    def hPref(self):
        return 2.*(self.__m1 + self.__m2)*MSUN_RSCALE*self.eta / (self.__dL*MPARSEC)
    def __str__(self):
        return f'<bbh: (m1,m2,chi1,chi2) = ({self.__m1}, {self.__m2}, {self.chi1Vec}, {self.chi2Vec})>'
    def __repr__(self):
        return self.__str__()
    @property
    def MT(self):
        return (self.m1+self.m2)*MSUN_TSCALE
    @property
    def e0(self):
        return self.__e0
    @property
    def m1(self):
        return self.__m1
    @property
    def m2(self):
        return self.__m2
    @property
    def eta(self):
        return self.__m1*self.__m2 / (self.__m1 + self.__m2)**2
    @property
    def dm(self):
        return np.sqrt(1. - 4.*self.eta)
    @property
    def X1(self):
        return 0.5*(1. + self.dm)
    @property
    def X2(self):
        return 0.5*(1. - self.dm)
    @property
    def chi1(self):
        return self.__chi1
    @property
    def chi2(self):
        return self.__chi2
    @property
    def s1(self):
        return self.__chi1*self.X1**2
    @property
    def s2(self):
        return self.__chi2*self.X2**2
    @property
    def chi1Vec(self):
        return self.__chi1*make_unit_vec(self.__chi1theta, self.__chi1phi)
    @property
    def chi2Vec(self):
        return self.__chi2*make_unit_vec(self.__chi2theta, self.__chi2phi)
    @property
    def s1Vec(self):
        return self.chi1Vec*self.X1**2
    @property
    def s2Vec(self):
        return self.chi2Vec*self.X2**2
    @property
    def Theta(self):
        return self.__Theta
    @property
    def Phi(self):
        return self.__Phi
    @property
    def sVec(self):
        return self.s1Vec + self.s2Vec
    @property
    def sigmaVec(self):
        return self.chi2Vec*self.X2 - self.chi1Vec*self.X1
    @property
    def chisVec(self):
        return 0.5*(self.chi1Vec + self.chi2Vec)
    @property
    def chiaVec(self):
        return 0.5*(self.chi1Vec - self.chi2Vec)
    @property
    def vmin(self):
        return np.power(np.pi * self.MT*self.__fmin, 1./3.)
    @property
    def fmin(self):
        return self.__fmin
    def Calc_v_form_DynVec(self):
        raise NotImplementedError('Calc_v_form_DynVec is not Implemented')
    
    def NDCalc_Prec_Newt_InitialCondition(self):
        raise NotImplementedError('NDCalc_Prec_Newt_InitialCondition is not implemented') 
    
    def NDCalc_Prec_Newt_DynWaveform(self, fmax:float = None, 
                                     t_bound = np.inf, 
                                     rtol=1e-13, 
                                     atol=1e-13, 
                                     max_step = np.inf):
        vmax = DEFAULT_vommax if fmax is None else np.power(np.pi * self.MT*fmax, 1./3.)
        def dydt(t, yy):
            x, y, z, vx, vy, vz, c1x, c1y, c1z, c2x, c2y, c2z = \
                yy[0], yy[1], yy[2], \
                yy[3], yy[4], yy[5], \
                yy[6], yy[7], yy[8]
            dynVec = pyPrecDynVec(self.eta, yy)
            raise NotImplementedError('dydt is not implemented')
        ICVec = self.NDCalc_Prec_Newt_InitialCondition()
        t0 = 0.0
        t = []
        y = []
        solver = RK45(dydt, t0, ICVec, t_bound = t_bound, rtol=rtol, atol=atol, max_step = max_step)
        while solver.status == 'running':
            t.append(solver.t)
            y.append(solver.y)
            if self.Calc_v_form_DynVec(solver.y) > vmax:
                break
            solver.step()
        tVec = np.asarray(t)
        DynArray = np.asarray(y)
        return pyPrecDyn(tVec, DynArray)
