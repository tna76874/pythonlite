#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
modules
"""
import sympy as sp
from sympy.core import evaluate
from sympy.utilities.lambdify import lambdify, implemented_function, lambdastr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
import sys
import random
from scipy.odr import *
import scipy.odr.odrpack as odrpack
import scipy.constants as cm
import pint as pn
import uncertainties
from uncertainties import ufloat, unumpy

import schemdraw
import schemdraw.elements as elm
from schemdraw.segments import *


# Definitions
## sympy
sf = sp.sympify

## units
ureg = pn.UnitRegistry(system='mks',autoconvert_offset_to_baseunit=True)
pe = ureg.parse_expression

ureg.define('fraction = [] = frac')
ureg.define('percent = 1e-2 frac = pct')
ureg.define('ppm = 1e-6 fraction')
ureg.default_format = "~P"

def uODR(func,beta0,xdata,ydata,**kwargs):
    """
    An ODR fit routine with standard deviations in ufloats
    Example:
    def func(p, x):
        a,b = p
        return a*x+b
    
    up , p, yfit, out = uODR(func,[0,0],X,Y)
    """
    xfit = kwargs.get('xfit',xdata)
    sx = kwargs.get('sx',None)
    sy = kwargs.get('sy',None)
    model = Model(func)
    data = RealData(xdata,ydata, sx=sx, sy=sy)
    odr = ODR(data, model, beta0=beta0)
    out = odr.run()
    # Estimated parameter values, of shape q
    p = out.beta
    # Standard errors of the estimated parameters, of shape p
    perr = out.sd_beta
    up =unumpy.uarray(p, perr)
    yfit=func(p, xfit) 
    return up , p, yfit, out

##### physical quantities
CM = pd.DataFrame(cm.physical_constants).T
CM['c'] = CM.index
CM = CM.reset_index(drop=True)

def getpc(quant):
    DF = CM
    DF = DF[DF['c']==quant].reset_index(drop=True)
    quan = DF[0][0] * ureg(DF[1][0])
    return quan

g = getpc('standard acceleration of gravity')
m_e = getpc('electron mass')
e = getpc('elementary charge')
c = getpc('speed of light in vacuum')
hp = getpc('Planck constant')
e0 = getpc('electric constant')
#####

