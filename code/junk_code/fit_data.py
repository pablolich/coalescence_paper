#!/usr/bin/env python3

__appname__ = '[fit_data.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from lmfit import Minimizer, Parameters, report_fit
import matplotlib.pylab as plt
import pandas as pd
import numpy as np

## CONSTANTS ##


## FUNCTIONS ##

data = pd.read_csv('../data/fit_data.csv')

params = Parameters()

params.add('a', value = 1.2,max = 10)
params.add('b', value = 0.58, max = 10)

def residuals(params, x, y):

    v = params.valuesdict()

    model = v['a']*np.exp(v['b']*x) + 1


    return model - y

def residuals_cubic(params, x, y):

    v = params.valuesdict()

    model = v['b']*x**2
    
    return  model - y

minner = Minimizer(residuals, params, fcn_args = ( data['x'], data['y'] )  )

fit = minner.minimize()

report_fit(fit)

#Plot
x = data['x']
y = data['y']
result = y + fit.residual
plt.scatter(x, y)
x_vec = np.linspace(0,1,1000)
y_vec = np.ones(len(x_vec))
residual_smooth = residuals(fit.params, x_vec, y_vec)

#Plot x transformation


plt.plot(x_vec, y_vec + residual_smooth)
plt.show()

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

