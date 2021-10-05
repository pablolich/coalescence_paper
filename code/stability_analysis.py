#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions_clean import *
from model import *
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Environment with s strains and m metabolites
    m = 60
    s = 60
    n_fam = 4
    n_memb = 15
    K = 0

    #Exponential rate to draw preferences
    beta = 5
    #Create a dictionary of parameters
    params = {'g':np.ones(s).reshape(s,1),
              's':s,
              'm':m,
              'K':20*np.ones(m).reshape(m,1),
              't':0.5*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    teval = np.arange(1, 500)
    tspan = np.array([min(teval), max(teval)])
    #Set initial conditions
    abundances = np.ones(s)
    concentrations = 2*np.ones(m)
    z0 = list(abundances)+list(concentrations)
    #Set parameter values
    kc = 0
    kf = 0    
    l_A = np.array([0.1, 0.9])
    Mc = class_matrix(m, m)
    nsim = 5
    for j in range(len(l_A)):
        for i in range(nsim):
            #Sample metablism (c and many D's at once)
            cA = preference_matrix(m, s, kc, beta = beta, n_r = None)
            maint_vecA = maintenance(cA, l_A[j])
            #Compute demands as the sum of the rows of c
            demands_A = np.sum(cA, axis = 0)
            DA= general_metabolic_matrix(0.05, kf, demands_A, K, Mc, m, n_memb)
            l_A_vec = l_A[j]*np.ones(s).reshape(s, 1)
            a = Community(kf, kc, s, cA, abundances, maint_vecA, l_A_vec,
                          DA, concentrations, 0)
            #Assembly
            kappa = 20*np.ones(m)
            A = a.assembly(kappa)
            a.n = np.delete(a.n, a.indext)
            Jac = jacobian(l = l_A[j], g = 1, n = a.n.reshape(a.r, 1),
                           C =  a.c, R = a.R.reshape(a.m, 1), D = a.D, 
                           tau = 0.5*np.ones(m).reshape(m,1))
            eig_vals = np.linalg.eigvals(Jac)
            plt.scatter(eig_vals.real, eig_vals.imag)
            plt.show()

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
