#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from functions_clean import *
from model import *
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##
def make_feasible_R(s, m, l):
    '''
    Create a feasible matrix c, such that (1) it has an inverse, and (2) the 
    resources concentration at equilibrium are always positive
    '''
    
    #Initialize counter
    counter = 0
    #Create a non-singular c matrix
    singular = 1
    #Initialite negative vector of concentrations
    R_star = -1*np.ones(m)
    #Check that all resource concentrations are negative
    while any(R_star < 0):
        singular = 1
        while singular:
            #Sample a metabolic preferences matrix
            c = preferences_matrix(0, m, 1, m, s)
            try: 
                print(counter, end = '\r')
                counter += 1
                #Invert it
                c_inv = np.linalg.inv(c)
                singular = 0
                #Sample maintenance cost
                z = maintenance(c, l).reshape(m, 1)
                #Compute solution
                R_star = 1/(1-l)*c_inv @ z
            except: 
                pass
    return R_star.reshape(m)


def main(argv):
    #Initialize parameters
    s = 5
    m = 5
    n_l = 200
    l = np.linspace(0.1, 0.9, n_l)
    R_star_mat = np.zeros(shape = (n_l, m))
    for i in range(len(l)):
        R_star_mat[i,:] = make_feasible_R(s, m, l[i])
    #Save matrix of resources for each leakage
    df = using_multiindex(R_star_mat, 'R', columns = ['l', 'R_star'])
    df['l_val'] = l[df['l']]
    df = df.drop(['l', 'R_star'], axis = 1)
    df.to_csv('../data/feasible_space_new_cost.csv')
    #Average over resources
    R_star_vec = np.mean(R_star_mat, axis = 1)
    plt.plot(l, R_star_vec)
    plt.ylim(0, 5)
    plt.show()
        
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

