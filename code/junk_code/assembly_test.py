#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from model import *
from functions import *
import numpy as np
import pandas as pd

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #First, integrate one community.
    s = 7
    m = 15
    l = 0.5
    n_sim = 1
    #Create a dictionary of parameters
    params = {'g':np.ones(s).reshape(s,1),
              's':s,
              'm':m,
              'K':20*np.ones(m).reshape(m,1),
              't':0.5*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    tspan = tuple([1, 1e4])
    #Set initial conditions
    z0 = list(np.ones(s))+list(2*np.ones(m))
    #Set identity matrix as consumer preference matrix
    c =  np.zeros(shape=(s, m))   
    np.fill_diagonal(c, np.ones(s))
    #Calculate demands as the sum of the rows of c
    demands = np.sum(c, axis = 0)
    #Sample metabolic matrix with kf = 0
    D = crossfeeding_matrix(demands, m, 0)
    #Compute costs
    maint_vec = maintenance(c) 
    #Calculate competition and facilitation matrices of the community 
    #before the assembly with a leakage factor of 0.5
    C = interaction_matrix(l, c, D, 
                           interaction = 'competition')
    F = interaction_matrix(l, c, D, 
                           interaction = 'facilitation')
    #Average non-zero elements to get community level facilitation and 
    #competition
    C0 = np.mean(C[np.triu_indices(len(C))])
    F0 = np.mean(F[np.triu_indices(len(F))])
    #Calculate cost of each species
    maint_vec = maintenance(c)
    #Add sampled strategies, metabolic, costs, and leakage
    #to the parameter dictionary
    params['D'] = D
    params['c'] = c
    params['x'] = maint_vec.reshape(s,1)
    params['l'] = l*np.ones(m).reshape(m,1)
    #Solve diferential equations
    sol = solve_ivp(lambda t,z: equations(t,z, params),
                    tspan, z0,
                    method = 'BDF', atol = 0.0001 )
    #Get abundance of species at stable state
    stable_abundance = sol.y[0:s,-1]
    #Integrate dynamics for the second community
    s = 8
    n_sim = 2
    #Create a dictionary of parameters
    params = {'g':np.ones(s).reshape(s,1),
              's':s,
              'm':m,
              'K':20*np.ones(m).reshape(m,1),
              't':0.5*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    tspan = tuple([1, 1e4])
    #Set initial conditions
    z0 = list(np.ones(s))+list(2*np.ones(m))
    #Set identity matrix as consumer preference matrix
    a = np.zeros(shape=(s, m-s))   
    b = np.zeros(shape=(s, s))
    np.fill_diagonal(b, np.ones(s))
    c1 = np.hstack([a, b])
    #Calculate demands as the sum of the rows of c
    demands = np.sum(c, axis = 0)
    #Sample metabolic matrix with kf = 0
    D1 = crossfeeding_matrix(demands, m, 0)
    #Compute costs
    maint_vec1 = maintenance(c1) 
    #Calculate competition and facilitation matrices of the community 
    #before the assembly with a leakage factor of 0.5
    C = interaction_matrix(l, c1, D1, 
                           interaction = 'competition')
    F = interaction_matrix(l, c1, D1, 
                           interaction = 'facilitation')
    #Average non-zero elements to get community level facilitation and 
    #competition
    C01 = np.mean(C[np.triu_indices(len(C))])
    F01 = np.mean(F[np.triu_indices(len(F))])
    #Calculate cost of each species
    maint_vec = maintenance(c)
    #Add sampled strategies, metabolic, costs, and leakage
    #to the parameter dictionary
    params['D'] = D1
    params['c'] = c1
    params['x'] = maint_vec1.reshape(s,1)
    params['l'] = l*np.ones(m).reshape(m,1)
    #Solve diferential equations
    sol1 = solve_ivp(lambda t,z: equations(t,z, params),
                    tspan, z0,
                    method = 'BDF', atol = 0.0001 )
    #Get abundance of species at stable state
    stable_abundance1 = sol1.y[0:s,-1]
    #Merge results from both communities
    #Cmatrices
    c_mats = pd.DataFrame(np.vstack([c, c1]), index=np.hstack([np.zeros(7), 
                                                               np.ones(8)]))
    #Dmatrices
    D_mats = pd.DataFrame(np.vstack([D, D1]), index=np.hstack([np.zeros(m), 
                                                               np.ones(m)]))
    #Abundances
    Ab_vecs = pd.DataFrame(np.vstack([np.hstack([stable_abundance, 0]),
                                      stable_abundance1]))
    #simulation results
    col_names = ['beta', 'kc', 'kf', 'l', 'n_sim', 'C0', 'F0', 'C', 'F', 'r']
    relleno = np.zeros(shape = (2, len(col_names)))
    #Create DataFrame
    df = pd.DataFrame(relleno, columns = col_names)
    #Populate dataframe
    df.iloc[0,:] = np.array([5, 0, 0, 0.5, 1, C0, F0, C0, F0, 7])
    df.iloc[1,:] = np.array([5, 0, 0, 0.5, 2, C01, F01, C01, F01, 8])
    #Save all results
    df.to_csv('../data/simulation_results_test.csv', index = False)
    c_mats.to_csv('../data/c_matrices_test.csv')
    D_mats.to_csv('../data/D_matrices_test.csv')
    Ab_vecs.to_csv('../data/abundances_test.csv')
    
    

    

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

