#!/usr/bin/env python3

__appname__ = '[coalescence_test.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
from model import *
from functions import *

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Load data from assembly simulations
    assembly_data = pd.read_csv('../data/simulation_results_test.csv')
    #Load rest of data from simulations
    D_mats = pd.read_csv('../data/D_matrices_test.csv', 
                         index_col = 0)
    c_mats = pd.read_csv('../data/c_matrices_test.csv', 
                         index_col = 0)
    abundances = pd.read_csv('../data/abundances_test.csv', 
                             index_col = 0)
    m = len(D_mats.columns)
    #Preallocate storage dataframe
    l_vec = np.unique(assembly_data['l'])
    r_vec = np.unique(assembly_data['r'])
    n_sim = 1
    df = pd.DataFrame({'l1':l_vec[0]*np.ones(n_sim),
                       'l2':l_vec[0]*np.ones(n_sim),
                       'r1':r_vec[0]*np.ones(n_sim, dtype = int),
                       'r2':r_vec[1]*np.ones(n_sim, dtype = int),
                       'O1':np.zeros(n_sim), #Cohesion of community 1
                       'O2':np.zeros(n_sim), #"                   " 2
                       'S':np.zeros(n_sim) #Similarity (Cr, C1)
                      })
    #Get indices of coalescence communities 1, and 2
    ind_1 = 0 
    ind_2 = 1
    df.loc[0,'O1'] = assembly_data.loc[ind_1, 'F'] - \
                     assembly_data.loc[ind_1, 'C']
    df.loc[0,'O2'] = assembly_data.loc[ind_2, 'F'] - \
                     assembly_data.loc[ind_2, 'C']
    #Set parameters for community 1 
    n1 = np.array(abundances.loc[ind_1,:])
    c1 = c_mats.loc[ind_1].reset_index(drop = True)
    #Get species abundance and preference vectors of present species 
    #only
    present_rows = np.where(n1 > 1)[0]
    s1 = len(present_rows)
    n1_present = n1[present_rows]
    c1_present = c1.iloc[present_rows,:]
    #Keep getting parameters...
    x1 = np.array(maintenance(c1_present)).reshape(s1, 1)
    l1 = l_vec[0]*np.ones(m).reshape(m, 1)
    D1 = D_mats.iloc[D_mats.index == ind_1, :]
    #Set parameters for community 2
    n2 = np.array(abundances.loc[ind_2,:])
    c2 = c_mats.loc[ind_2].reset_index(drop = True)
    #Get species abundance and preference vectors of present species 
    #only
    present_rows = np.where(n2 > 1)[0]
    s2 = len(present_rows)
    n2_present = n2[present_rows]
    c2_present = c2.iloc[present_rows,:]
    #Keep getting parameters...
    x2 = np.array(maintenance(c2_present)).reshape(s2, 1)
    l2 = l_vec[0]*np.ones(m).reshape(m, 1)
    D2 = D_mats.iloc[D_mats.index == ind_2, :]
    #Create joint system
    ext_system = joint_system(c1_present, D1, n1_present, l1, x1, 
                              c2_present, D2, n2_present, l2, x2)
    #Create dictionary of parameters
    params = {'g':np.ones(s1 + s2).reshape(s1 + s2, 1),
              's':s1 + s2,
              'm':m,
              'K':20*np.ones(m).reshape(m, 1),
              't':0.5*np.ones(m).reshape(m, 1),
              'coal':1,#This is a coalescence event
              'D':ext_system['D'],
              'c':ext_system['C'],
              'x':ext_system['x'],
              'l':ext_system['l']
              }
    #Initial conditions
    z0 = list(ext_system['N']) + list(2*np.ones(m)) 
    #Create time vector
    tspan = tuple([1, 1e4])
    #Integrate diferential equations
    sol = solve_ivp(lambda t,z: equations(t,z, params),
                    tspan, z0,
                    method = 'BDF', atol = 0.0001 )
    #Get steady state abundance results
    stable_abundance = sol.y[0:params['s'],-1]
    #Concatenate preference matrices of c1 and c2
    c_both = np.vstack([c1_present, c2_present])
    #Get indices of surviving species
    ind_surv = np.where(stable_abundance > 1)[0]
    #Get preference matrix of mixed community
    c_mixed = c_both[ind_surv,:]
    #Get vector of species presence in the mixed space for C1 
    spec_1 = np.hstack([n1_present, np.zeros(s2)])
    #Calculate similarity
    S = np.dot(spec_1, stable_abundance)/ \
        (np.linalg.norm(spec_1)*np.linalg.norm(stable_abundance))
    #Store similarity
    df.loc[0,'S'] = S
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

