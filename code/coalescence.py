#!/usr/bin/env python3

__appname__ = '[coalescence.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
import itertools
import warnings
import random
import progressbar
from model import maintenance, equations
from functions import joint_system
from scipy.integrate import solve_ivp


## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Load data from assembly simulations
    assembly_data = pd.read_csv('../data/simulation_results.csv')
    #Count number of communities for every richness value
    counts_richness = assembly_data['r'].value_counts()
    #Get richness that have more than 500 counts
    r_vec = np.array([20, 30, 37])#np.array(counts_richness[counts_richness > 500].index)
    #Load rest of data from simulations
    D_mats = pd.read_csv('../data/D_matrices.csv', index_col = 0)
    c_mats = pd.read_csv('../data/c_matrices.csv', index_col = 0)
    abundances = pd.read_csv('../data/abundances.csv', index_col = 0)
    #Number of resources
    m = len(D_mats.columns)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(assembly_data['l'])#np.unique(np.array(assembly_data['l']))
    #Fraction of communities that will be coalesced
    results = pd.DataFrame(columns = ['l1', 'l2', 'r1', 'r2', 'O1', 'O2', 'S'])
    for i in range(len(l_vec)):
        #Get data only with those levels of l
        data_i = assembly_data.loc[assembly_data['l'] == l_vec[i]]
        #Get vector of richnesses so I can iterate over it
        for k in range(len(r_vec)):
            #Print message
            print('Coalescing communities for l = ',l_vec[i], \
                  'and richness = ', int(r_vec[k]))
            #Get data with l and r values
            data_ik = data_i.loc[data_i['r'] == r_vec[k]]
            #Get abundance vectors of these communities
            abundances_i = abundances.iloc[data_ik.index,:]
            comb = np.array(np.triu_indices(len(data_ik), k = 1)).transpose()
            #Avoid computational intensity
            if len(comb) > 5000:
                #Sample community pairs
                ind = np.random.randint(0, len(comb), 5000)
                #Choose thos indices
                comb = comb[ind]
            n_sim = len(comb)
            #Preallocate storage dataframe
            df = pd.DataFrame({'l1':l_vec[i]*np.ones(n_sim),
                               'l2':l_vec[i]*np.ones(n_sim),
                               'r1':r_vec[k]*np.ones(n_sim, dtype = int),
                               'r2':r_vec[k]*np.ones(n_sim, dtype = int),
                               'O1':np.zeros(n_sim), #Cohesion of community 1
                               'O2':np.zeros(n_sim), #"                   " 2
                               'S':np.zeros(n_sim) #Similarity (Cr, C1)
                               })
            #Perform all possible coalescence experiments between selected
            #communities.
            for j in progressbar.progressbar(range(n_sim)):
                #Get indices of coalescence communities 1, and 2
                ind_1 = data_ik.index[comb[j][0]]
                ind_2 = data_ik.index[comb[j][1]]
                #Get cohesion levels for each of the coalescing communities
                df.loc[j,'O1'] = data_ik.loc[ind_1, 'F'] - \
                                 data_ik.loc[ind_1, 'C']
                df.loc[j,'O2'] = data_ik.loc[ind_2, 'F'] - \
                                 data_ik.loc[ind_2, 'C']
                #Set parameters for community 1 
                n1 = np.array(abundances_i.loc[ind_1,:])
                c1 = c_mats.loc[ind_1].reset_index(drop = True)
                #Get species abundance and preference vectors of present species 
                #only
                present_rows = np.where(n1 > 1)[0]
                s1 = len(present_rows)
                n1_present = n1[present_rows]
                c1_present = c1.iloc[present_rows,:]
                #Keep getting parameters...
                x1 = np.array(maintenance(c1_present)).reshape(s1, 1)
                l1 = l_vec[i]*np.ones(m).reshape(m, 1)
                D1 = D_mats.iloc[D_mats.index == ind_1, :]
                #Set parameters for community 2
                n2 = np.array(abundances_i.loc[ind_2,:])
                c2 = c_mats.loc[ind_2].reset_index(drop = True)
                #Get species abundance and preference vectors of present species 
                #only
                present_rows = np.where(n2 > 1)[0]
                s2 = len(present_rows)
                n2_present = n2[present_rows]
                c2_present = c2.iloc[present_rows,:]
                #Keep getting parameters...
                x2 = np.array(maintenance(c2_present)).reshape(s2, 1)
                l2 = l_vec[i]*np.ones(m).reshape(m, 1)
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
                df.loc[j,'S'] = S
            #Add to the total dataframe 
            results = pd.concat([results, df])
    #Save dataframe for analysis
    results.to_csv('../data/coalescence_results.csv', index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

