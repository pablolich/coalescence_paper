#!/usr/bin/env python3

__appname__ = '[coalescence.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
import itertools
import progressbar
from model import maintenance, equations
from functions import joint_system
from scipy.integrate import solve_ivp


## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Load data from assembly simulations
    assembly_data = pd.read_csv('../data/simulation_results_sh.csv')
    D_mats = pd.read_csv('../data/D_matrices_sh.csv', index_col = 0)
    c_mats = pd.read_csv('../data/c_matrices_sh.csv', index_col = 0)
    abundances = pd.read_csv('../data/abundances_sh.csv', index_col = 0)
    #Number of resources
    m = len(D_mats.columns)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(np.array(assembly_data['l']))
    #Fraction of communities that will be coalesced
    fraction = 0.05
    for i in range(len(l_vec)):
        #Get data only with those levels of l
        data_i = assembly_data.loc[assembly_data['l'] == l_vec[i]]
        #Get vector of richnesses so I can iterate over it
        r_vec = np.unique(np.array(data_i['r']))
        for k in range(len(r_vec)):
            #Print message
            print('Coalescing simulations for l = ',l_vec[i], \
                  'and richness = ', int(r_vec[k]))
            #Get data with l and r values
            data_ik = data_i.loc[data_i['r'] == r_vec[k]]
            #Get rows in data frame have their l and r equal to i 
            #Get abundances of these communities
            abundances_i = abundances.iloc[data_ik.index,:]
            #Get closest integer to the fraction of coalescing communities
            n_comm = round(len(data_ik)*fraction)
            #Get cohesion vector
            cohesion = data_ik['F'] - data_ik['C']
            #Get community indices of those with top n_comm values of cohesion
            ind_top = cohesion.index[cohesion.argsort()[-n_comm:][::-1]]
            #Get community indices of those with lowest n_comm values of cohesion
            ind_low = cohesion.index[cohesion.argsort()[:n_comm]]
            #Concatenate these vectors
            ind_coal = np.hstack([np.array(ind_top), np.array(ind_low)])
            #Get indices corresponding to these communitie in the D_matrices and
            #c_matrices data frames
            D_index_pd = np.where(np.in1d(np.array(D_mats.index), ind_coal))[0]
            #Note that np.in1d(A, B) returns a boolean array indicating whether 
            #each value of A is found in B. np.where returns the indices of the 
            #True values.
            c_index_pd = np.where(np.in1d(np.array(c_mats.index), ind_coal))[0]    
            #Mask matrices with this indices
            D_i = D_mats.iloc[D_index_pd, :]
            c_i = c_mats.iloc[c_index_pd, :]
            #Get a vector of all possible pairwise combinations of selected
            #communities
            if r_vec[k] == 20:
                import ipdb; ipdb.set_trace(context = 20)
            comb = np.array(np.triu_indices(len(data_ik.index), k = 1)).transpose()
            n_sim = len(comb)
            #Perform all possible coalescence experiments between selected
            #communities.
            for j in progressbar.progressbar(range(n_sim)):
                #Get indices of coalescence communities 1, and 2
                ind_1 = data_ik.index[comb[j][0]]
                ind_2 = data_ik.index[comb[j][1]]
                #Set parameters for community 1 
                n1 = np.array(abundances_i.loc[ind_1,:])
                c1 = (c_i.iloc[c_i.index == ind_1, :]).reset_index(drop = True)
                #Get species abundance and preference vectors of present species 
                #only
                present_rows = np.where(n1 > 1)[0]
                s1 = len(present_rows)
                n1_present = n1[present_rows]
                c1_present = c1.iloc[present_rows,:]
                #Keep getting parameters...
                x1 = np.array(maintenance(c1_present)).reshape(s1, 1)
                l1 = l_vec[i]*np.ones(m).reshape(m, 1)
                D1 = D_i.iloc[D_i.index == ind_1, :]
                #Set parameters for community 2
                n2 = np.array(abundances_i.loc[ind_2,:])
                c2 = (c_i.iloc[c_i.index == ind_2, :]).reset_index(drop = True)
                #Get species abundance and preference vectors of present species 
                #only
                present_rows = np.where(n2 > 1)[0]
                s2 = len(present_rows)
                n2_present = n2[present_rows]
                c2_present = c2.iloc[present_rows,:]
                #Keep getting parameters...
                x2 = np.array(maintenance(c2_present)).reshape(s2, 1)
                l2 = l_vec[i]*np.ones(m).reshape(m, 1)
                D2 = D_i.iloc[D_i.index == ind_2, :]
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
                #Solve diferential equations
                sol = solve_ivp(lambda t,z: equations(t,z, params),
                                tspan, z0,
                                method = 'BDF', atol = 0.0001 )
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

