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
from model import maintenance
from functions import joint_system


## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Load data from assembly simulations
    assembly_data = pd.read_csv('../data/simulation_results.csv')
    D_mats = pd.read_csv('../data/D_matrices.csv', index_col = 0)
    c_mats = pd.read_csv('../data/c_matrices.csv', index_col = 0)
    abundances = pd.read_csv('../data/abundances.csv', index_col = 0)
    #Number of resources
    m = len(D_mats.columns)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(np.array(assembly_data['l']))
    r_vec = np.unique(np.array(assembly_data['r']))
    #Create parameter grid 
    l_r_vec = np.array(list(itertools.product(l_vec, r_vec)))
    ncol = len(l_r_vec)
    #fraction of communities that will be coalesced
    fraction = 0.05
    for i in range(ncol):
        print('Coalescing simulations for l = ',l_r_vec[i][0], \
              'and richness = ', int(l_r_vec[i][1]))
        #Get rows in data frame have their l and r equal to i 
        data_i = assembly_data.loc[(assembly_data['l'] == l_r_vec[i][0]) & \
                                   (assembly_data['r'] == l_r_vec[i][1])]
        #Get closest integer to the fraction of coalescing communities
        n_comm = round(len(data_i)*fraction)
        #Get cohesion vector
        cohesion = np.array(data_i['F'] - data_i['C'])
        #Get community indices of those with top n_comm values of cohesion
        ind_top = cohesion.argsort()[-n_comm:][::-1]
        #Get community indices of those with lowest n_comm values of cohesion
        ind_low = cohesion.argsort()[:n_comm]
        #Concatenate these vectors
        ind_coal = np.hstack([ind_top, ind_low])
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
        comb = np.array(np.triu_indices(len(data_i.index), k = 1)).transpose()
        n_sim = len(comb)
        #Perform all possible coalescence experiments between selected
        #communities.
        for j in progressbar.progressbar(range(n_sim)):
            #Get indices of coalescence communities 1, and 2
            ind_1 = comb[j][0]
            ind_2 = comb[j][1]
            #if l_r_vec[i][1] == 5:
            #Set parameters for community 1 
            n1 = abundances.iloc[ind_1,:]
            c1 = (c_i.iloc[c_i.index == ind_1, :]).reset_index(drop = True)
            #Get species abundance and preference vectors of present species 
            #only
            present_rows = np.where(n1 > 1)[0]
            s1 = len(present_rows)
            n1_present = n1[present_rows]
            c1_present = c1.iloc[present_rows,:]
            #Keep getting parameters...
            x1 = maintenance(c1_present).reshape(s1, 1)
            l1 = l_r_vec[i][0]*np.ones(m).reshape(m, 1)
            D1 = D_i.iloc[D_i.index == ind_1, :]
            import ipdb; ipdb.set_trace(context = 20)
            #Set parameters for community 2
            n2 = abundances.iloc[ind_2,:]
            c2 = (c_i.iloc[c_i.index == ind_2, :]).reset_index(drop = True)
            #Get species abundance and preference vectors of present species 
            #only
            present_rows = np.where(n2 > 1)[0]
            s2 = len(present_rows)
            n2_present = n2[present_rows]
            c2_present = c2.iloc[present_rows,:]
            #Keep getting parameters...
            x2 = maintenance(c_2_present).reshape(s2, 1)
            l2 = l_r_vec[i][0]*np.ones(m).reshape(m, 1)
            D2 = D_i.iloc[D_i.index == ind_2, :]
            #Create joint system
            #full = joint_system(c_




            #Choose top and worst n_comm in terms of cohesion
            ##Get what indices of D_mats and c_mats are equal to these indices
            #data = assembly_data.iloc[ind, :]
            ##Set number of coalescing communities from each cohesion group
            #n = 10
            #ind_top_D = np.where(np.in1d(np.array(D_data.index), ind))[0]
            #ind_top_c = np.where(np.in1d(np.array(c_data.index), ind))[0]
            #ind_low_D = np.where(np.in1d(np.array(D_data.index), ind))[0]
            #ind_low_c = np.where(np.in1d(np.array(c_data.index), ind))[0]
            #import ipdb; ipdb.set_trace(context = 20)
            ##Concatenate to mask all of them in the dataframes
            #ind_coal = np.hstack([ind_top, ind_low])
            #ind_D = np.hstack([ind_top_D, ind_low_D]) 
            #ind_c = np.hstack([ind_top_c, ind_low_c]) 
            ##Filter data according with these indices
            #data_coal = data.iloc[ind_coal,:]
            #D_coal = D_data.iloc[ind_D, :]
            #c_coal = c_data.iloc[ind_c, :]


    
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

