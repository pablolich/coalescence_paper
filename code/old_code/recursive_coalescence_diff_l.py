#!/usr/bin/env python3

__appname__ = '[App_name_here]'
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
from functions import joint_system, boolean_similarity, facilitation_in,\
                      Community, community_facilitation, extinctions_one_two, \
                      remove_chunk, full_c, remove_stripe, \
                      community_competition, remove_all, extinctions, \
                      cohesion_looser
from scipy.integrate import solve_ivp

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Load data
    if len(sys.argv) > 1:
        #Load data from assembly simulations
        assembly_data = pd.read_csv('../data/simulation_results_' + \
                                    sys.argv[1] + '.csv')
        #Load rest of data from simulations
        D_mats = pd.read_csv('../data/D_matrices_'+sys.argv[1]+'.csv', 
                             index_col = 0)
        c_mats = pd.read_csv('../data/c_matrices_'+sys.argv[1]+'.csv', 
                             index_col = 0)
        c_mats = c_mats.iloc[:,0:-1]
        abundances = pd.read_csv('../data/abundances_'+sys.argv[1]+'.csv', 
                                 index_col = 0)
    else:
        #Load data from assembly simulations
        assembly_data = pd.read_csv('../data/simulation_results.csv')
        #Load rest of data from simulations
        D_mats = pd.read_csv('../data/D_matrices.csv', 
                             index_col = 0)
        c_mats = pd.read_csv('../data/c_matrices.csv', 
                             index_col = 0)
        #Get rid of the last column
        c_mats = c_mats.iloc[:,0:-1]
        abundances = pd.read_csv('../data/abundances.csv', 
                                 index_col = 0)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(np.array(assembly_data['l']))
    #l_vec = np.insert(l_vec, 10, 0.5)
    #Number of coalescence events, replicates, and leakage values
    n_l = len(l_vec)
    n_sim = 20000
    n_coal = 1
    inter_l = False
    if inter_l:
        #Get all posible combination between parameter values
        product = itertools.product(np.arange(n_l, dtype = int), 
                                    np.arange(n_l, dtype = int), 
                                    np.arange(n_sim, dtype = int), 
                                    np.arange(n_coal, dtype = int))
        #Initialize results storage dataframe
        results = pd.DataFrame(data = np.array(list(product)),
                               columns = ['l1', 'l2', 'n_sim', 'n_coal'],
                               dtype = int)
    else:
        #Get all posible combination between parameter values
        product = itertools.product(np.arange(n_l, dtype = int), 
                                    np.arange(n_sim, dtype = int), 
                                    np.arange(n_coal, dtype = int))
        #Initialize results storage dataframe
        comb = np.array(list(product))
        ncol = np.shape(comb)[0]
        data = np.hstack([comb[:, 0].reshape(ncol, 1), comb])
        results = pd.DataFrame(data = data,
                               columns = ['l1', 'l2', 'n_sim', 'n_coal'],
                               dtype = int)
    #Create  vector indexing each row 
    index_vec = np.array(results.index)
    #Total number of rows in dataframe
    ncol = len(index_vec)
    #Preallocate
    results['lav'] = np.zeros(ncol) #Community average leakage
    results['F1'] = np.zeros(ncol) #Community 1 total facilitation
    results['F2'] = np.zeros(ncol) #Community 2 total facilitation
    results['Ftot'] = np.zeros(ncol) #Community mixed total facilitation
    results['Fav'] = np.zeros(ncol) #Average mixed community facilitation
    results['C1'] = np.zeros(ncol) #Community 1 total competition
    results['C2'] = np.zeros(ncol) #Community 2 total competition
    results['Ctot'] = np.zeros(ncol) #Community mixed total competititon
    results['Cav'] = np.zeros(ncol) #Average mixed community competition
    results['inv1'] = np.zeros(ncol) #Percent of invasions to community 1
    results['inv2'] = np.zeros(ncol) #Percent of invasions to community 2
    results['ext1'] = np.zeros(ncol) #Percent of extinctions in community 1
    results['ext2'] = np.zeros(ncol) #Percent of extinctions in community 2
    results['Ext'] = np.zeros(ncol) #Community extinctions
    results['Fext'] = np.zeros(ncol) #Facilitation of extinctions of looser 
    results['Cext'] = np.zeros(ncol) #Competition of extinctions of looser 
    results['Inv'] = np.zeros(ncol) #Community invasions
    results['r1'] = np.zeros(ncol) #Community 1 species-richness
    results['r2'] = np.zeros(ncol) #Community 2 species-richness
    results['r'] = np.zeros(ncol) #Mixed community species-richness 
    results['R'] = np.zeros(ncol) #Mixed community community-Richness 
    results['path'] = np.zeros(ncol) #Mixed community average "              "
    results['S'] = np.zeros(ncol) #Similarity C1-C2
    results['looser'] = np.zeros(ncol) #Loosing community (1, -1 or 0)
    #Get a dictionary where keys are column names and values are column numbers
    name_number = {c: i for i, c in enumerate(results.columns)}
    #Perform simulations across all parameter space
    for ind in progressbar.progressbar(index_vec):
        #Every n_l*n_sim*n_coal, change leakage of target communities
        if (ind%(n_l*n_sim*n_coal) == 0) & (inter_l == True):
            #Get target communities with the value of leakage corresponding to 
            #the simulation we are in
            l1 = l_vec[int(results.iloc[ind, name_number['l1']])]
            data_i = assembly_data.loc[assembly_data['l'] == l1]
            #Get abundance vectors of these communities
            abundances_i = abundances.iloc[data_i.index,:]
            #Store leakage of target communities
            results.iloc[ind:(ind+n_l*n_sim*n_coal), name_number['l1']] = l1
            #Select sample communities
            A_inds = np.random.choice(a = data_i.index, size = n_sim)
            #Every n_sim*n_coal change leakage of weapon communities
            if ind%(n_sim*n_coal) == 0:
                #Get weapon communities with the value of leakage corresponding to 
                #the simulation we are in
                l2 = l_vec[int(results.iloc[ind, name_number['l2']])]
                data_j = assembly_data.loc[assembly_data['l'] == l2]
                #Get abundance vectors of these communities
                abundances_j = abundances.iloc[data_j.index,:]
                #Store leakage of weapon communities
                results.iloc[ind:(ind+n_sim*n_coal), \
                             name_number['l2']] = l2
        elif (ind%(n_sim*n_coal) == 0) & (inter_l != True):
            #Get target communities with the value of leakage corresponding to 
            #the simulation we are in
            l1 = l_vec[int(results.iloc[ind, name_number['l1']])]
            data_i = assembly_data.loc[assembly_data['l'] == l1]
            #Get abundance vectors of these communities
            abundances_i = abundances.iloc[data_i.index,:]
            #Store leakage of target communities
            results.iloc[ind:(ind+n_sim*n_coal), name_number['l1']] = l1
            #Select sample communities
            A_inds = np.random.choice(a = data_i.index, size = n_sim)
            #Get weapon communities with the value of leakage corresponding to 
            #the simulation we are in
            l2 = l_vec[int(results.iloc[ind, name_number['l2']])]
            data_j = assembly_data.loc[assembly_data['l'] == l2]
            #Get abundance vectors of these communities
            abundances_j = abundances.iloc[data_j.index,:]
            #Store leakage of weapon communities
            results.iloc[ind:(ind+n_sim*n_coal), \
                         name_number['l2']] = l2
        #Every n_coal events, select new target community and sample n_coal 
        #weapon communities
        if ind%n_coal == 0:
            #Initialize vector of present communities in the mix
            all_comms = np.array([1])
            #Select communities to coalesce with the target 
            B_inds = np.random.choice(a = data_j.index, size = n_coal)
        #Start repeated bombardment
        #Select target community corresponding to the index we are in
        A_ind = A_inds[results.loc[ind, 'n_sim']]
        #Set parameters for A
        #Abundance vector
        nA = np.array(abundances_i.loc[A_ind,:])
        #Metabolic preference matrix
        cA = c_mats.loc[A_ind].reset_index(drop = True)
        #Get species abundance and preference vectors of present species only
        present_rows = np.where(nA > 1)[0]
        #Number of species in community A
        sA= len(present_rows)
        #Abundance vector of extant species
        nA_present = nA[present_rows]
        #Metabolic preference matrix of extant species
        cA_present = cA.iloc[present_rows,:]
        #Costs of extant species
        xA = np.array(maintenance(cA_present, l1))
        #Instantiate an object A of class Community
        A = Community(kf = data_i.loc[A_ind, 'kf'],
                      kc = data_i.loc[A_ind, 'kc'],
                      r = data_i.loc[A_ind, 'r'],
                      c = np.array(cA_present),
                      n = nA_present,
                      x = xA,
                      l = l1,
                      D = np.array(D_mats.iloc[D_mats.index == A_ind, :]),
                      n_coal = 1)
        #Get index of B_inds corresponding to the ind we are in
        B_ind = B_inds[results.loc[ind, 'n_coal']]
        #Avoid coalescence of the same community
        if B_ind == A_ind:
            continue
        #Initialize community B
        nBn = np.array(abundances_j.loc[B_ind,:])
        #Metabolic preference matrix
        cBn = c_mats.loc[B_ind].reset_index(drop = True)
        #Get species abundance and preference vectors of present species only
        present_rows = np.where(nBn > 1)[0]
        #Number of species in community A
        sBn= len(present_rows)
        #Abundance vector of extant species
        nBn_present = nBn[present_rows]
        #Metabolic preference matrix of extant species
        cBn_present = cBn.iloc[present_rows,:]
        #Costs of extant species
        xBn = np.array(maintenance(cBn_present, l2))
        #Initialite community B
        B = Community(kf = data_j.loc[B_ind, 'kf'],
                      kc = data_j.loc[B_ind, 'kc'],
                      r = data_j.loc[B_ind, 'r'],
                      c = np.array(cBn_present),
                      n = nBn_present,
                      x = xBn,
                      l = l2,
                      D = np.array(D_mats.loc[D_mats.index == B_ind]),
                      n_coal = 1)
        #Store in results the richness of the target comunity 
        results.iloc[ind:ind+n_coal, name_number['r1']] = A.r
        #Get species richness of the weapon communities
        weapon_r = np.array(data_j.loc[B_inds, 'r'])
        results.iloc[ind:ind+n_coal, name_number['r2']] = weapon_r
        #Store competition and facilitation of the sampled communities
        CA = data_i.loc[A_ind, 'Cav']
        results.iloc[ind:ind+n_coal, name_number['C1']] = CA
        CB = np.array(data_j.loc[B_inds, 'Cav'])
        results.iloc[ind:ind+n_coal, name_number['C2']] = CB
        FA = data_i.loc[A_ind, 'Fav']
        results.iloc[ind:ind+n_coal, name_number['F1']] = FA
        FB = np.array(data_j.loc[B_inds, 'Fav'])
        results.iloc[ind:ind+n_coal, name_number['F2']] = FB
        import ipdb; ipdb.set_trace(context = 20)
        all_comms = np.append(all_comms, results.loc[ind, 'n_coal'] + 2)
        #Choose the resources that will be supplied
        supply = np.arange(A.m)
        #Coalesce community A with community B
        params, sol = A.coalescence(B, supply)
        if not sol:
            #The solution has diverged
            import ipdb; ipdb.set_trace(context = 20)
            results.loc[ind, 'lav'] = None 
            break
        abundance_mix = sol.y[0:params['s'],-1] 
        remaining = extinctions_one_two(abundance_mix, A.s)
        present_rows = np.where(abundance_mix > 1)[0]
        #Create binary vectos of abundances of communities A and B in the 
        #space of the mix
        bool_A = np.hstack([np.ones(A.r, dtype = int), 
                            np.zeros(B.r, dtype = int)])
        bool_B = np.hstack([np.zeros(A.r, dtype = int), 
                            np.ones(B.r, dtype = int)])
        #Create binary vector of abundances after coalescence
        bool_AB = np.zeros(A.r + B.r, dtype = int)
        bool_AB[present_rows] = 1
        #Situate mixed community in the C1 C2 spectrum
        S = boolean_similarity(bool_AB, bool_A, bool_B)
        #Get number of extinctions in each community
        ext_comms = extinctions(bool_AB, bool_A, bool_B)
        #Determine looser community
        if S == 0:
            looser = 0
        else:
            looser =  -1*np.sign(S)
        F_looser = cohesion_looser(bool_AB, \
                                       A.r, A.D, A.c, A.l, \
                                       B.r, B.D, B.c, B.l, \
                                       looser)[0]
        C_looser = cohesion_looser(bool_AB, \
                                       A.r, A.D, A.c, A.l, \
                                       B.r, B.D, B.c, B.l, \
                                       looser)[1]
        #Update community A after it has coalesced with B
        A.s = np.copy(A.r)
        A.r = len(present_rows)
        A.c = params['c'][present_rows,:]
        A.n = abundance_mix[abundance_mix > 1]
        A.x = params['x'][present_rows]
        A.coal += 1
        com_before = len(np.unique(A.source))
        A.source = np.hstack([A.source, np.repeat(A.coal, B.r)])
        A.path = np.sum(A.c)
        #Remove extinctions 
        A.source = A.source[present_rows]
        com_present = np.unique(A.source)
        com_after = len(com_present)
        #Get intersection between present communities and all coalesced
        #communities
        o, o2, ind1 = np.intersect1d(A.source, all_comms, return_indices = True)
        #Get indices communities that are completely extinct
        extinct = np.delete(all_comms, ind1)
        #Remove the extinct communty from all_comm
        o, o1, ind2 = np.intersect1d(extinct, all_comms, return_indices = True)
        all_comms = np.delete(all_comms, ind2)
        if any(extinct):
            #Get rid of D matrices of communities that have gone compeltely 
            #extinct
            params['D'] = remove_chunk(ind2, params['D'], A.m) 
            #Get rid of c columns that are all 0
            A.c = remove_stripe(ind2, A.c, A.m)
        #Update D matrix in the mixed community
        A.D = params['D']
        #Obtain all the matrices D in the mix
        all_D = remove_all(A.D, A.m)
        #Get matrix of c in the off diagonals
        c_mix = full_c(A)
        weights = A.n
        values =  results.loc[ind, 'l2']*np.ones(len(A.source))
        original = np.where(A.source == 1)[0]
        values[original] = results.loc[ind, 'l1']
        #Get weighted average of l, by the abundance of each species
        lav = np.average(values, weights = weights)
        #Measure the total facilitation of the community
        F_tot, F_inter = community_facilitation(A.c, c_mix, A.D, A.l, A.l)
        C_tot, C_inter = community_competition(A.c, c_mix, A.l, A.l)
        #Store results after coalescence event
        results.loc[ind, 'Ftot'] = F_tot
        results.loc[ind, 'Fav'] = F_inter
        results.loc[ind, 'Ctot'] = C_tot
        results.loc[ind, 'Cav'] = C_inter
        results.loc[ind, 'Ext'] = len(extinct)
        results.loc[ind, 'Inv'] = com_after - com_before + len(extinct)
        results.loc[ind, 'lav'] = lav 
        results.loc[ind, 'r'] = A.r
        results.loc[ind, 'R']  = com_after
        results.loc[ind, 'path'] = np.sum(A.c)
        results.loc[ind, 'S'] = S
        results.loc[ind, 'ext1'] = ext_comms[0]
        results.loc[ind, 'ext2'] = ext_comms[1]
        results.loc[ind, 'looser'] = looser
        results.loc[ind, 'Fext'] = F_looser
        results.loc[ind, 'Cext'] = C_looser

    #Save results
    if len(sys.argv) > 1:
        results.to_csv('../data/recursive_diff_coalescence_results'+sys.argv[1]+\
                       '.csv', index = False)
    else:
        df.to_csv('../data/recursive_diff_coalescence_results.csv', index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

