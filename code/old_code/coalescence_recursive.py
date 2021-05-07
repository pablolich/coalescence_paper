#!/usr/bin/env python3

__appname__ = '[succesive_coalescence.py]'
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
                      community_competition, remove_all
from scipy.integrate import solve_ivp

## CONSTANTS ##

global seed; seed = 64

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

    #Number of resources
    m = len(D_mats.columns)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(np.array(assembly_data['l']))
    #Number of coalescence events for each leakage value
    n_sim = 100
    n_coal = 20
    #Initialize results storage dataframe
    product = itertools.product(l_vec, l_vec, np.arange(n_sim, dtype = int), 
                                np.arange(n_coal + 1, dtype = int))
    results = pd.DataFrame(data = np.array(list(product)),
                           columns = ['l', 'n_sim', 'n_coal'])
    ncol = len(results)
    results['Ftot'] = np.zeros(ncol) #Community total facilitation
    results['Fav'] = np.zeros(ncol) #Average community facilitation
    results['Ctot'] = np.zeros(ncol) #Community total competititon
    results['Cav'] = np.zeros(ncol) #Average community facilitation
    results['ext'] = np.zeros(ncol) #Community extinctions
    results['inv'] = np.zeros(ncol) #Community invasions
    results['R'] = np.zeros(ncol) #Community Richness (# of communities)
    results['r'] = np.zeros(ncol) #Species Richness (# number of species)
    results['path'] = np.zeros(ncol) #Average value of # pathways accross 
                                     #communities
    results['sd_met'] = np.zeros(ncol)
    #Perform coalescence for each value of l
    ijk = 0
    bar = progressbar.ProgressBar(max_value= ncol)
    for i in progressbar.progressbar(range(len(l_vec))):
        m = 0
        while m < n_sim:
            #Print status
            #Get data only with those levels of l
            data_i = assembly_data.loc[assembly_data['l'] == l_vec[i]]
            #Get abundance vectors of these communities
            abundances_i = abundances.iloc[data_i.index,:]
            #Select sample community  
            random.seed(seed)
            A_ind = int(np.random.choice(a = data_i.index, size = 1))
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
            xA = np.array(maintenance(cA_present))
            #Instantiate an object A of class Community
            A = Community(kf = data_i.loc[A_ind, 'kf'],
                          kc = data_i.loc[A_ind, 'kc'],
                          r = data_i.loc[A_ind, 'r'],
                          c = np.array(cA_present),
                          n = nA_present,
                          x = xA,
                          l = l_vec[i],
                          D = np.array(D_mats.iloc[D_mats.index == A_ind, :]),
                          n_coal = 1)
            #Calculate total facilitation
            F0, F0_i = community_facilitation(A.c, A.c, A.D, A.l)
            C0, C0_i = community_competition(A.c, A.c, A.l)
            results.loc[ijk, 'Ftot'] = F0 
            results.loc[ijk, 'Fav'] = F0_i 
            results.loc[ijk, 'Ctot'] = C0 
            results.loc[ijk, 'Cav'] = C0_i
            results.loc[ijk, 'ext'] = 0
            results.loc[ijk, 'inv'] = 0
            results.loc[ijk, 'R'] = 1
            results.loc[ijk, 'r'] = A.r
            results.loc[ijk, 'path'] = np.sum(A.c)
            results.loc[ijk, 'sd_met'] = np.std(A.D)
            ijk += 1
            bar.update(ijk)
            #Select communities to coalesce with the one 
            ind = np.random.choice(a = data_i.index, size = n_coal, 
                                   replace = False)
            #Start recursive coalescence
            n = 0
            all_comms = np.array([1])
            while n < n_coal:
                #Initialize community B
                #Abundance vector
                nBn = np.array(abundances_i.loc[ind[n],:])
                #Metabolic preference matrix
                cBn = c_mats.loc[ind[n]].reset_index(drop = True)
                #Get species abundance and preference vectors of present species only
                present_rows = np.where(nBn > 1)[0]
                #Number of species in community A
                sBn= len(present_rows)
                #Abundance vector of extant species
                nBn_present = nBn[present_rows]
                #Metabolic preference matrix of extant species
                cBn_present = cBn.iloc[present_rows,:]
                #Costs of extant species
                xBn = np.array(maintenance(cBn_present))
                #Initialite community B
                B = Community(kf = data_i.loc[ind[n], 'kf'],
                              kc = data_i.loc[ind[n], 'kc'],
                              r = data_i.loc[ind[n], 'r'],
                              c = np.array(cBn_present),
                              n = nBn_present,
                              x = xBn,
                              l = l_vec[i],
                              D = np.array(D_mats.iloc[D_mats.index == ind[n], :]),
                              n_coal = 1)
                all_comms = np.append(all_comms, n+2)
                #Choose the resources that will be supplied
                #supply = np.random.choice(a = np.arange(A.m),
                #                          size = np.random.randint(A.m))
                supply = np.arange(A.m)
                #Coalesce community A with community B
                params, sol = A.coalescence(B, supply)
                #Update community A
                abundance_mix = sol.y[0:params['s'],-1] 
                remaining = extinctions_one_two(abundance_mix, A.s)
                present_rows = np.where(abundance_mix > 1)[0]
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
                #Measure the total facilitation of the community
                F_tot, F_inter = community_facilitation(A.c, c_mix, A.D, A.l)
                C_tot, C_inter = community_competition(A.c, c_mix, A.l)
                results.loc[ijk, 'Ftot'] = F_tot
                results.loc[ijk, 'Fav'] = F_inter
                results.loc[ijk, 'Ctot'] = C_tot
                results.loc[ijk, 'Cav'] = C_inter
                results.loc[ijk, 'ext'] = len(extinct)
                results.loc[ijk, 'inv'] = com_after - com_before + len(extinct)
                results.loc[ijk, 'R'] = com_after
                results.loc[ijk, 'r'] = A.r
                results.loc[ijk, 'path'] = np.sum(A.c)
                results.loc[ijk, 'sd_met'] = np.std(all_D)
                ijk += 1
                bar.update(ijk)
                n += 1
            m += 1

    #Save results
    if len(sys.argv) > 1:
        results.to_csv('../data/recursive_coalescence_results'+sys.argv[1]+\
                       '.csv', index = False)
    else:
        df.to_csv('../data/recursive_coalescence_results.csv', index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
      
