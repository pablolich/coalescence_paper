#!/usr/bin/env python3

__appname__ = '[coalescence_clean.py]'
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
from functions_clean import joint_system, boolean_similarity, \
                            community_facilitation
from scipy.integrate import solve_ivp



## CONSTANTS ##

#global seed; seed = 69

## FUNCTIONS ##

def main(argv):
    '''Main function'''
    if len(sys.argv) > 1:
        #Load data from assembly simulations
        assembly_data = pd.read_csv('../data/simulation_results_'+sys.argv[1]+'.csv')
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
    #Get values of K (which are values for Kc and Kf)
    K_vec = np.unique(np.array(assembly_data['K']))
    K_vals = pd.DataFrame(np.triu_indices(len(K_vec))).T
    index_vec = np.array(K_vals.index)
    for struc in index_vec:
        results = pd.DataFrame(columns = ['l1', 'l2', 'S'])
        #Get vectors of parameters over which coalescence experiments will be 
        #run
        l_vec = np.unique(np.array(assembly_data['l']))
        for i in range(len(l_vec)):
            #Get values of K for each group 
            K1 = K_vec[K_vals.iloc[struc, 0]]
            K2 = K_vec[K_vals.iloc[struc, 1]]
            #Get data only with those levels of l
            data_iK = assembly_data.loc[(assembly_data['l'] == l_vec[i]) & \
                                        (assembly_data['K'] == K1)]
            #Get the other one
            data_ixl = assembly_data.loc[(assembly_data['l'] == l_vec[i]) & \
                                         (assembly_data['K'] == K2)]
            #Print message
            print('Coalescing communities for l = ', l_vec[i])
            #Get abundance vectors of these communities
            abundances_iK = abundances.iloc[data_iK.index,:]
            abundances_ixl = abundances.iloc[data_ixl.index,:]
            #Randomly sample n_sim pairs of communities to coalesce
            n_comm = len(data_iK)
            n_sim = 20000
            ind = np.random.choice(a = np.arange(n_comm), size = 2*n_sim)
            comb = ind.reshape(n_sim, 2)
            #Preallocate storage dataframe
            df = pd.DataFrame({'l1':l_vec[i]*np.ones(n_sim),
                               'l2':l_vec[i]*np.ones(n_sim),
                               'r1':np.ones(n_sim, dtype = int), #Richness of 1 
                               'r2':np.ones(n_sim, dtype = int), #"         " 2
                               'S':np.zeros(n_sim), #Similarity measure
                               'Fav1':np.zeros(n_sim), #Community fac 1
                               'Fav2':np.zeros(n_sim),  #Community fac 2
                               'Cav1':np.zeros(n_sim), #Community comp 1
                               'Cav2':np.zeros(n_sim), #Community comp 2
                               'Fmix':np.zeros(n_sim), #Facilitation mix
                               'F_looser':np.zeros(n_sim), #Facilitation 
                                                           #looser
                               'F_loos_ext':np.zeros(n_sim) #Facilitation of 
                                                            #looser and extinct
                               })
            #Perform all possible coalescence experiments between selected
            #communities.
            for j in progressbar.progressbar(range(n_sim)):
                #Get indices of coalescence communities 1, and 2
                ind_1 = data_iK.index[comb[j][0]]
                ind_2 = data_ixl.index[comb[j][1]]
                #Assign richness of each community
                df.loc[j, 'r1'] = data_iK.loc[ind_1, 'r']
                df.loc[j, 'r2'] = data_ixl.loc[ind_2, 'r']
                #Get competition and facilitation levels of communities 1 and 2
                Fav1 = data_iK.loc[ind_1, 'Fav']
                df.loc[j, 'Fav1'] = Fav1
                Cav1 = data_iK.loc[ind_1, 'Cav']
                Cbav1 = data_iK.loc[ind_1, 'Cbav']
                Ct1 = Cav1 + Cbav1
                df.loc[j, 'Cav1'] = Ct1
                Fav2 = data_ixl.loc[ind_2, 'Fav']
                df.loc[j, 'Fav2'] = Fav2
                Cav2 = data_ixl.loc[ind_2, 'Cav']
                Cbav2 = data_ixl.loc[ind_2, 'Cbav']
                Ct2 = Cav2 + Cbav2
                df.loc[j, 'Cav2'] = Ct2
                #Set parameters for community 1 
                n1 = np.array(abundances_iK.loc[ind_1,:])
                c1 = c_mats.loc[ind_1].reset_index(drop = True)
                #Get species abundance and preference vectors of present 
                #species only
                present_rows = np.where(n1 > 1)[0]
                s1 = len(present_rows)
                n1_present = n1[present_rows]
                c1_present = c1.iloc[present_rows,:]
                #Keep getting parameters...
                x1 = np.array(maintenance(c1_present))
                l1 = l_vec[i]*np.ones(m)
                D1 = D_mats.iloc[D_mats.index == ind_1, :]
                #Set parameters for community 2
                n2 = np.array(abundances_ixl.loc[ind_2,:])
                c2 = c_mats.loc[ind_2].reset_index(drop = True)
                #Get species abundance and preference vectors of present 
                #species only
                present_rows = np.where(n2 > 1)[0]
                s2 = len(present_rows)
                n2_present = n2[present_rows]
                c2_present = c2.iloc[present_rows,:]
                #Keep getting parameters...
                x2 = np.array(maintenance(c2_present))
                l2 = l_vec[i]*np.ones(m)
                D2 = D_mats.iloc[D_mats.index == ind_2, :]
                #Create joint system
                ext_system = joint_system(c1_present, D1, n1_present, l1, x1, 
                                          c2_present, D2, n2_present, l2, x2)
                #Get new demand profile due to being in the mix
                demands_tot = np.array(np.sum(c1_present, axis = 0) + \
                              np.sum(c2_present, axis = 0))
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
                z0 = (list(ext_system['N']) + #Abundances after assembly
                      list(2*np.ones(m))) #Fresh media
                #Create time vector
                tspan = tuple([1, 1e4])
                #Integrate diferential equations
                try:
                    sol = solve_ivp(lambda t,z: equations(t,z, params),
                                    tspan, z0,
                                    method = 'BDF', atol = 0.0001 )
                except:
                    continue
                #Get steady state abundance results
                stable_state = sol.y # length = r1 + r2 + m (r1, r2 < s)
                pop_abundance = stable_state[0:params['s'],-1]
                #Get vector of extinct species
                ext_sp = np.where(pop_abundance < 1)[0]
                #Get preference matrix of these species
                C_extinct = ext_system['C'][ext_sp]
                #Compute average facilitation of this group
                F_extinct = (community_facilitation(C_extinct, C_extinct,
                                                   ext_system['D'], 
                                                   l_vec[i], l_vec[i]))[1]
                #Store
                df.loc[j, 'Fext'] = F_extinct
                #Get index of present species after coalescence
                Cr = np.where(pop_abundance > 1, 1, 0)
                #Compute facilitation of the mixed community
                C_present = ext_system['C'][Cr]
                F_mix = (community_facilitation(C_present, C_present, 
                                                ext_system['D'],
                                                l_vec[i], l_vec[i]))[1]
                #Store
                df.loc[j, 'Fmix'] = F_mix
                #Create binary vector for the other abundances
                C1 = np.hstack([np.ones(s1, dtype = int), 
                                np.zeros(s2, dtype = int)])
                C2 = np.hstack([np.zeros(s1, dtype = int), 
                                np.ones(s2, dtype = int)])
                #Situate Cr in thee C1-C2 spectrum (compute similarity)
                S = boolean_similarity(Cr, C1, C2)
                #Does community 1 win?
                winner = S > 0
                if winner:
                    #Get presence absence vector of community 2 after coal
                    C2_pre = np.array(Cr[np.where(C2 !=0)[0]], dtype = bool)
                    #Reverse, to get absences instead of presences
                    C2_abs = np.invert(C2_pre)
                    #Add s1 preceding zeros to vector of presences
                    C2_pre = np.concatenate([np.zeros(s1, dtype = bool), 
                                             C2_pre])
                    #Add s1 preceding zeros to vector of absences
                    C2_abs = np.concatenate([np.zeros(s1, dtype = bool), 
                                             C2_abs])
                    #Select present species from the rows of C matrix
                    c2_surv = ext_system['C'][C2_pre]
                    #Select absent species from the rows of C matrix
                    c2_extinct = ext_system['C'][C2_abs]
                    #Record facilitation of community 2 (looser)
                    F_looser_surv = (community_facilitation(c2_surv, 
                                                            c2_surv, 
                                                            ext_system['D'],
                                                            l_vec[i],
                                                            l_vec[i]))[1]
                    #Record facilitation of extinctions in community 2 (looser)
                    F_looser_extinct = (community_facilitation(c2_extinct, 
                                                               c2_extinct, 
                                                               ext_system['D'],
                                                               l_vec[i],
                                                               l_vec[i]))[1]
                    #Store
                    df.loc[j, 'F_looser'] = F_looser_surv
                    df.loc[j, 'F_loos_ext'] = F_looser_extinct
                else:
                    #Get presence absence vector of community 1 after coal
                    C1_pre = np.array(Cr[np.where(C1 !=0)[0]], dtype = bool)
                    #Reverse, to get absences instead of presences
                    C1_abs = np.invert(C1_pre)
                    #Add s2 following zeros to vector of presences
                    C1_pre = np.concatenate([np.zeros(s2, dtype = bool), 
                                             C1_pre])
                    #Add s2 following zeros to vector of absences
                    C1_abs = np.concatenate([np.zeros(s2, dtype = bool), 
                                             C1_abs])
                    #Select present species from the rows of C matrix
                    c1_surv = ext_system['C'][C1_pre]
                    #Select absent species from the rows of C matrix
                    c1_extinct = ext_system['C'][C1_abs]
                    #Record facilitation of community 1 (looser)
                    F_looser_surv = (community_facilitation(c1_surv, 
                                                            c1_surv, 
                                                            ext_system['D'],
                                                            l_vec[i],
                                                            l_vec[i]))[1]
                    #Record facilitation of extinctions in community 1 (looser)
                    F_looser_extinct = (community_facilitation(c1_extinct, 
                                                               c1_extinct, 
                                                               ext_system['D'],
                                                               l_vec[i],
                                                               l_vec[i]))[1]
                    #Store
                    df.loc[j, 'F_looser'] = F_looser_surv
                    df.loc[j, 'F_loos_ext'] = F_looser_extinct

                #Get indices of surviving species
                ind_surv = np.where(pop_abundance > 1)[0]
                #Store similarities
                df.loc[j,'S'] = S
            #Add to the total dataframe 
            results = pd.concat([results, df])
        #Save dataframe for analysis
        if len(sys.argv) > 1:
            results.to_csv('../data/coalescence_clean_results_' + \
                           str(K1) + str(K2) + sys.argv[1] + '.csv', 
                           index = False)
        else:
            results.to_csv('../data/coalescence_clean_results_' + \
                           str(K1) + str(K2) + '.csv', index = False)

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
