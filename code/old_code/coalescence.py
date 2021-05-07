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
from functions import joint_system, boolean_similarity, facilitation_in, \
                      extinctions
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

    #Get values of kc and kf
    kc_vec = np.unique(assembly_data['kc'])
    kf_vec = np.unique(assembly_data['kf'])
    #Number of resources
    m = len(D_mats.columns)
    #Get vectors of parameters over which coalescence experiments will be run
    l_vec = np.unique(np.array(assembly_data['l']))
    #ind = [True,True,False,False, False, False, True]
    #l_vec = l_vec[ind]
    #Fraction of communities that will be coalesced
    results = pd.DataFrame(columns = ['l1', 'l2', 'r1', 'r2', 'O1', 'O2', 'S'])
    for i in range(len(l_vec)):
        #Get data only with those levels of l
        data_i = assembly_data.loc[assembly_data['l'] == l_vec[i]]
        ##Count number of communities for every richness value
        #counts_richness = data_i['r'].value_counts()
        #rmax = counts_richness.idxmax()
        #if rmax == 1:
        #    counts_richness = counts_richness[1:]
        #    rmax = counts_richness.idxmax()
        #Get vector of richnesses so I can iterate over it
        #r_vec = np.array([rmax])
        #Coalescence only between communities with most frequent richness
        #for k in range(len(r_vec)):
        #Print message
        print('Coalescing communities for l = ', l_vec[i])
        #Get data with l and r values
        #data_ik = data_i.loc[data_i['r'] == r_vec[0]]
        ##Eliminate some kc and kf values to avoid computational intensity
        #data_ik = data_ik.loc[((data_ik['kc'] == kc_vec[0]) |
        #                       (data_ik['kc'] == kc_vec[3]) |
        #                       (data_ik['kc'] == kc_vec[6])) &
        #                      ((data_ik['kf'] == kf_vec[0]) |
        #                       (data_ik['kf'] == kf_vec[3]) |
        #                       (data_ik['kf'] == kf_vec[6]))]
        #Get abundance vectors of these communities
        abundances_i = abundances.iloc[data_i.index,:]
        #Sample n_sim of communities to perform coalescence
        n_comm = len(data_i)
        n_sim = 15000
        ind = np.random.choice(a = np.arange(n_comm), size = 2*n_sim)
        comb = ind.reshape(n_sim, 2)
        #comb = np.array(np.triu_indices(len(data_ik), k = 1)).transpose()
        #Avoid computational intensity
        #if len(comb) > 5000:
        #    #Sample community pairs
        #    ind = np.random.randint(0, len(comb), 5000)
        #    #Choose thos indices
        #    comb = comb[ind]
        #n_sim = len(comb)
        #Preallocate storage dataframe
        df = pd.DataFrame({'l1':l_vec[i]*np.ones(n_sim),
                           'l2':l_vec[i]*np.ones(n_sim),
                           'r1':np.ones(n_sim, dtype = int), #Richness of 1 & 2
                           'r2':np.ones(n_sim, dtype = int),
                           'kc1':np.zeros(n_sim), #Constants 
                           'kf1':np.zeros(n_sim),
                           'kc2':np.zeros(n_sim),
                           'kf2':np.zeros(n_sim),
                           'Kf1':np.zeros(n_sim),
                           'Kc1':np.zeros(n_sim),
                           'Kf2':np.zeros(n_sim),
                           'Kc2':np.zeros(n_sim),
                           'C1':np.zeros(n_sim), #Competition levels
                           'F1':np.zeros(n_sim), #Facilitation levels
                           'C2':np.zeros(n_sim),
                           'F2':np.zeros(n_sim),
                           'O1':np.zeros(n_sim), #Cohesion of community 1
                           'O2':np.zeros(n_sim), #"                   " 2
                           'S':np.zeros(n_sim), 
                           'Fav1':np.zeros(n_sim), #Community facilitation 1
                           'Fav2':np.zeros(n_sim),  #Community facilitation 2
                           'Cav1':np.zeros(n_sim), #Community competition 1
                           'Cav2':np.zeros(n_sim) #Community competition 2
                           })
        #Perform all possible coalescence experiments between selected
        #communities.
        for j in progressbar.progressbar(range(n_sim)):
            #Get indices of coalescence communities 1, and 2
            ind_1 = data_i.index[comb[j][0]]
            ind_2 = data_i.index[comb[j][1]]
            #Assign richness of each community
            df.loc[j, 'r1'] = data_i.loc[ind_1, 'r']
            df.loc[j, 'r2'] = data_i.loc[ind_2, 'r']
            #Get competition and facilitation levels of communities 1 and 2
            Fav1 = data_i.loc[ind_1, 'Fav']
            df.loc[j, 'Fav1'] = Fav1
            Cav1 = data_i.loc[ind_1, 'Cav']
            Cbav1 = data_i.loc[ind_1, 'Cbav']
            Ct1 = Cav1 + Cbav1
            df.loc[j, 'Cav1'] = Ct1
            F2 = data_i.loc[ind_2, 'F']
            df.loc[j, 'F2'] = F2
            C2 = data_i.loc[ind_2, 'C']
            df.loc[j, 'C2'] = C2
            Fav2 = data_i.loc[ind_2, 'Fav']
            df.loc[j, 'Fav2'] = Fav2
            Cav2 = data_i.loc[ind_2, 'Cav']
            Cbav2 = data_i.loc[ind_2, 'Cbav']
            Ct2 = Cav2 + Cbav2
            df.loc[j, 'Cav2'] = Ct2
            #Get competition and facilitation factors of communities 1 and 2
            kf1 = data_i.loc[ind_1, 'kf']
            df.loc[j, 'kf1'] = kf1 
            kc1 = data_i.loc[ind_1, 'kc']
            df.loc[j, 'kc1'] = kc1
            kf2 = data_i.loc[ind_2, 'kf']
            df.loc[j, 'kf2'] = kf2
            kc2 = data_i.loc[ind_2, 'kc']
            df.loc[j, 'kc2'] = kc2
            #Get competition and facilitation factors of communities 1 and 2
            Kf1 = data_i.loc[ind_1, 'Kf']
            df.loc[j, 'Kf1'] = Kf1 
            Kc1 = data_i.loc[ind_1, 'Kc']
            df.loc[j, 'Kc1'] = Kc1
            Kf2 = data_i.loc[ind_2, 'Kf']
            df.loc[j, 'Kf2'] = Kf2
            Kc2 = data_i.loc[ind_2, 'Kc']
            df.loc[j, 'Kc2'] = Kc2
            #Get cohesion levels for each of the coalescing communities
            df.loc[j,'O1'] = F1 - C1 
            df.loc[j,'O2'] = F2 - C2 
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
            x1 = np.array(maintenance(c1_present))
            l1 = l_vec[i]*np.ones(m)
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
            x2 = np.array(maintenance(c2_present))
            l2 = l_vec[i]*np.ones(m)
            D2 = D_mats.iloc[D_mats.index == ind_2, :]
            #Create joint system
            ext_system = joint_system(c1_present, D1, n1_present, l1, x1, 
                                      c2_present, D2, n2_present, l2, x2)
            #Get new demand profile due to being in the mix
            demands_tot = np.array(np.sum(c1_present, axis = 0) + \
                          np.sum(c2_present, axis = 0))
            #Get facilitation of each community in the mix
            F1_mix = facilitation_in(df.loc[j, 'l1'], np.array(c1_present), 
                                     np.array(c2_present), np.array(D1), 
                                     demands_tot)
            F2_mix = facilitation_in(df.loc[j, 'l2'], np.array(c2_present), 
                                     np.array(c1_present), np.array(D2), 
                                     demands_tot)
            #Get inter-community facilitation terms
            F21 = facilitation_in(df.loc[j, 'l2'], np.array(c1_present), 
                                  np.array(c2_present), np.array(D2),
                                  demands_tot, inter = True)
            F12 = facilitation_in(df.loc[j, 'l1'], np.array(c2_present), 
                                  np.array(c1_present), np.array(D1),
                                  demands_tot, inter = True)
            #Get effective harvests of each community in the mix
            df.loc[j, 'H1'] = data_i.loc[ind_1, 'B'] + F1_mix + F21
            df.loc[j, 'H2'] = data_i.loc[ind_2, 'B'] + F2_mix + F12
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
            try:
                sol = solve_ivp(lambda t,z: equations(t,z, params),
                                tspan, z0,
                                method = 'BDF', atol = 0.0001 )
            except:
                continue
            #Get steady state abundance results
            stable_abundance = sol.y[0:params['s'],-1]
            #Create binary vector of abundaces
            Cr = np.where(stable_abundance > 1, 1, 0)
            #Create binary vector for the other abundances
            C1 = np.hstack([np.ones(s1, dtype = int), 
                            np.zeros(s2, dtype = int)])
            C2 = np.hstack([np.zeros(s1, dtype = int), 
                            np.ones(s2, dtype = int)])
            #Situate Cr in thee C1 C2 spectrum
            S = boolean_similarity(Cr, C1, C2)
            #Get indices of surviving species
            ind_surv = np.where(stable_abundance > 1)[0]
            #Get vector of species presence in the mixed space for C1 
            spec_1 = np.hstack([n1_present, np.zeros(s2)])
            #Calculate similarity
            S2 = np.dot(spec_1, stable_abundance)/ \
                (np.linalg.norm(spec_1)*np.linalg.norm(stable_abundance))
            #Store similarities
            df.loc[j,'S'] = S
            df.loc[j,'S2'] = S2
            ext_comms = extinctions(Cr, C1, C2)
            df.loc[j, 'ext1'] = ext_comms[0]
            df.loc[j, 'ext2'] = ext_comms[1]
        #Add to the total dataframe 
        results = pd.concat([results, df])
    #Save dataframe for analysis
    if len(sys.argv) > 1:
        results.to_csv('../data/coalescence_results_'+sys.argv[1]+'.csv',
                       index = False)
    else:
        results.to_csv('../data/coalescence_results.csv', index = False)


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

