#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions_clean import *
from model import *
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Environment with s strains and m metabolites
    m = 60
    s = 60
    n_fam = 4
    n_memb = 15
    K = 0

    #Exponential rate to draw preferences
    beta = 5
    #Create a dictionary of parameters
    params = {'g':np.ones(s).reshape(s,1),
              's':s,
              'm':m,
              'K':2*np.ones(m).reshape(m,1),
              't':0.25*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    teval = np.arange(1, 500)
    tspan = np.array([min(teval), max(teval)])
    #Set initial conditions
    abundances = np.ones(s)
    concentrations = 2*np.ones(m)
    z0 = list(abundances)+list(concentrations)
    #Set parameter values
    kc = 0
    kf = 0    
    n_l = 2
    n_coal = 30
    l = np.linspace(0.1, 0.9, n_l)
    Mc = class_matrix(m, m)
    n_sim = 10
    #Create parameter values
    l_vec = np.linspace(0.1, 0.9, n_l)
    coal_vec = np.arange(n_coal, dtype = int)
    sim_vec = np.arange(n_sim, dtype = int)
    #Create N-D parameter grid 
    product = itertools.product(sim_vec, l_vec, l_vec, coal_vec)
    #Create column names of data frame
    col_names = ['n_sim', 'l1', 'l2', 'n_coal']
    #Create dataframe for storing parameter values and simulation results
    df = pd.DataFrame(data = np.array(list(product)), columns = col_names)
    pref_vec = np.arange(m, dtype = int)
    product_pref = itertools.product(sim_vec, l_vec, l_vec, coal_vec, pref_vec)
    col_names_pref = ['n_sim', 'l1', 'l2', 'n_coal', 'n_pref']
    df_pref = pd.DataFrame(data = np.array(list(product_pref)), 
                           columns = col_names_pref)
    #Add columns for storing numbers
    df['Cav'] = 0
    df['Fav'] = 0
    df_pref['N_n_pref'] = 0
    row = 0
    #Loop over simulations
    for i in range(n_sim):
        #Draw target community for lakage j
        for j in range(n_l):
            print('Sample metabolism of target community with l = ', l_vec[j])
            #Create vector of species classes
            sc_vec = np.repeat(np.arange(n_fam), n_memb)
            #Sample metabolism
            cA = preference_matrix(m, s, kc, beta = beta, n_r = None)
            #Calculate cost of each species
            maint_vecA = maintenance(cA, l_vec[j])
            #Compute demands as the sum of the rows of c
            demands_A = np.sum(cA, axis = 0)
            DA= general_metabolic_matrix(0.05, kf, demands_A, K, Mc, m, n_memb)
            for k in range(n_l):
                print('Assemble target community with l = ', l_vec[j])
                #Instantiate an object a and b of class Community
                a = Community(kf, kc, s, cA, abundances, maint_vecA, 
                              np.ones(s)*l_vec[j], DA, concentrations, 0)
                #Assembly target community
                kappa = 2*np.ones(m)
                T = a.assembly(kappa)
                #Get rid of extinct species
                T.n = np.delete(T.n, T.indext, axis = 0)
                n_bomb = 0
                all_comms = np.array([0])
                #Repeatedly bombard the target community with the weapon
                #communities
                print('Bombarding with communities of l = ', l_vec[k])
                for l in range(n_coal):
                    cB = preference_matrix(m, s, kc, beta = beta, n_r = None)
                    #Draw weapon communities for leakage k
                    demands_B = np.sum(cB, axis = 0)
                    DB = general_metabolic_matrix(0.05, kf, demands_B, K, Mc, 
                                                  m, n_memb)
                    maint_vecB = maintenance(cB, l_vec[k])
                    b = Community(kf, kc, s, cB, abundances, maint_vecB, 
                                  np.ones(s)*l_vec[k], DB, concentrations, 0)
                    W = b.assembly(kappa) 
                    #Get rid of extinct species
                    W.n = np.delete(W.n, W.indext, axis = 0)
                    supply = np.arange(m)
                    #Bombard A with B n_coal times
                    params, sol = T.coalescence(W, supply)
                    #Record average resource abundance
                    R = sol.y[params['s']:params['s'] + params['m'],-1]
                    #Plus one in number of bombardments
                    print(n_bomb, end = '\r')
                    n_bomb += 1
                    #Add community W to all_comms
                    all_comms = np.append(all_comms, n_bomb)
                    #Get species abundance and preference vectors of present 
                    #species only
                    abundance_mix = sol.y[0:params['s'],-1] 
                    present_rows = np.where(abundance_mix > 1)[0]
                    #Update community T after it has coalesced with W
                    T.s = np.copy(T.r)
                    T.r = len(present_rows)
                    T.c = params['c'][present_rows,:]
                    T.n = abundance_mix[abundance_mix > 1]
                    T.x = params['x'][present_rows]
                    T.l = params['l'][present_rows].reshape(1, T.r)[0]
                    T.source = np.hstack([T.source, np.repeat(T.coal, W.r)])
                    T.path = np.sum(T.c)
                    #Remove extinctions 
                    T.source = T.source[present_rows]
                    #Number of invasions
                    n_inv = len(np.where(T.source == T.coal)[0])
                    T.coal += 1
                    com_present = np.unique(T.source)
                    #Get missing communities
                    missing = list(set(all_comms).difference(com_present))
                    o, ind_missing, o = np.intersect1d(all_comms, missing, 
                                                       return_indices = True)
                    #Delete from allcomms
                    all_comms = np.delete(all_comms, ind_missing)
                    if ind_missing.size != 0:
                        #Get rid of D matrices of communities that have gone 
                        #compeltely extinct
                        params['D'] = remove_chunk(ind_missing, params['D'], 
                                                   T.m) 
                        #Get rid of c columns that are all 0
                        T.c = remove_stripe_revisited(T.c, T.m)
                    #Update D matrix in the mixed community
                    T.D = params['D']
                    #Obtain all the matrices D in the mix
                    all_D = remove_all(T.D, T.m)
                    #Get matrix of c in the off diagonals
                    c_mix = full_c(T, all_comms)
                    #Calculate competition and facilitation matrices
                    C = competition_matrix(c_mix)
                    F = facilitation_matrix(T.l, T.D, c_mix)
                    #df.loc[row, 'Cav'] = np.mean(C)
                    df.loc[row, 'Cav_i'] = (np.sum(C)-np.trace(C))/ \
                                           (C.size - len(C))
                    #df.loc[row, 'Fav'] = np.mean(F)
                    df.loc[row, 'Fav_i'] = (np.sum(F)-np.trace(F))/ \
                                           (F.size - len(F))
                    df.loc[row, 'n_inv'] = n_inv
                    df.loc[row, 'lav'] = np.mean(T.l)
                    df.loc[row, 'Rav'] = np.mean(R[abs(R) < 1])
                    df.loc[row, 'cost'] = np.mean(T.x)
                    if np.isnan(df.loc[row, 'Rav']):
                        df.loc[row, 'Rav'] = 0
                    df.loc[row, 'r'] = T.r
                    #Get abundances of group of species with n preferences 
                    group_abund = abundances_n_preferences(T.n, T.c, T.m)
                    df_pref.loc[m*row:m*(row+1)-1,'N_n_pref'] = group_abund/\
                                                            sum(group_abund)
                    
                    ##Add 1 to the row number
                    row += 1
    df.to_csv('../data/recursive_results.csv')
    df_pref.to_csv('../data/recursive_n_pref_results.csv')
    return 0

# CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

