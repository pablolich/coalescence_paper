#!/usr/bin/env python3

__appname__ = '[assembly.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from model import *
import itertools
from functions import *
import pandas as pd
import progressbar

## CONSTANTS ##

global seed; seed = 69

## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Environment with s strains and m metabolites
    m = 60
    s = 60
    #Partition families
    n_fam = 4
    n_memb = 15
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
    #Create parameter vectors
    beta = np.linspace(5, 5,  num = 1, dtype = int)
    nr = np.linspace(1, 5, 2, dtype = int)
    kc = np.linspace(0, 0.9, num = 3)
    kf = np.linspace(0.01, 0.99, num = 3)
    Kc = np.linspace(0.1, 0.9, num = 3)
    Kf = np.linspace(0.1, 0.9, num = 3)
    #l = np.array([0.5])
    l = np.array([0.2, 0.5, 0.9])
    #l = np.array([0.1, 0.2, 0.3, 0.4, 0.50,
    #              0.50, 0.6, 0.7, 0.8, 0.9])
    #l = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.50,
    #              0.50, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
    nsim = range(2)
    #Create N-D parameter grid 
    product = itertools.product(beta, kc, kf, Kc, Kf, l, nsim)
    #Create column names of data frame
    col_names = ['beta', 'kc', 'kf', 'Kc', 'Kf', 'l', 'n_sim']
    #Create dataframe for storing parameter values and simulation results
    df = pd.DataFrame(data = np.array(list(product)), columns = col_names)
    #Preallocate columns for C0 and F0 calculations
    ncol = len(df)
    df['C0'] = np.zeros(ncol)
    df['F0'] = np.zeros(ncol)
    df['C0tot'] = np.zeros(ncol)
    df['C0av'] = np.zeros(ncol)
    df['C0bav'] = np.zeros(ncol)
    df['C0b1av'] = np.zeros(ncol)
    df['F0tot'] = np.zeros(ncol)
    df['F0av'] = np.zeros(ncol)
    #Preallocate columns for C and F calculations (after community assembly)
    df['C'] = np.zeros(ncol)
    df['F'] = np.zeros(ncol)
    df['Ctot'] = np.zeros(ncol)
    df['Cav'] = np.zeros(ncol)
    df['Cbav'] = np.zeros(ncol)
    df['Cb1av'] = np.zeros(ncol)
    df['Ftot'] = np.zeros(ncol)
    df['Fav'] = np.zeros(ncol)
    #Preallocate column for richness
    df['r'] = np.zeros(ncol, dtype = int)
    #Preallocate harvest factors
    df['B'] = np.zeros(ncol)
    df['Fin'] = np.zeros(ncol)
    #Preallocate average abundance and average pathway number vectors 
    df['av_ab'] = np.zeros(ncol)
    df['av_path'] = np.zeros(ncol)
    df['std_dem'] = np.zeros(ncol)
    #Preallocate columns for preference and metabolic matrices
    all_c = pd.DataFrame(data = np.zeros(shape = (s*ncol, m + 1), dtype = int))
    #Note that I add one column at the end to specify wether that species is 
    #extinct or not after assembly.
    all_D = pd.DataFrame(data = np.zeros(shape = (m*ncol, m)))
    all_abundances = pd.DataFrame(data = np.zeros(shape = (ncol, s)))
    #Set indices to number each community
    all_c = all_c.set_index(np.repeat(np.arange(ncol), s))
    all_D = all_D.set_index(np.repeat(np.arange(ncol), m))
    #Check the signs
    plus = 0
    minus = 0
    #Run simulations across the parameter grid
    for i in progressbar.progressbar(range(ncol)):
        #c = preference_matrix(m, s, df['kc'][i],
        #                      beta = df['beta'][i], n_r = None)
        #c = preferences_matrix(df['kc'][i], m, n_memb, n_fam,  s)
        #Create vector of species classes
        sc_vec = np.repeat(np.arange(n_fam), n_memb)
        #Draw preference matrix and metabolic matrix
        c = gen_pref_matrix(0, df['kc'][i], m, n_memb, sc_vec, sc_vec)
        #Compute demands as the sum of the rows of c
        demands = np.sum(c, axis = 0)
        #Sample metabolic matrix
        Dant = crossfeeding_matrix(demands, m, df['kf'][i])
        #Get matrix of classes
        Mc = class_matrix(m, n_fam)
        D = metabolic_matrix(kf = df['kf'][i], s = 0.05, M = Mc, Mc = n_memb, 
                             m = m)
        #D = general_metabolic_matrix(0.05, df['kf'][i], demands, 
        #                             0, Mc, m, n_memb)
        #Store in dataframe
        all_D[m*i:m*(i+1)]= D
        #Compute costs
        maint_vec = maintenance(c) 
        #Calculate facilitation cycling
        F_cy = community_facilitation(c, c, D, df['l'][i], df['l'][i])
        C_cy = community_competition(c, c, D, df['l'][i], df['l'][i])
        #Calculate competition and facilitation matrices of the community 
        #before the assembly
        C = interaction_matrix(df['l'][i], c, D, 
                               interaction = 'competition')
        F = interaction_matrix(df['l'][i], c, D, 
                               interaction = 'facilitation')
        ##Check the frequency of Cii' > Fii'
        #check = (C - F).reshape(1, s**2)
        #plus += sum((check  > 0)[0])
        #minus += sum((check < 0)[0])
        #if i > 8000:
        #    import ipdb; ipdb.set_trace(context = 20)
        #Average non-zero elements to get community level facilitation and 
        #competition leaving out the 0 of the diagonal
        df.loc[i, 'C0'] = np.sum(C)/(np.size(C)-len(np.diag(C)))
        df.loc[i, 'F0'] = np.sum(F)/(np.size(F)-len(np.diag(F)))
        df.loc[i, 'F0tot'] = F_cy[0]
        df.loc[i, 'C0tot'] = C_cy[0]
        df.loc[i, 'F0av'] = F_cy[1]
        df.loc[i, 'C0av'] = C_cy[1]
        df.loc[i, 'C0bav'] = C_cy[2]
        #Calculate cost of each species
        maint_vec = maintenance(c)
        #Add sampled strategies, metabolic, costs, and leakage
        #to the parameter dictionary
        params['D'] = D
        params['c'] = c
        params['x'] = maint_vec.reshape(s,1)
        params['l'] = df['l'][i]*np.ones(m).reshape(m,1)
        #Solve diferential equations
        sol = solve_ivp(lambda t,z: equations(t,z, params),
                        tspan, z0,
                        method = 'BDF', atol = 0.0001 )
        #Get abundance of species at stable state
        stable_abundance = sol.y[0:s,-1]
        #Store in dataframe
        all_abundances.iloc[i] = stable_abundance
        #Get indices of extant species
        ind_extant = np.where(stable_abundance > 1)[0]
        #Create a vector of extant species
        extant = -1*np.ones(shape = (s, 1))
        #Flip to 1 the extant ones
        extant[ind_extant] = 1
        c_tot = np.hstack([c, extant])
        #Store in data frame
        all_c[s*i:s*(i+1)] = c_tot
        #Get average abundance of extant species
        av_abundance = np.mean(stable_abundance[ind_extant])
        #Store average abundance
        df.loc[i, 'av_ab'] = av_abundance
        #Store community richness
        df.loc[i, 'r'] = len(ind_extant)
        #Get rid of these rows  in the matrix of preferences
        c_assembly = c[ind_extant,:]
        #Get average number of pathways in the community
        av_pathways = np.mean(np.sum(c_assembly, axis = 1))
        #Store
        df.loc[i, 'av_path'] = av_pathways
        #Get rid of these rows  in the vector of costs
        maint_assembly = maint_vec[ind_extant]
        #Get total demand of each resource
        demands = np.sum(c_assembly, axis = 0)
        #Store its standard deviation
        df.loc[i, 'std_dem'] = np.std(demands)
        #Calculate average terms of the effective harvest
        B = potential_harvest(df['l'][i], c_assembly, demands)
        Fplus = facilitation_in(df['l'][i], c_assembly, None, D, demands)
        #Store in data frame
        df.loc[i,'B'] = B
        df.loc[i,'Fin'] = Fplus
        #Recalculate facilitation and competititon community-level indices
        F = interaction_matrix(df['l'][i], c_assembly, D, 
                               interaction = 'facilitation')
        C = interaction_matrix(df['l'][i], c_assembly, D, 
                               interaction = 'competition')
        #Recalculate facilitation cycling
        F_cy = community_facilitation(c_assembly, c_assembly, D, 
                                      df['l'][i], df['l'][i])
        C_cy = community_competition(c_assembly, c_assembly, D,
                                     df['l'][i], df['l'][i])
        #Average non-zero elements to get community level facilitation and 
        #competition leaving out the 0 of the diagonal
        df.loc[i, 'C'] = np.sum(C)/(np.size(C)-len(np.diag(C)))
        df.loc[i, 'F'] = np.sum(F)/(np.size(F)-len(np.diag(F)))
        df.loc[i, 'Ftot'] = F_cy[0]
        df.loc[i, 'Ctot'] = C_cy[0]
        df.loc[i, 'Fav'] = F_cy[0]
        df.loc[i, 'Cav'] = C_cy[1]
        df.loc[i, 'Cbav'] = C_cy[2]
        if df['l'][i] == 0.9:
            import ipdb; ipdb.set_trace(context = 20)

    #Save results
    if len(sys.argv) > 1:
        df.to_csv('../data/simulation_results_'+sys.argv[1]+'.csv',
                  index = False)
        all_c.to_csv('../data/c_matrices_'+sys.argv[1]+'.csv')
        all_D.to_csv('../data/D_matrices_'+sys.argv[1]+'.csv')
        all_abundances.to_csv('../data/abundances_'+sys.argv[1]+'.csv')
    else:
        df.to_csv('../data/simulation_results.csv', index = False)
        all_c.to_csv('../data/c_matrices.csv')
        all_D.to_csv('../data/D_matrices.csv')
        all_abundances.to_csv('../data/abundances.csv')
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

