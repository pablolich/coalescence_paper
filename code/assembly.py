#!/usr/bin/env python3

__appname__ = '[marsland_assembly.py]'
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

    #Environment with 50 strains and 50 metabolites
    s = 10
    m = 15 
    #Create a dictionary of parameters
    params = {
              'g':np.ones(s).reshape(s,1),
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
    kc = np.linspace(0, 0.9999,  10)
    kf = np.linspace(0, 1, num = 10)
    l = np.linspace(0.1, 0.9, num = 5)
    nsim = range(100)
    #Create N-D parameter grid 
    product = itertools.product(beta, kc, kf, l, nsim)
    #Create column names of data frame
    col_names = ['beta', 'kc', 'kf', 'l', 'n_sim']
    #Create dataframe for storing parameter values and simulation results
    df = pd.DataFrame(data = np.array(list(product)), columns = col_names)
    #Preallocate columns for C0 and F0 calculations
    ncol = len(df)
    df['C0'] = np.zeros(ncol)
    df['F0'] = np.zeros(ncol)
    #Preallocate columns for C and F calculations (after community assembly)
    df['C'] = np.zeros(ncol)
    df['F'] = np.zeros(ncol)
    #Preallocate column for richness
    df['r'] = np.zeros(ncol)
    #Run simulations across the parameter grid
    for i in progressbar.progressbar(range(ncol)):
        #Draw preference matrix and metabolic matrix
        c = preference_matrix(m, s, df['kc'][i],
                              beta = df['beta'][i], n_r = None)
        #Calculate demands as the sum of the rows of c
        demands = np.sum(c, axis = 0)
        #Sample metabolic matrix
        D = crossfeeding_matrix(demands, m, df['kf'][i])
        #Compute costs
        maint_vec = maintenance(c) 
        #Calculate competition and facilitation matrices of the community 
        #before the assembly
        C = interaction_matrix(df['l'][i], c, D, 
                               interaction = 'competition')
        F = interaction_matrix(df['l'][i], c, D, 
                               interaction = 'facilitation')
        #Average non-zero elements to get community level facilitation and 
        #competition
        df['C0'].iloc[i] = np.mean(C[np.triu_indices(len(C))])
        df['F0'].iloc[i] = np.mean(F[np.triu_indices(len(F))])
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
        #Get indices of extant species
        ind_extant = np.where(stable_abundance > 1)[0]
        #Store community richness
        df['r'].iloc[i] = len(ind_extant)
        #Get rid of these rows  in the matrix of preferences
        c_assembly = c[ind_extant,:]
        #Recalculate competition and facilitation community-level indices
        C = interaction_matrix(df['l'][i], c_assembly, D, 
                               interaction = 'competition')
        F = interaction_matrix(df['l'][i], c_assembly, D, 
                               interaction = 'facilitation')
        #Average non-zero elements to get community level facilitation and 
        #competition
        df['C'].iloc[i] = np.mean(C[np.triu_indices(len(C))])
        df['F'].iloc[i] = np.mean(F[np.triu_indices(len(F))])

        
    #Save results
    df.to_csv('../data/simulation_results.csv', index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
