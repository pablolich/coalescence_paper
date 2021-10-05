
#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions_clean import *
from model import *
import matplotlib.pylab as plt
import progressbar

## CONSTANTS ##


## FUNCTIONS ##

def group_abundances(C, abundances):
    '''
    Sum abundances of those species with the same number of metabolic
    preferences

    Parameters:
        C (sxm): Metabolic matrix of abundances.
        abundances (1xs): Vector of species abundances at steady state.
    '''
    #Get maximum number of metabolic preferences
    m = len(C.T)
    #Get total abundance
    tot_abund = sum(abundances)
    #Get vector of number of preferences of each species
    n_r = np.sum(C, axis = 1)
    #Preallocate vector of abundances
    group_abund = np.zeros(m)
    for i in range(m):
        #Get indices of species with number of preferences i 
        ind_i = np.where(n_r == i+1)[0]
        #Get normalized abundances of this species
        group_abund[i] = sum(abundances[ind_i])/tot_abund
    return(group_abund)
        
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
              't':0.5*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    teval = np.arange(1, 500)
    tspan = np.array([min(teval), max(teval)])
    #Set initial conditions
    abundances = np.ones(s)
    concentrations = 2*np.ones(m)
    z0 = list(abundances)+list(concentrations)
    Mc = class_matrix(m, m)
    kf = 0    
    n_l = 1
    n_sim = 10000
    #n_kc = 5
    #Set parameter values
    kc_vec = np.array([0,  0.6, 0.8, 0.85])
    n_kc = len(kc_vec)
    #kc_vec = np.linspace(0, 0.9, num = n_kc)
    l_vec = np.linspace(start = 0.1, stop = 0.9, num = n_l)
    nr_vec = np.arange(m, dtype = int) + 1 
    sim_vec = np.arange(n_sim, dtype = int)
    #Create N-D parameter grid 
    product = itertools.product(sim_vec, l_vec, kc_vec, nr_vec)
    #Create column names of data frame
    col_names = ['sim', 'l', 'kc', 'n_pref']
    #Create dataframe for storing parameter values and simulation results
    df = pd.DataFrame(data = np.array(list(product)), columns = col_names)
    #Add columns for storing numbers
    #df['n0'] = 0
    #df['N0'] = 0
    #df['n'] = 0
    #df['N'] = 0
    df['N0nr'] = 0
    df['Nnr'] = 0
    df['r'] = 0
    df['rel_ab'] = 0
    for sim_k in progressbar.progressbar(sim_vec):
        for kc_i in range(n_kc):
            for l_j in range(n_l):
                #Create vector of species classes
                sc_vec = np.repeat(np.arange(n_fam), n_memb)
                #Sample metabolism
                cA = preference_matrix(m, s, kc = kc_vec[kc_i], beta = beta, 
                                       n_r = None)
                #Calculate cost of each species
                maint_vecA = maintenance(cA, l_vec[l_j])
                #Compute demands as the sum of the rows of c
                demands_A = np.sum(cA, axis = 0)
                DA = general_metabolic_matrix(0.05, kf, demands_A, K, Mc, m, 
                                              n_memb)
                #Instantiate an object a and b of class Community
                a = Community(kf, kc_vec[kc_i], s, cA, abundances, maint_vecA, 
                              np.ones(s)*l_vec[l_j], DA, concentrations, 0)
                #Get counts of each n preferences
                count = np.bincount(np.sum(cA, axis = 1))
                ii = range(len(count))
                #Create a frequency table
                freq_tab = np.vstack([ii, count[ii]]).T
                #Fill with zeros until length m is reached
                rest = m - (len(count)-1)
                fill = np.hstack([np.arange(len(count), m+1).reshape(rest, 1),
                                  np.zeros(shape = (rest, 1))])
                #Get complete table of frequencies
                freq_tab = np.vstack([freq_tab, fill])[1:, :]
                #Add frequencies to dataframe
                df.loc[(df['l'] == l_vec[l_j]) &
                       (df['kc'] == kc_vec[kc_i]) &
                       (df['sim'] == sim_vec[sim_k]), 'N0nr'] = freq_tab[:, 1]
                #Select complete set of count values
                #vals = freq_tab[:, 1]
                ##Select rows and add where appropriate
                #freqs = df.loc[(df['l'] == l_vec[l_j]) & \
                #               (df['kc'] == kc_vec[kc_i]), 'n0']
                #df.loc[(df['l'] == l_vec[l_j]) & \
                #       (df['kc'] == kc_vec[kc_i]), 'n0'] = freqs + vals
                ##Keep track of total number of preferences
                #df.loc[(df['l'] == l_vec[l_j]) & \
                #       (df['kc'] == kc_vec[kc_i]), 'N0'] = sum(freqs + vals)
                #Assembly community
                kappa = 2*np.ones(m)
                T = a.assembly(kappa)
                #Calculate normalize abundance of each group
                gr_ab_vec = group_abundances(T.c, T.n)
                #Store in dataframe
                df.loc[(df['l'] == l_vec[l_j]) &
                       (df['kc'] == kc_vec[kc_i]) &
                       (df['sim'] == sim_vec[sim_k]), 'rel_ab'] = gr_ab_vec
                #Store richness in dataframe
                df.loc[(df['l'] == l_vec[l_j]) &
                       (df['kc'] == kc_vec[kc_i]) &
                       (df['sim'] == sim_vec[sim_k]), 'r'] = T.r
                #Count number of consumers with n preferences after assembly
                #Get counts of each n preferences
                count = np.bincount(np.sum(T.c, axis = 1))
                ii = range(len(count))
                #Create a frequency table
                freq_tab = np.vstack([ii, count[ii]]).T
                #Fill with zeros until length m is reached
                rest = m - (len(count)-1)
                fill = np.hstack([np.arange(len(count), m+1).reshape(rest, 1),
                                  np.zeros(shape = (rest, 1))])
                #Get complete table of frequencies
                freq_tab = np.vstack([freq_tab, fill])[1:, :]
                #Add frequencies to dataframe
                df.loc[(df['l'] == l_vec[l_j]) &
                       (df['kc'] == kc_vec[kc_i]) &
                       (df['sim'] == sim_vec[sim_k]), 'Nnr'] = freq_tab[:, 1]

                ##Select complete set of count values
                #vals = freq_tab[:, 1]
                ##Select rows and add where appropriate
                #freqs = df.loc[(df['l'] == l_vec[l_j]) & \
                #               (df['kc'] == kc_vec[kc_i]), 'n']
                #df.loc[(df['l'] == l_vec[l_j]) & \
                #       (df['kc'] == kc_vec[kc_i]), 'n'] = freqs + vals
                ##Keep track of total number of preferences
                #df.loc[(df['l'] == l_vec[l_j]) & \
                #       (df['kc'] == kc_vec[kc_i]), 'N'] = sum(freqs + vals)

    df.to_csv('../data/n_preferences_before_after_new.csv')
    return 0

# CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

     


