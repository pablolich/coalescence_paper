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
from scipy.stats import pearsonr

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
              'K':20*np.ones(m).reshape(m,1),
              't':0.5*np.ones(m).reshape(m,1),
              'coal':0 #This is not a coalescence event
             }
    #Create time vector
    teval = np.arange(1, 5000)
    tspan = np.array([min(teval), max(teval)])
    #Set initial conditions
    abundances = np.ones(s)
    concentrations = 2*np.ones(m)
    z0 = list(abundances)+list(concentrations)
    #Set parameter values
    kc = 0
    kf = 0    
    l_B = np.linspace(start = 0.1, stop = 0.9, num = 100)
    l_A = np.array([0.1, 0.5, 0.9])
    Mc = class_matrix(m, m)
    nsim = 20
    #Preallocate vector of similarities
    S = np.zeros(shape = (len(l_A), nsim, len(l_B)))
    #Preallocate vector of facilitations
    Fav = np.zeros(shape = (len(l_A), nsim, len(l_B)))
    #Preallocate matrix of autocorrelation values
    autocorr = np.zeros(shape = (len(l_A)*nsim, len(l_B)))
    for j in range(len(l_A)):
        print('Leakage at assembly : ', j)
        #Loop over simulations
        for sim in range(nsim):
            print('Simulation number: ', sim)
            #Create vector of species classes
            sc_vec = np.repeat(np.arange(n_fam), n_memb)
            ##sample metablism (c and many D's at once)
            #spar = 0.05
            #c, D = sample_metabolism(m, s, kc, beta, kf, K, Mc, n_memb, l, spar)
            cA = preference_matrix(m, s, kc, beta = beta, n_r = None)
            cB = preference_matrix(m, s, kc, beta = beta, n_r = None)
            #Calculate cost of each species
            maint_vecA = maintenance(cA, l_A[j])
            #Compute demands as the sum of the rows of c
            demands_A = np.sum(cA, axis = 0)
            DA= general_metabolic_matrix(0.05, kf, demands_A, K, Mc, m, n_memb)
            demands_B = np.sum(cB, axis = 0)
            DB = general_metabolic_matrix(0.05, kf, demands_B, K, Mc, m, n_memb)
            #Instantiate an object a and b of class Community
            a = Community(kf, kc, s, cA, abundances, maint_vecA, 
                          np.ones(s)*l_A[j], DA, concentrations, 0)
            #Assembly
            kappa = 20*np.ones(m)
            A = a.assembly(kappa)
            #Get rid of extinct species
            A.n = np.delete(A.n, A.indext, axis = 0)
            #All resources are supplied
            supply = np.arange(A.m)
            ext_second_assembly = 0
            try:
                del n_vec
            except:
                pass
            for i in progressbar.progressbar(range(len(l_B))):
                maint_vecB = maintenance(cB, l_B[i], seed = sim+1)
                b = Community(kf, kc, s, cB, abundances, maint_vecB, 
                              np.ones(s)*l_B[i], DB, concentrations, 0)
                try:
                    n_vec_prev = n_vec
                except:
                    n_vec_prev = abundances
                B = b.assembly(kappa) 
                n_vec = B.n
                #Calculate autocorrelation between abundance vectors
                autocorr[nsim*j+sim, i] , _  = pearsonr(n_vec_prev, n_vec) 
                #Get rid of extinct species
                B.n = np.delete(B.n, B.indext, axis = 0)
                params, sol = A.coalescence(B, supply)
                abundance_mix = sol.y[0:params['s'],-1] 
                present_rows = np.where(abundance_mix > 1)[0]
                #Calculate facilitation matrix of community B
                Fmat = facilitation_matrix(B.l, B.D, B.c)
                Fmean = np.mean(Fmat)
                #Store
                Fav[j, sim, i] = Fmean
                #Create binary vectos of abundances of communities A and B in 
                #the space of the mix
                bool_A = np.hstack([np.ones(A.r, dtype = int), 
                                    np.zeros(B.r, dtype = int)])
                bool_B = np.hstack([np.zeros(A.r, dtype = int), 
                                    np.ones(B.r, dtype = int)])
                #Create binary vector of abundances after coalescence
                bool_AB = np.zeros(A.r + B.r, dtype = int)
                bool_AB[present_rows] = 1
                #Situate mixed community in the C1 C2 spectrum
                S[j, sim, i] = boolean_similarity(bool_AB, bool_A, bool_B)
    
    #Produce a data frame to save and analyze with R
    df = using_multiindex(S,'S', columns = list(['l_a', 'sim', 'l_b']))
    df_F = using_multiindex(Fav,'Fav', columns = list(['l_a', 'sim', 'l_b']))
    import ipdb; ipdb.set_trace(context = 20)
    df['Fav'] = df_F['Fav']
    df.to_csv('../data/effect_facilitation_new_cost.csv')
    df_autocorr = pd.DataFrame(autocorr)
    df_autocorr.to_csv('../data/autocorr.csv')
    #Average each group
    S_mean = np.zeros(shape = (len(S), len(l_B)))
    for i in range(len(S)):
        plt.plot(l_B, S[i,:,:].T, linewidth = 0.1, alpha = 0.3)
        S_mean[i,:] = np.mean(S[i, :, :], axis = 0)
        plt.plot(l_B, S_mean[i,:], linewidth = 1.5, 
                 label = r'$l_A = {}$'.format(l_A[i]))

    plt.xlabel('Leakage (l)')
    plt.ylabel('Similarity to prents ' + r'$(S_{1, 2})$')
    plt.legend()
    plt.show()

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

