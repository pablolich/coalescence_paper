#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np

## CONSTANTS ##


## FUNCTIONS ##
def pairwise_competition(l, c_a, c_b):
    return((1 - l) * c_a @ c_b.T)

def pairwise_facilitation(D, l, c_a, c_b):
    return(l * c_a @ D @ c_b.T)

def interaction_matrix(l, c, D, interaction = 'both'):
    '''
    Calculate level of competition/facilitation/cohesion between all species 
    pairs
    '''
    #Get number of species and metabolites
    s = np.shape(c)[0]
    m = np.shape(c)[1]
    #Since the interaction matrix is symmetric, get the indices of the non-zero
    #elements (those on the upper triangle with offset of 1 over the main 
    #diagonal)
    ind = np.triu_indices(s, k = 1)
    #Prealocate vector of interaction
    interaction_vec = np.zeros(len(ind[0]))
    #Calculate interaction indices
    for i in range(len(ind[0])):
        #Get species indices
        ith = ind[0][i]
        jth = ind[1][i]
        #Find preferences of each species
        c_ith = c[ith,:]
        c_jth = c[jth,:]
        #Find shared preferences and normalize
        common = len(np.where(c_ith == c_jth)[0])/len(c_ith)
        if interaction == 'both':
            #Compute and store facilitation level
            interaction_vec[i] = pairwise_competition(l, c_ith, c_jth) - \
                                 pairwise_facilitation(D, l, c_ith, c_jth)
        elif interaction == 'competition':
            interaction_vec[i] = pairwise_competition(l, c_ith, c_jth) 
        elif interaction == 'facilitation':
            interaction_vec[i] = pairwise_facilitation(D, l, c_ith, c_jth)
        else:
            print('Unknown type of interaction')
            raise
    #Create cohesion matrix
    interaction_matrix = np.zeros(shape = (s, s))
    #Populate cohesion matrix
    interaction_matrix[ind] = interaction_vec
    return interaction_matrix

def joint_system(c_1, D_1, N_1, c_2, D_2, N_2, l1, l2, x1, x2):
    '''Create vectors and matrix of the joint system'''

    #Create off-diagonal zeros of joint system
    O = np.zeros( shape = (m, m), dtype = int )
    #Create matrices and vectors of the joint system
    D = np.asarray( np.bmat([[D_1, O], [O, D_2]]) )
    C = np.asarray( np.bmat([[c_1, O], [O, c_2]]) )
    N = np.vstack( [N_1, N_2] )
    l = np.vstack( [l1, l2] ) 
    x = np.vstack( [x1, x2] )
    #Create dictionary to return system
    di = {'D':D, 'C':C, 'N':N, 'l':l, 'x':x}
    return(di)

def competition_test(l, k1, k2):
    '''
    Test wether communities with higher kc have higher levels of total
    competition
    '''

    #Set parameters
    m = s = 50
    n_com = 100 
    n_k1 = len(k1)
    n_k2 = len(k2)
    #Prealocate vector of community-level competition
    C_ten = np.zeros(shape = (n_k1, n_com, n_k2))
    F_ten = np.zeros(shape = (n_k1, n_com, n_k2))
    #Generate preferences matrix with varying kc
    for k_1 in range(n_k1):
        for k_2 in range(n_k2):
            for sim in  range(n_com):
                #Sample a metabolic matrix
                c = preference_matrix(m, s, kc = k1[k_1],  beta = 5, n_r = None)
                C = interaction_matrix(l, c, D = 2, interaction = 'competition') 
                C_ten[k_1, sim, k_2] = np.mean(C[np.triu_indices(len(C))])
                #Calculate facilitation
                cum_demand = np.sum(c, 0)
                D = crossfeeding_matrix(cum_demand, m, k2[k_2])
                F = interaction_matrix(l, c, D = D, interaction = 'facilitation') 
                F_ten[k_1, sim,  k_2] = np.mean(F[np.triu_indices(len(F))])
    return(C_ten, F_ten)


##Test
#import matplotlib.pylab as plt
#from model import preference_matrix, crossfeeding_matrix
#n_k1 = 10
#n_k2 = 10
#k1 = np.linspace(0, 0.999999, n_k1)
#k2 = np.linspace(0, 0.999999, n_k2)
#C_ten, F_ten = competition_test(0.5, k1, k2)
##Average across simulations
#C_av = np.mean(C_ten, axis = 1)
#F_av = np.mean(F_ten, axis = 1)
##Plot
#for i in range(n_k1):
#    plt.scatter(C_av[i,:], F_av[i,:], s = 100, label = r'$k_c = %.2f$' %k1[i])
#plt.ylabel(r'$\overline{F_0}$', fontsize = 25)
#plt.xlabel(r'$\overline{C_0}$', fontsize = 25)
#plt.title(r'$k_c \in [0,\dots, 1) \quad k_f \in [0 ... 1)$', fontsize = 25)
#plt.legend(fontsize = 10)
#plt.show()
        
        
