#!/usr/bin/env python3

__appname__ = '[sampling_no_structure.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import itertools
from functions import class_matrix

## FUNCTIONS ##

def rescaling(v):
    '''
    Rescaling vector v to be between 0 and 1
    '''

    return((v - min(v) * np.heaviside(-min(v), 0)) / \
           (sum(v - min(v) * np.heaviside(-min(v), 0) )))

def crossfeeding_matrix(demands, m, kf):
    '''
    Create metabolic matrix
    '''

    #Base level cross-feeding is given by a Dirichlet distribution with low
    #variance
    D_base = np.random.dirichlet(100*np.ones(m), m)
    #Base level cross-feeding is given by a uniform distribution. 
    #D_base = 1/m*np.ones(shape = (m, m))
    #Initialite matrix of metabolic cross-feeding based on demands
    D_dem = np.ones(shape = (m, m), dtype = int)
    #Create iterator with non-diagonal indices of matrix
    ind = itertools.permutations(np.arange(m), 2)
    for (i, j) in ind:
        D_dem[i, j] = (demands[j]-demands[i] ) \
        * np.heaviside((demands[j] - demands[i]), 0.5)
    
    #Get normalization constants for each column
    D_resc = (rescaling(D_dem.reshape(m*m, 1))).reshape(m,m)
    norm = np.sum(D_dem, axis = 1)
    D_norm = D_dem / norm[:, None]
    D_norm[np.isnan(D_norm)] = 0
    #Add base metabolic matrix and demand-based metabolic matrix weighted by k
    D = (1 - kf) * D_base + kf * D_norm
    #Rescale  rows of D so that they are all positive and add up to 1
    for i in range(m):
        D[i,:] = rescaling(D[i,:])
    return(D)

def metabolic_matrix_nostruct(kf, spar, demands, m, s):
    '''
    Generage vector of concentration parameters of a Dirichlet distribution

    Parameters: 
        kf (float): Facilitation constant
        s  (float): Sparsity parameter
        demands (1xm): Number of consumers of each resource
        m (int: Number of resources
    '''
    n_it = m 
    #Sample D according to the demands vector
    D = np.random.dirichlet((1 + kf*demands)/spar, size = m)
    return D

def sampling_probability(kc, cumulative_demand, m):
    '''
    Compute the sampling probability vector
    '''

    if not any(cumulative_demand):
        return(1/m*np.ones(m))
    else:
        base = (1 - kc) * 1 / m
        regime = kc * cumulative_demand/sum(cumulative_demand)
        #Transform nans to 0
        regime[np.isnan(regime)] = 0
        #Rescale to get probabilities
        score = base + regime
        #Rescale to get probabilities
        #probability =  (score - min(score) * np.heaviside(-min(score), 0)) / \
        #               (1 - m * min(score) * np.heaviside(-min(score), 0))
        return(score)

def preference_matrix(m, s, kc, beta, n_r):
    '''
    Create metabolic preferences matrix
    '''

    #Prealocate matrix
    c = np.zeros(shape = (s, m), dtype = int)
    for i in range(s):
        #Sample number of imported metabolites
        n = 0
        #All species having the same number of preferences
        if n_r:
            n = round(n_r)
        #Number of preferences follows an exponential distribution with rate
        #beta
        elif beta:
            #Make sure that it is greater than 1, but smaller than the total number
            #of resources
            while n < 1 or n > m/3:
                n = round(np.random.exponential(beta))
        #Compute the cumulative demand
        cum = np.sum(c[0:i,:], axis = 0)
        #Get sampling probability
        prob = sampling_probability(kc, cum, m)
        #Draw indices with probability p
        ind = np.random.choice(range(m), n, replace = False, p = prob)
        #Flip the ind positions in the vector of strategies for species sp
        c[i,ind] = 1
    return c

def main(argv):
    '''Main function'''

    #Number of resources
    m = 60
    #Number of species
    s = 60
    #Set vectors of kf and kc
    kf_vec = np.array([0.3, 0.5, 0.6])
    n_kf = len(kf_vec)
    Kc_vec = np.array([0.4, 0.5, 0.6])
    n_kc = len(Kc_vec)
    c_mat = preference_matrix(m, s, kf_vec[2], beta = 5, n_r = None)
    #Calculate demands as the sum of the rows of c
    demands = np.sum(c_mat, axis = 0)
    D = metabolic_matrix_nostruct(kf_vec[0], 0.05, demands, m, s)
    CD = c_mat@D
    plt.imshow(c_mat@D)
    np.savetxt('../data/effective_D.csv', CD, delimiter = ',')
    plt.show()
    #Plot
    fig, axs = plt.subplots(2, n_kc, figsize=(10, 6))
    #Sample metabolic preferences for all species
    #for k in range(n_kc): 
    #    c_mat = preference_matrix(m, s, kf_vec[k], beta = 5, n_r = None)
    #    #Calculate demands as the sum of the rows of c
    #    demands = np.sum(c_mat, axis = 0)
    #    D = metabolic_matrix_nostruct(kf_vec[k], 0.05, demands, m, s)
    #    np.savetxt('../data/c_mat' + str(k) + '.csv', c_mat, delimiter = ',')
    #    np.savetxt('../data/D_mat' + str(k) + '.csv', D, delimiter = ',')
    #    axs[1][k].imshow(D)
    #    axs[1][k].set_title('kf = ' + str(kf_vec[k]))
    #    axs[0][k].imshow(c_mat)
    #    axs[0][k].set_title('kc = ' + str(Kc_vec[k]))
    ##Sample metabolic matrix
    #plt.tight_layout()
    #plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
