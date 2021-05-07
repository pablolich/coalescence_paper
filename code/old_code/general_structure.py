#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from functions import class_matrix

## FUNCTIONS ##

def constant_a(T, m, nc, kc, Tc, Kc):
    '''Return constant a'''
    return (Kc + nc/m*(1 - Kc))/((Kc + 1)*(nc/m*(1 - kc) + Tc/T*kc)) 

def prob_piecewise(Kc, kc, m, nc, sc, mc, dnorm, a):
    '''
    Generate probability vector of preference j for species s
    Parameters:
        Kc (float): Taxonomic heterogeneity constant
        kc (float): Competitive factor
        m (int): Total number of resources
        nc (int): Number of resources per class
        sc (int): Class to which species s belongs
        mc (int): Class to which resource j belongs
        dnorm (1xm): Normalized demand of resource j 
    '''
    if sc == mc:
        #resource is in the class to which the species belongs
        return ((1-kc)*1/m + dnorm*kc)*a*(1 + Kc)
    else:
        #resource is not in the class to which the species belongs
        return 1/m*(1-Kc)

def gen_prob_vec(Kc, kc, d, m, nc, sc, mc_vec):
    '''
    Compute sampling probability for the m preferences and species s
    Parameters:
        Kc (float): Taxonomic heterogeneity constant
        kc (float): Competitive factor
        d (1xm): Number of consumers consumer each resource
        m (int): Number of resources
        nc (int): Number of resources per class
        sc (int): Class to which species s belongs
        mc_vec (1xm): Vector of classes to which each resource belong
    '''
    #Calculate total number of demands
    T = np.sum(d)
    #Calculate demands within current class
    Tc = np.array([d[nc*i:nc*(i+1)].sum() for i in range(m//nc)])
    #Calculate normalization constants a and b
    a = constant_a(T, m, nc, kc, Tc, Kc)
    a = np.repeat(a, nc)
    #Preallocate vector of probabilities of each metabolite
    p_j = np.zeros(m)
    #Assign probability to each metabolite
    for j in range(m):
        p_j[j] = prob_piecewise(Kc, kc, m, nc, sc, mc_vec[j], d[j]/T, a[j])
    #Note that even though p_j sums to one, due to precision erros it doesn't
    #therefore, I normalize it.
    return p_j

def gen_pref_matrix(Kc, kc, m, nsp, sc_vec, mc_vec, beta = 3):
    '''
    Samplekpreference matrix
    Parameters:
        Kc (float): Taxonomic heterogeneity constant
        kc (float): Competitive factor
        m (int): Number of resources
        nsp (int): Number of resources per class
        sc_vec (1xm): Vector of classes to which each species belong
        mc_vec (1xm): Vector of classes to which each resource belong
    '''
    #Number of species
    s = len(sc_vec)
    #Preallocate matrix
    C = np.zeros(shape = (s, m))
    #Initialize vector of demands 
    d = np.ones(m)
    for alpha in range(s):
        #Sampling probability of each metabolite
        prob = gen_prob_vec(Kc, kc, d, m, nsp, sc_vec[alpha], mc_vec)
        #Sample number of consumed metabolites
        n = 0
        while n < 1 or n > m/3:
            n = round(np.random.exponential(beta))
        #Sample indices of the matrix
        ind = np.random.choice(range(m), n, replace = False, p = prob)
        C[alpha, ind] = 1
        #Calculate demands at iteration alpha
        d = np.sum(C, axis = 0)
    return C

def concentration_vector(kf, spar, M, d, Kf, m, alpha, nsp):
    '''
    Create concentration vector for each row (substrate)
    '''
    #Preallocate vector
    q = np.zeros(m) 
    for prod in range(m):
        #Assign value to each element of the vector
        if M[alpha, prod]: 
            #When product and substrate belong to the same class assign a 
            #smaller value
            q[prod] = (1 + kf*d[prod])*(1 - Kf)/(spar*nsp)
        else:
            #When product and substrate belong to a different class assign a
            #smaller value
            q[prod] = (1 + kf*d[prod])*(1 + Kf)/(spar*nsp)
        #Note that the values of d match independently of  the class of 
        #product and substrate when kf = 0.
    return q

def general_metabolic(s, kf, d, Kf, M, m, nsp):
    '''
    Generage vector of concentration parameters of a Dirichlet distribution

    Parameters: 
        kf (float): Facilitation constant
        s  (float): Sparsity parameter
        nsp (int): Number of individuals in each class
        M (mxm array): Matrix of classes
    '''
    n_it = len(M)
    #Preallocate metabolic matrix
    D = np.zeros(shape = (n_it, n_it))
    #Sample rows of D according to the concentration vector
    for i in range(n_it):
        #Sample concentration vector 
        c_vec = concentration_vector(kf, s, M, d, Kf, m, i, nsp)
        D[i,:] = np.random.dirichlet(c_vec)
    return D


def main(argv):
    '''Main function'''

    #Total number of metabolites
    m = 60
    #Number of classes
    n_class = 4
    #Get matrix of classes
    M = class_matrix(m, n_class)
    #Species per class
    nsp = m//n_class
    #Create vector of species classes
    sc_vec = np.repeat(np.arange(n_class), nsp)
    #gen_pref_matrix(Kc, kc, m, nsp, sc_vec, mc_vec, beta = 3)
    C = gen_pref_matrix(0.9, 0.1, m, nsp, sc_vec, sc_vec)
    d = np.sum(C, axis = 0)
    #general_metabolic(s, kf, d, Kf, M, m, nsp)
    D = general_metabolic(0.05, 0.5, d, 0.9, M, m, nsp)
    #both
    CD = C@D
    np.savetxt('../data/effective_D_struct.csv', CD, delimiter = ',')
    plt.imshow(CD)
    plt.show()

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

