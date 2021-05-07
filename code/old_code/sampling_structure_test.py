#Test dirichlet distribution sampling

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from functions_clean import class_matrix

## FUNCTIONS ##

def concentration_vector(kf, spar, M, alpha, m, Mc):
    '''
    Create concentration vector for each row (substrate)
    '''
    #Preallocate vector
    d = np.zeros(m)
    for prod in range(m):
        #Assign value to each element of the vector
        if M[alpha, prod]: 
            #When product and substrate belong to the same class assign a 
            #smaller value
            d[prod] = (1 - kf)/(spar*Mc)
        else:
            #When product and substrate belong to a different class assign a
            #smaller value
            d[prod] = 1/(spar*Mc)
        #Note that the values of d match independently of  the class of 
        #product and substrate when kf = 0.
    return d

def metabolic_matrix(kf, s, Mc, M, m):
    '''
    Generage vector of concentration parameters of a Dirichlet distribution

    Parameters: 
        kf (float): Facilitation constant
        s  (float): Sparsity parameter
        Mc (int): Number of individuals in each class
        M (mxm array): Matrix of classes
    '''
    n_it = len(M)
    #Preallocate metabolic matrix
    D = np.zeros(shape = (n_it, n_it))
    #Sample rows of D according to the concentration vector
    for i in range(n_it):
        #Sample concentration vector 
        c_vec = concentration_vector(kf, s, M, i, m, Mc)
        D[i,:] = np.random.dirichlet(c_vec)
    return D

def preferences_number(beta, m):
    '''
    Sample the number of preferences from an exponential distribution
    '''
    #Make sure that it is greater than 1, but smaller than the total number
    #of resources
    n = 0
    while n < 1 or n > m/3:
        n = round(np.random.exponential(beta))
    return n

def preferences_vector(Kc, m, nc, sc):
    '''
    Generate vector of metabolic preferences for species s 
    Parameters:
        Kc (float): Taxonomic heterogeneity constant
        m (int): Number of resources
        nc (int): Number of species per resource class
        sc (int): Class to which the species belongs
    '''
    #Preallocate preferences vector
    c_vec = np.zeros(m)
    #Get high and low probabilities
    high = 1/m*(1 + Kc*(m/nc - 1))
    low = 1/m*(1-Kc)
    #Sample number of preferences
    n = preferences_number(beta = 5, m = m)
    #What is the class of species s?
    prob = np.zeros(m)
    prob[nc*sc:nc*(sc+1)] = 1
    ind_0 = np.where(prob == 0)[0]
    prob[ind_0] = low
    ind_1 = np.where(prob == 1)[0]
    prob[ind_1] = high
    #Draw indices with probability p
    ind = np.random.choice(range(m), n, replace = False, p = prob)
    c_vec[ind] = 1
    return(c_vec)

def main(argv):
    '''Main function'''

    #Number of resources
    m = 60
    #Number of species
    s = 60
    #Number of classes
    c = 4
    #Get matrix of classes
    M = class_matrix(m, c)
    #Set value for facilitation constant
    kf = 0.5
    #Equal number of species in each class
    Mc = int(np.floor(s/c))
    #Set vectors of kf and kc
    kf_vec = np.array([0, 0.5, 0.8])
    n_kf = len(kf_vec)
    Kc_vec = np.array([0, 0.5, 0.8])
    n_kc = len(Kc_vec)
    fig, axs = plt.subplots(3, n_kc, figsize=(10, 6))
    #Sample metabolic matrix
    for i in range(n_kf):
        D = metabolic_matrix(kf_vec[i], 0.05, Mc, M, m)  
        np.savetxt('../data/D' + str(i) + '_struct.csv', D, delimiter = ',')
        axs[1][i].imshow(D)
        axs[1][i].set_title('Kf = ' + str(kf_vec[i]))
    #Preallocate matrix of preferences
    c_mat = np.zeros(shape = (s, m))
    #Create vector of classes
    class_vec = np.repeat(np.arange(c), Mc)
    #Sample metabolic preferences for all species
    for k in range(n_kc): 
        for i in range(s):
            c_mat[i,:] = preferences_vector(Kc_vec[k], m, nc = Mc, 
                                            sc = class_vec[i])

        np.savetxt('../data/c_mat' + str(k) + '_struct.csv', c_mat, delimiter = ',')
        axs[0][k].imshow(c_mat)
        axs[0][k].set_title('Kc = ' + str(Kc_vec[k]))
    plt.tight_layout()
    plt.show()

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

