#!/usr/bin/env python3

__appname__ = '[functions_cleaned.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from model import equations
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

## CONSTANTS ##


## FUNCTIONS ##

def joint_system(c_1, D_1, N_1, l1, x1, c_2, D_2, N_2, l2, x2):
    '''Create vectors and matrix of the joint system'''
    #Get number of metabolites and species in each community
    m = np.shape(c_2)[1]
    s1 = len(N_1)
    s2 = len(N_2)
    #Create off-diagonal zeros of joint system
    dO1 = np.zeros( shape = (np.shape(D_1)[0], np.shape(D_2)[1]), dtype = int )
    dO2 = np.zeros( shape = (np.shape(D_2)[0], np.shape(D_1)[1]), dtype = int )
    cO1 = np.zeros( shape = (np.shape(c_1)[0], np.shape(c_2)[1]), dtype = int )
    cO2 = np.zeros( shape = (np.shape(c_2)[0], np.shape(c_1)[1]), dtype = int )
    #Create matrices and vectors of the joint system
    D = np.asarray( np.bmat([[D_1, dO1], [dO2, D_2]]) )
    C = np.asarray( np.bmat([[c_1, cO1], [cO2, c_2]]) )
    N = np.hstack( [N_1, N_2] ) #These will be reshaped inside model equations
    l = np.hstack([l1, l2]).reshape(s1+s2, 1)
    x = np.vstack( [x1.reshape(s1, 1), x2.reshape(s2, 1)] )
    #Create dictionary to return system
    system = {'D':D, 'C':C, 'N':N, 'l':l, 'x':x}
    return(system)

def boolean_similarity(v, v1, v2):
    '''
    Calculate how more (less) similar is vector v to v1 than to v2

    Parameters: 
        v (1:s1+s2): Vector of absence-precence (1-0) values for species in the
                     resulting community
        v1 (1:s1+s2): Vector of absence-precnce values for species in 
                      community 1 with s2 zeros at the end
        v2 (1:s2+s1): Vector of absence-precnce values for species in 
                      community 2 with s1 preceeding zeros
    '''

    #Calculate normalized similarities of v with original communities
    s1 = np.dot(v, v1)/sum(v1)
    s2 = np.dot(v, v2)/sum(v2)
    #Calculate the difference between these similarities
    return(s2-s1)

def community_facilitation(c1, c2,  D, l1, l2):
    '''
    Compute the community-level facilitation, i.e., the total amount of energy
    that is leaked to the environment and captured by other species (or the 
    same one)
    '''
    
    #Calculate the amount of facilitation in the community
    F = c1 @ D @ np.transpose(c2)
    #Merge
    both = np.array([l1*(np.sum(F)-np.trace(F))/(F.size - len(F)), 
                     l1*np.mean(F),
                     l1*(1-l1)*np.mean(F)])
    return (both)

def competition_matrix(C):
    return(C @ np.transpose(C))

def facilitation_matrix(l, D, C):
    return( np.diag(l) @ C @ D @ C.T ) 

def community_competition(c1, c2, D, l1, l2):
    '''
    Compute the community/communities-level competition, i.e., the total amount
    of energy that is overlappingly required by all species pairs in the 
    community/mix
    '''
    shap1 = np.shape(c1)
    shap2 = np.shape(c2)
    #Calculate amount of competition in the community
    C = c1@ np.transpose(c2)
    #Get upper diagonal (diag included) indices of matrix of competition.
    #Note that these are the only indices that we want because Cb is symmetric
    ind = np.triu_indices(np.shape(C)[0], k = 1)
    n_sp = len(ind[0])
    #Preallocate matrix of biotic competition
    Cb = np.zeros(shape = np.shape(C))
    Cb1 = np.zeros(shape = np.shape(C))
    for i in range(n_sp):
        a = ind[0][i]
        b = ind[1][i]
        c_a = c1[a,:]
        c_b = c2[b,:]
        Sab = c_a + c_b
        Pab  = c_a * c_b
        Cb[a, b] = Sab @ D @ np.transpose(Pab)
    #Merge
    both = np.array([np.sum(C[ind]), 
                     np.mean((1-l1)*C[ind]),
                     np.mean(l1*Cb[ind]),
                     np.mean((1-l1)**2*C[ind])])
    return(both)

def class_matrix(m, c):
    '''
    Create a square matrix with 1 for elements of the c-sized box diagonal and 
    0 for elements in the off box diagonals.

    Parameters: 
        m (int): size of the matrix
        c (int): size of the boxes
    '''
    #Initialize matrix
    M = np.zeros(shape = (m, m))
    #Create vector of indices for upper diagonal
    inds = np.triu_indices(m)
    n_it = len(inds[0])
    #Loop over all index pairs of the upper diagonal
    for ind in range(n_it):
        i = inds[0][ind]
        j = inds[1][ind]
        #Number of resources per class 
        mc = m/c
        if np.floor(i/mc) == np.floor(j/mc):
            #The element is in a class
            M[i, j] = 1
            M[j, i] = 1
    return(M)

def preferences_number(beta, m):
    '''
    Sample the number of preferences from an exponential distribution
    '''
    #Make sure that it is greater than 1, but smaller than the total number
    #of resources
    n = 0
    while n < 1 or n > m:
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

def preferences_matrix(Kc, m, nc, c, s):
    '''
    Construct preference matrix
    Parameters:
        Kc (float): Taxonomic heterogeneity constant
        m (int): Number of resources
        s (int): Number of species
        nc (int): Number of species per resource class
        c (int): Number of classes 
    '''
    #Preallocate matrix of preferences
    c_mat = np.zeros(shape = (s, m))
    #Equal number of species in each class
    Mc = int(np.floor(s/c))
    #Create vector of classes
    class_vec = np.repeat(np.arange(c), Mc)
    for i in range(s):
        c_mat[i,:] = preferences_vector(Kc, m, Mc, 
                                        class_vec[i])
    return c_mat

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
            #bigger value
            d[prod] = 1/(spar*Mc)
        #Note that the values of d match independently of  the class of 
        #product and substrate when kf = 0.
    return d

#Functions for general metabolic sampling


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
    return p_j

def sampling_probability(kc, cumulative_demand, m):
    '''
    Compute the sampling probability vector
    '''

    if not any(cumulative_demand):
        return(1/m*np.ones(m))
    else:
        base = (1 - kc) * 1 / m
        #Make positive 
        cumulative_pos = cumulative_demand + min(cumulative_demand) 
        regime = kc * (cumulative_pos)/sum(cumulative_pos)
        #Transform nans to 0
        regime[np.isnan(regime)] = 0
        #Sum up both terms
        score = base + regime
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

def gen_pref_matrix(Kc, kc, m, nsp, sc_vec, mc_vec, beta = 5):
    '''
    Sample preference matrix
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

def general_metabolic_matrix(s, kf, d, Kf, M, m, nsp):
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
        c_vec = general_concentration_vector(kf, s, M, d, Kf, m, i, nsp)
        D[i,:] = np.random.dirichlet(c_vec)
    return D

def general_concentration_vector(kf, spar, M, d, Kf, m, alpha, nsp):
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

def competition_term(N_a, N_b, c_aj, c_bj):
    return(N_a/N_b*c_aj*c_bj)

def competition_flux(N, C, s, m):
    J = np.zeros(m)
    for j in range(m):
        for a in range(s):
            for b in range(s):
                if b != a:
                    J[j] += competition_term(N[a], N[b], C[a, j], C[b, j])
    return(J)

def competition_flux_vec(N, C, s):
    return(np.ones(s) @ np.diag(N) @ C @ (np.diag(1/N) @ C).T)

def facilitation_flux(l, N, C, D, s):
    return(l * (np.diag(np.ones(s) @ np.diag(N) @ C @ D) @ \
                (np.diag(1/N) @ C).T) @ np.ones(s))
                
def one_preference_vector(m, kc, p, beta, n_r):
    '''
    Sample one preference vector
    '''
    #Preallocate vector of preferences of species alpha
    c_vec = np.zeros(m)
    #Sample number of imported metabolites
    #If all species having the same number of preferences
    if n_r:
        n = round(n_r)
    #If Number of preferences follows an exponential distribution with rate
    #beta
    elif beta:
        #Make sure that it is greater than 1, but smaller than the total number
        #of resources
        n = 0
        while n < 1 or n > m/3:
            n = round(np.random.exponential(beta))
    #Draw n indices with probability p
    ind = np.random.choice(range(m), n, replace = False, p = p)
    #Flip the ind positions in the vector of strategies for species sp
    c_vec[ind] = 1
    return c_vec

def extract_diag_blocks(a, N, M):
    '''
    Extract matrices in a block diagonal matrix into a list of arrays
        N: Number of rows of blocks
        M: Number of columns of block
    '''
    if len(a) == 1:
        return(a)
    else:
        return([a[i*N:(i+1)*N,i*M:(i+1)*M] for i in range(a.shape[0]//N)])

def block_diagonal_bool(m, n):
    '''
    Create a square block diagonal matrix of booleans
        m: Size of the matrix
        n: Size of the blocks
    '''
    

def demands(c_vecs, D_mats, l, m):
    '''
    Calculate effective demands
    '''
    c_mat = np.vstack(c_vecs)
    D_block_mat = block_diag(*D_mats)
    harvest = np.sum(c_mat, axis = 0)
    c_block_mat = block_diag(*c_vecs)
    secretion_vecs = extract_diag_blocks(c_block_mat @ D_block_mat, int(1), m)
    secretion = np.sum(np.vstack(secretion_vecs), axis = 0)
    return(harvest - l*secretion)

def sample_metabolism(m, s, kc, beta, kf, K, M, nsp, l, spar):
    '''
    Compute metabolic matrix C and secretion matrices D_alpha 
    Parameters: 
        m (int): Number of metabolites
        s (int): Number of species
        kc (float): Competitivenes factor
        beta (float): Exponential parameter
        kf (float): Facilitation constant
        K (float): Structure constant
        M (mxm array): Matrix of classes
        nsp (int): Number of individuals in each class
        l (float): Leakage
        spar (float): Sparsity parameter
    '''
    #Initialize demands
    d = np.zeros(m)
    #Compute sampling probability of each metabolite
    p = sampling_probability(kc, d, m)
    #Preallocate list of preference vectors
    c_vecs = list()
    D_mats = list()
    for i in range(s):
        #Sample preference vector of species alpha
        c_vec = one_preference_vector(m, kc, p, beta, n_r = None)
        #Add metabolic vector to the list
        c_vecs.append(c_vec)
        #Sample metabolic matrix of species alpha
        D = general_metabolic_matrix(spar, kf, d, K, M, m, nsp)
        #Add metabolic matrix to the list
        D_mats.append(D)
        #Update vector of demands
        d = demands(c_vecs, D_mats, l, m)
    #Convert to block diagonal matrix
    c_mat = block_diag(*c_vecs)
    D_mat = block_diag(*D_mats)

    return(list([c_mat, D_mat]))

class Community:
    def __init__(self, kf, kc, r, c, n, x, l, D, R, n_coal):
        self.kf = kf #Facilitation factor
        self.kc = kc # Competition factor
        self.r = r #Community richness
        self.c = c #Metabolic preferences matrix
        self.n = n #Vector of species abundances
        self.x = x #Maintenance vector
        self.l = l #Leakage 
        self.D = D #Metabolic matrix
        self.R = R #Concentration of resources before/after assembly/coalescence
        self.s = np.shape(c)[0] #Number of species before assembly
        self.m = np.shape(c)[1] #Number of metabolites
        self.coal = n_coal #Number of coalescence events to which exposed
        self.source = np.repeat(n_coal, r) #Community-level composition
        self.path = np.mean(np.sum(self.c, axis = 1)) #Facilitation factor
        self.indext = 0 #Indices of extinctions

    def assembly(self, kappa):
        #Create dictionary of parameters
        params = {'g':np.ones(self.r).reshape(self.r, 1),
                  's': self.r,
                  'm':self.m,
                  'K':2*kappa.reshape(self.m, 1),
                  't':0.5*np.ones(self.m).reshape(self.m, 1),
                  'coal':self.coal,#This is not a coalescence event
                  'D':self.D,
                  'c':self.c,
                  'x':self.x.reshape(self.r, 1),
                  'l':self.l.reshape(self.r, 1)}
        #Initial conditions
        z0 = list(self.n) + list(self.R) 
        #Create time vector
        tspan = tuple([1, 1e4])
        #Integrate diferential equations
        sol = solve_ivp(lambda t,z: equations(t,z, params),
                        tspan, z0,
                        method = 'BDF', atol = 0.0001 )
        #Record and population abundances at each time point
        N_t = sol.y[0:self.r, :] #sxt matrix
        #Record resource abundances at equilibrium
        R = sol.y[self.r:self.r + self.m, -1]
        #Get indices of extinctions
        ind_ext = np.where(N_t[:, -1] < 1)[0]
        self.indext = ind_ext
        #Record stable state
        self.n = N_t[:,-1]
        #np.delete(N_t, ind_ext, axis = 0)[:,-1]
        self.x = np.delete(self.x, ind_ext) 
        #Get new number of species after dynamics
        self.r = len(self.n) - len(self.indext)
        #Eliminate rows from matrix c
        self.c = np.delete(self.c, ind_ext, axis = 0)
        #Eliminate rows from vector l
        self.l = np.delete(self.l, ind_ext)
        #Eliminate rows from vector of source
        self.source = np.delete(self.source, ind_ext)
        #We count assembly as one coalescence process
        self.coal = 1
        self.R = R
        return(self)

    
    def coalescence(self, comm, resources):
        '''Perform coalescence between 2 communities'''
     
        #Create joint system
        #Remember to take away the extinct species
        ext_system = joint_system(self.c, self.D, self.n, self.l, self.x, 
                                  comm.c, comm.D, comm.n, comm.l, comm.x)
        
        #Create the vector of supplied resources
        kappa = np.zeros(self.m)
        kappa[resources] = 1
        #Create dictionary of parameters
        params = {'g':np.ones(self.r + comm.r).reshape(self.r + comm.r, 1),
                  's':self.r + comm.r,
                  'm':self.m,
                  'K':2*kappa.reshape(self.m, 1),
                  't':0.25*np.ones(self.m).reshape(self.m, 1),
                  'coal':self.coal,#This is a coalescence event
                  'D':ext_system['D'],
                  'c':ext_system['C'],
                  'x':ext_system['x'],
                  'l':ext_system['l']}
        #Initial conditions
        z0 = list(ext_system['N']) + list(2*np.ones(self.m)) 
        #Create time vector
        tspan = tuple([1, 1e4])
        #Integrate diferential equations
        try:
            sol = solve_ivp(lambda t,z: equations(t,z, params),
                            tspan, z0,
                            method = 'BDF', atol = 0.0001 )
        except:
            #The solution diverges
            sol = None  
        return(params, sol)

def using_multiindex(A, name, columns):
    shape = A.shape
    index = pd.MultiIndex.from_product([range(s)for s in shape], 
                                       names=columns)
    df = pd.DataFrame({name: A.flatten()}, index=index).reset_index()
    return df

def jacobian(l, g, n, C, R, D, tau):
    '''
    Calculate the Jacobian numerically for each equilibrium
    Parameters:
        l (float): leakage term 
        g (float): energy-biomass conversion factor
        n (sx1): equilibrium population abundance
        C (sxm): metabolic preferences matrix 
        R (mx1): equilibrium resource abundance
        D (mxm): Metabolic matrix
        tau (mx1): Dilution rate vector
    '''

    s = len(C)
    #Calculate each block 
    dndn = np.zeros(shape = (s, s))
    dndr = (1-l)*np.diag((g*n).T[0]) @ C
    drdn = - np.diag(R.T[0]) @ C.T + l * D.T @ np.diag(R.T[0]) @ C.T 
    drdr = -np.diag((1/tau).T[0]) - np.diag((n.T @ C)[0]) + l * D.T @ \
           np.diag(n.T[0] @ C)
    #Put blocks together
    J = np.bmat([[dndn, dndr], [drdn, drdr]])
    return(J)

def remove_all(A, s):
    '''
    Collect all square matrices a of size m from matrices  separatelly from 
    the bigger matrix A
    '''
    #Get number of matrices to collect
    n_mat = np.shape(A)[0]//s
    #Preallocate tensor to store all  D matrices
    all_a = np.zeros(shape = (n_mat, s, s))
    for i in range(n_mat):
        all_a[i,:,:] = remove_square(i, s, A)
    return all_a


def remove_square(k, s, A):
    '''
    Remove the kth square matrix of size s inside matrix A
    '''
    return(A[s*k:s*(k+1), s*k:s*(k+1)])

def remove_chunk(extint_ind, D, m):
    '''
    Remove a matrix inside a bigger matrix
    '''
    #Create matrix of bools of the size of Dtot
    bool_mat = np.ones(shape = np.shape(D), dtype = bool)
    shape = np.array(np.shape(D))
    #Create crosses of False to be removed from D
    for i in extint_ind:
        #Create chunks of False on the rows
        bool_mat[i*m:(i+1)*m, :] = False
        #Create chunks of False on the columns
        bool_mat[:, i*m:(i+1)*m] = False
    #Subtract the chunks from D
    D = D[bool_mat]
    #Reshape D to be square again
    D = D.reshape(tuple(shape - len(extint_ind)*m))
    return D

def full_c(mix_comm, all_comms):
    '''
    Create the mixed  matrix with elements in the off diagonals. Specifically,
    each row is a copy of the matrix in the box-diagonal
    '''
    #Get vector of distinct communities and how many species each community has 
    values, counts = np.unique(mix_comm.source, return_counts = True)
    #Get cummulative sum of cunts, add 0 at the begginning
    cum_counts = np.hstack([0, np.cumsum(counts)])
    #Initialize matrix
    c_mix = np.zeros(np.shape(mix_comm.c))
    m = mix_comm.m
    #Loop over each communtiy
    for i in range(len(values)):
        #Repeat c_mat of the ith community, as many times as different
        #communities are
        c_chunk = mix_comm.c[cum_counts[i]:cum_counts[i+1], i*m:(i+1)*m]
        c_mix[cum_counts[i]:cum_counts[i+1],:] = np.tile(c_chunk, len(values))
    return c_mix

def remove_stripe(extinct, c, m):
    '''
    Remove a vertical stripe along a matrix
    '''
    #Create a boolean matrix to mask c
    shape = np.shape(c)
    bool_mat = np.ones(shape, dtype = bool)
    #Make a vertical stripe of Falses in the place of extinct
    #communities
    for i in extinct:
        bool_mat[:,(i)*m:(i+1)*m] = False
    cm = c[bool_mat]
    cm = cm.reshape(shape[0], shape[1] - m*len(extinct))
    return cm

def remove_stripe_revisited(C, m):
    '''
    Find empty chunks and recursively remove them
    '''
    zero_stripes = check_block_zeros(C, m)
    Cor = C
    if np.all(zero_stripes == 0):
        return(C)
    else:
        #Get number of consecutive stripes of width m
        n_stripes = len(C.T)//m
        for i in range(n_stripes):
            #Select ith block
            stripe_i = C[:, i*m:m*(i+1)]
            #Check if all elements are 0
            if np.all(stripe_i == 0):
                C = np.delete(C, slice(m*i, m*(i+1)), axis = 1)
        return remove_stripe_revisited(C, m)

def check_block_zeros(C, m):
    '''
    Check if there are stripes of zeros of width m
    '''
    #Get number of consecutive stripes of width m
    n_stripes = len(C.T)//m
    zeros = np.zeros(n_stripes)
    for i in range(n_stripes):
        #Select ith block
        stripe_i = C[:, i*m:m*(i+1)]
        #Check if all elements are 0
        if np.all(stripe_i == 0):
            zeros[i] = 1
    return zeros
         
def abundances_n_preferences(abundances, preferences, m):
    '''
    Calculate total abundance of group of species with i number of preferences
    '''
    n_pref = np.sum(preferences, axis = 1)
    #Preallocate
    N_n_pref = np.zeros(m)
    for i in range(m):
        #Get indices of species with i number of preferences
        ind_i = np.where(n_pref == i+1)[0]
        N_n_pref[i] = sum(abundances[ind_i])
    return N_n_pref


