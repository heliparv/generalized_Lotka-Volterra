import numpy as np
from scipy.linalg import solve as slv
from scipy.stats import bernoulli

"""steady_state_glv: Function finds a single steady state for a generalized Lotka-Volterra model with only
pairwise interactions when model parameters are given as an input
loo_steady_state_glv: Function finds single steady states for all leave-one-out cultures for
a generalized Lotka-Volterra model with only pairwise interactions when model parameters are given

Inputs:
n = number of species
ri = intrinsic rate of growth for each species, a numpy vector
ineractions_matrix = matrix of pairwise interactions in the system
"""

def steady_state_glv(n, ri, interactions_matrix):
    abundance_n_member_community = slv(interactions_matrix, -ri*np.ones(n))
    return abundance_n_member_community

def loo_steady_state_glv(n, ri, interactions_matrix):
    loo_abundances = np.zeros((n,n))
    for i in range(0, n):
        temp_interactions = np.delete(interactions_matrix, i, 0)
        temp_interactions = np.delete(temp_interactions, i, 1)
        abund = slv(temp_interactions, -ri*np.ones((n-1,1)))
        abund = np.insert(abund, i, 0)
        loo_abundances[i] = abund
    
    return loo_abundances