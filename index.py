from __future__ import absolute_import
import numpy as np
from scipy.stats import bernoulli
from copy import deepcopy

def generate_growth_rates(n, std):
    ri = np.random.normal(loc=(1/n), scale=(1/n)*std, size=n)
    for i in range(0, len(ri)):
        if ri[i] < 0:
            ri[i] = ri[i]*-1
    return ri

def generate_pairwise_interactions(n, mean, std, sparcity):
    interactions_matrix = np.random.normal(loc=mean, scale=std, size=(n,n))

    draw = bernoulli(sparcity)
    sparcity_matrix = np.array(np.reshape(draw.rvs(n*n), (n, n)), dtype=bool)
    interactions_matrix[sparcity_matrix] = 0

    for i in range(0,n):
        interactions_matrix[i][i] = -1
    
    return interactions_matrix

def generate_starting_abundances(n, mean, std):
    abundances = np.random.normal(loc=mean, scale=std, size=n)
    for i in range(0, len(abundances)):
        if abundances[i] < 0:
            abundances[i] = abundances[i]*-1

def simulate(n, maxtime, interaction_mean=0.0, interaction_std=0.1, sparcity=0.1, growth_std=0.1, abundances_mean=10, abundances_std=0):
    interactions = generate_pairwise_interactions(n, interaction_mean, interaction_std, sparcity)
    ri = generate_growth_rates(n, growth_std)
    abundances = generate_starting_abundances(abundances_mean, abundances_std)
    time = 0
    new_abundances = np.zeros(n)

    while time < maxtime+1:
        for species in range(0, n):
            change_function = ri[species] #TODO + spesies interactions (A*N in gLV)
            change = abundances[species]*change_function
            new_abundances[species] = abundances[species]+change
        abundances = deepcopy(new_abundances)
        time += 1
    
