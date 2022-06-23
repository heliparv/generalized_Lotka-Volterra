import numpy as np
from scipy.stats import bernoulli
from math import comb

""" Collection of functions used in generating initial values
for the generalized Lotka-Volterra model.

functions:
generate_growth_rates: generates the intrinsic growth rates for species with the option of variation
between species. Values are drawn from a standard distribution with mean 1/n where n is number of
species and standard deviation can be chosen.

generate_starting_abundances: Generates the abundances of each species at start of simulation
by drawing from a normal distribution where the mean and standard deviation can be set.

add_sparcity: takes a 1- or 2-dimensional numpy array and a sparcity term. Uses sparcity term
to draw a boolean index from a bernouilli distribution, sets chosen indexes to zero on the
given array.

generate_pairwise_interactions: Generates a matrix of pairwise interactions where rows are the
affected species and columns the species affecting it. Interactions are drawn from a normal distribution
with the mean 0.0 and where standard deviation can be chosen. Sparcity is introduced to the interactions
by calling function add_sparcity.

adjust_selfinteractions: Takes generated pairwise interactions and adjusts the diagonal terms,
the pairwise interactions, and replaces them with values drawn from a normal distribution with
the given mean and standard deviation.

"""

def generate_growth_rates(n, mean, std = 0.1):
    ri = np.random.normal(loc=mean, scale=mean*std, size=n)
    for i in range(0, len(ri)):
        if ri[i] < 0:
            ri[i] = ri[i]*-1
    return ri

def generate_starting_abundances(n, mean=100, std=0):
    abundances = np.random.normal(loc=mean, scale=mean*std, size=n).astype(int)
    for i in range(0, len(abundances)):
        if abundances[i] < 0:
            abundances[i] = abundances[i]*-1
    return abundances

def add_sparcity(array, sparcity):
    draw = bernoulli(sparcity)
    if len(np.shape(array)) == 1:
        sparcity_array = np.array(draw.rvs(np.size(array)), dtype=bool)
    else:
        i, j = np.shape(array)
        sparcity_array = np.array(np.reshape(draw.rvs(i*j), (i, j)), dtype=bool)
        for i in range(0,i):
            sparcity_array[i][i] = 0
    array[sparcity_array] = 0
    return array

def generate_pairwise_interactions(n, mean=0, std=0.01, sparcity=0.1):
    interactions_matrix = np.random.normal(loc=mean, scale=std, size=(n,n))

    interactions_matrix = add_sparcity(interactions_matrix, sparcity)

    return interactions_matrix

def adjust_selfinteractions(n, interactions, mean=-0.1, std=0.1):
    selfinteractions = np.random.normal(loc=mean, scale=abs(mean*std), size=n)

    for i in range(0,n):
        if selfinteractions[i] >= 0:
            interactions[i][i] = -1*selfinteractions[i]
        else:
            interactions[i][i] = selfinteractions[i]
    
    return interactions
