from __future__ import absolute_import
import numpy as np
import pandas as pd
from scipy.stats import bernoulli

""" A rudimentary version of the generalized Lotka-Volterra model.
Currently the model tends to result in negative values of abundances at certain timepoints
when using default values, so those have to be tweaked.

functions:
generate_growth_rates: generates the intrinsic growth rates for species with the option of variation
between species. Values are drawn from a standard distribution with mean 1/n where n is number of
species and standard deviation can be chosen.

generate_pairwise_interactions: Generates a matrix of pairwise interactions where rows are the
affected species and columns the species affecting it. Interactions are drawn from a normal distribution
with the mean 0.0 and where standard deviation can be chosen. Sparcity to interactions map is
introduced by setting some interactions to zero. Boolean indexing for this is created by drawing
from a bernoulli distribution, where probability of setting a value to zero can be chosen by giving
value for "sparcity" in function call. Selfinteractions are then set to -1

generate_starting_abundances: Generates the abundances of each species at start of simulation
by drawing from a normal distribution where the mean and standard deviation can be set.

simulate: Calls previous functions to set up the simulation. This data is then used to simulate
species abundance trajectories with generalized Lotka-Volterra dynamics for all of the species
for as many timepoints as desired.

"""

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
    abundances = np.random.normal(loc=mean, scale=std, size=n).astype(int)
    for i in range(0, len(abundances)):
        if abundances[i] < 0:
            abundances[i] = abundances[i]*-1
    return abundances

def simulate(n, maxtime, interaction_mean=0.0, interaction_std=0.1, sparcity=0.1, growth_std=0.1, abundances_mean=10000, abundances_std=0):
    interactions = generate_pairwise_interactions(n, interaction_mean, interaction_std, sparcity)
    ri = generate_growth_rates(n, growth_std)
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = generate_starting_abundances(n, abundances_mean, abundances_std)
    
    time = 1
    while time < maxtime+1:
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            change = abundances[species]*change_per_capita
            abundances[time][species] = sum(abundances[species]+change)
        time += 1
    
    return abundances

print(simulate(3, 5))