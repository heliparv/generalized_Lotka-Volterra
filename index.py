from __future__ import absolute_import
import numpy as np
import pandas as pd
from scipy.stats import bernoulli
from math import comb

""" A rudimentary version of the generalized Lotka-Volterra model.
Currently the model tends to result in negative values of abundances at certain timepoints
when using default values, so those have to be tweaked.

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
by calling function add_sparcity. Selfinteractions are then set to negative values based on the
originally drawn values.

simulate: Calls previous functions to set up the simulation. This data is then used to simulate
species abundance trajectories with generalized Lotka-Volterra dynamics for all of the species
for as many timepoints as desired.

"""

#TODO separate birth and death rate in model
#TODO take lag times into account
#TODO add function for drawing trajectories

def generate_growth_rates(n, std):
    ri = np.random.normal(loc=(1/n), scale=(1/n)*std, size=n)
    for i in range(0, len(ri)):
        if ri[i] < 0:
            ri[i] = ri[i]*-1
    return ri

def generate_starting_abundances(n, mean, std):
    abundances = np.random.normal(loc=mean, scale=std, size=n).astype(int)
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

def generate_pairwise_interactions(n, mean, std, sparcity):
    interactions_matrix = np.random.normal(loc=mean, scale=std, size=(n,n))

    interactions_matrix = add_sparcity(interactions_matrix, sparcity)

    for i in range(0,n):
        if interactions_matrix[i][i] >= 0:
            interactions_matrix[i][i] = -1*interactions_matrix[i][i]

    return interactions_matrix

def generate_trajectories(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    unviable = False
    while time < maxtime+1:
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            change = abundances[species]*change_per_capita
            abundances[time][species] = sum(abundances[species]+change)
            if abundances[time][species] <= 0:
                unviable = True
                break
        if unviable:
            return(-1)
        time += 1
    return abundances

def simulate(n, maxtime, interaction_mean=0.0, interaction_std=0.01, sparcity=0.1, growth_std=0.1, abundances_mean=10000, abundances_std=0):
    interactions = generate_pairwise_interactions(n, interaction_mean, interaction_std, sparcity)
    ri = generate_growth_rates(n, growth_std)
    starting_abundances = generate_starting_abundances(n, abundances_mean, abundances_std)
    return generate_trajectories(n, maxtime, interactions, ri, starting_abundances)
