import numpy as np
from scipy.stats import bernoulli
from math import comb
import random

""" Collection of functions used in generating initial values
for the generalized Lotka-Volterra model.

functions:
generate_growth_rates: generates the intrinsic growth rates for species with the option of variation
between species. Values are drawn from a standard distribution where mean and standard deviation can
be chosen freely.

generate_starting_abundances: Generates the abundances of each species at start of simulation
by drawing from a normal distribution where the mean and standard deviation can be set.

add_sparcity: takes a 1- or 2-dimensional numpy array and a sparcity term. Uses sparcity term
to draw a boolean index from a bernouilli distribution, sets chosen indexes to zero on the
given array.

generate_interactions: Generates a matrix of interactions where rows are the affected species and
columns the species effecting it. The order defines the interaction's order, 1 means pairwise
interactions, 2 means tertiary and so on. Interactions are drawn from a normal distribution.
Sparcity is introduced to the interactions by calling function add_sparcity.

adjust_selfinteractions: Takes generated pairwise interactions and adjusts the diagonal terms,
the pairwise interactions, and replaces them with values drawn from a normal distribution with
the given mean and standard deviation.

generate_loo_interaction_matrices: Takes a pairwise interactions matrix and return a list
of interactions matrices for leave-one-out simulations

generate_abundance_call: When modelling higher order interactions the calculated effect of species on
the effected species depends on the abundances of all the effector species. The call vector is created
for a convenient way to calculate the product of the appropriate abundances at each step so that
this vector can be used in the way a vector of single species abundances is used with pairwise interactions.
Doesn't allow for a species to be represented twice in an interaction as an effector, but the model doesn't
control for interactions to be zero if the affected species is one of the effectors

generate_abundance_call_with_repeats: Same as previous function, but allows for the same species to be
repeated in the call for an arbitrary number of times.

calculate_carrying_capacities: Calculates the maximum carrying capacities for bacteria in the gLV with K
model, relates intrinsic growth rate to self-interaction

generate_positive_vector: Draws a vector of numbers from a normal distribution with length, mean,
standard deviation, and seed specified in the call. Returns absolute values of drawn numbers rounded to
4 decimals. Can be used to, for example, draw Ki and gamma for nutrient dynamics

"""

def generate_growth_rates(n, mean, seed_growth,std = 0.1):
    np.random.seed(seed_growth)
    ri = np.random.normal(loc=mean, scale=abs(mean*std), size=n)
    for i in range(0, len(ri)):
        if ri[i] < 0:
            ri[i] = ri[i]*-1
    return ri

def generate_starting_abundances(n,seed_abundance, mean=100, std=0):
    np.random.seed(seed_abundance)
    abundances = np.random.normal(loc=mean, scale=abs(mean*std), size=n).astype(int)
    for i in range(0, len(abundances)):
        if abundances[i] < 0:
            abundances[i] = abundances[i]*-1
    return abundances

def add_sparcity(array, sparcity, seed_sparcity):
    np.random.seed(seed_sparcity)
    draw = bernoulli(sparcity)
    if len(np.shape(array)) == 1:
        sparcity_array = np.array(draw.rvs(np.size(array)), dtype=bool)
    else:
        i, j = np.shape(array)
        sparcity_array = np.array(np.reshape(draw.rvs(i*j), (i, j)), dtype=bool)
    array[sparcity_array] = 0
    return array

def add_pairwise_sparcity(n, pairwise_interactions, sparcity):
    draw = bernoulli(sparcity)
    spar_list = list(map(bool, draw.rvs(int(((n*n)-n)/2))))
    for i in range(n):
        for j in range(i+1, n):
            sparse = spar_list.pop()
            if sparse:
                pairwise_interactions[i][j] = 0
                pairwise_interactions[j][i] = 0
    return pairwise_interactions

def generate_interactions(n, order,seed_interactions, mean=0, std=0.1):
    np.random.seed(seed_interactions)
    columns = comb(n, order)
    interactions = np.around((np.random.normal(loc=mean, scale=std, size=(n,columns))), decimals=8)
    return interactions

def adjust_selfinteractions(n, interactions, seed_selfinter, mean=-0.1, std=0.1):
    np.random.seed(seed_selfinter)
    selfinteractions = np.around((np.random.normal(loc=mean, scale=abs(mean*std), size=n)), decimals=4)
    for i in range(0,n):
        if selfinteractions[i] >= 0:
            interactions[i][i] = -1*selfinteractions[i]
        else:
            interactions[i][i] = selfinteractions[i]
    return interactions

def generate_loo_interaction_matrices(interactions):
    n = len(interactions[0])
    loo_interactions = []
    for i in range(0, n):
        temp_interactions = np.delete(interactions, i, 0)
        temp_interactions = np.delete(temp_interactions, i, 1)
        loo_interactions.append(temp_interactions)
    return loo_interactions

def generate_abundance_call(species_list, available_species, choose):
    if choose == 0:
        return [species_list]
    else:
        table = []
        for i in range(0, len(available_species)-choose+1):
            table = table + generate_abundance_call(species_list+[available_species[i]], available_species[i+1:], choose-1)
        return table

def generate_abundance_call_with_repeats(species_list, available_species, choose):
    if choose == 0:
        return [species_list]
    else:
        table = []
        for i in range(0, len(available_species)-choose+1):
            table = table + generate_abundance_call_with_repeats(species_list+[available_species[i]], available_species[i:], choose-1)
        return table

def calculate_carrying_capacities(ri, interactions):
    carrying_capacities = []
    for i in range(0, len(ri)):
        carrying_capacities.append(np.around((-ri[i]/interactions[i][i]),decimals=4))
    return np.array(carrying_capacities)

def generate_sigma(n, seed_sigma, mean=0.1, std=1):
    np.random.seed(seed_sigma)
    sigma = np.abs(np.random.normal(loc=mean, scale=abs(mean*std), size=n))
    return np.around(sigma, decimals=4)

def generate_positive_vector(n, seed, mean, std):
    np.random.seed(seed)
    return abs(np.around(np.random.normal(loc=mean, scale=mean*std, size=n), decimals=4))
