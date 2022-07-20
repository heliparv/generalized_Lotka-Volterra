import random
import numpy as np
from scipy.stats import bernoulli
from math import comb
import random
from itertools import permutations
from random import shuffle, sample
#from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions

""" 
drop_species_from_rps: Add sparcity to the rps model.Takes input of number of species and desired
sparcity. If given sparcity is lower than 1, it is used as a fraction, otherwise it's assumed to
be number of species dropped.

even_groups_for_rps: Takes input of a list of species to choose from after sparcity function used on
the species group and number of groups, which is the cycle size in the cyclic dynamics. Randomly
divides the species to approximately equal group sizes

random_groups_for_rps: Takes input of a list of species to choose from after sparcity function used on
the species group and number of groups, which is the cycle size in the cyclic dynamics. Randomly assigns
each species to one of the groups.

"""

def rpsls(pairwise_interactions,mean=0, std=0.1): 
 
    n=pairwise_interactions.shape[1]
    choice=random.sample(range(n), 5)
    #Each species is assigned numbers from 0 to 5 accoridng to their index and we work with 0-5 from now on
    perm = list(permutations(range(5), 2))
    #Can use their number in module 5 because we could get 0 even if not same species
    diff= [((b - a) % 5) for a, b in perm]
    for i in range(len(diff)):
        j=choice[perm[i][0]]
        k=choice[perm[i][1]]
        if (diff[i]==1 or diff[i]==2):
            pairwise_interactions[k,j]=abs(pairwise_interactions[k,j])
            pairwise_interactions[j,k]=-(pairwise_interactions[k,j])
        if (diff==3 or diff==4):
            pairwise_interactions[k,j]=-abs(pairwise_interactions[k,j])
            pairwise_interactions[j,k]=-(pairwise_interactions[k,j])
    return pairwise_interactions                                           

def drop_species_from_rps(n, sparcity):
    if sparcity < 1:
        keep = int(n*(1-sparcity))
    else:
        keep = n-sparcity
    return sample(range(n), keep)

def even_groups_for_rps(species, groups):
    n = len(species)
    group_size = int(n/groups)
    shuffle(species)
    grouping = [species[i:i+group_size] for i in range(0, n, group_size)]
    group_dict = {}
    if n%groups != 0:
        extras = grouping.pop(-1)
        for i in extras:
            group_dict[i] = sample(range(groups),1)[0]
    for i in range(groups):
        for j in grouping[i]:
            group_dict[j] = i
    return group_dict

def random_groups_for_rps(species, groups):
    group_dict = {}
    for i in species:
        group_dict[i] = sample(range(groups),1)
    return group_dict

