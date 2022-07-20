import random
import numpy as np
from scipy.stats import bernoulli
from math import comb
import random
from itertools import permutations
from random import shuffle, sample
#from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions

""" 
generalized_rps: Creates cyclic dynamics in the given matrix of pairwise interactions according to
the given specifications. Input groups is the maximum size of cycles, distance is the maximum
distance along the cycle with which species can interact and sparcity describes the number or
fraction of species left outside of the simulation. grouping_function is input of either
even_groups_for_rps or random_groups_for_rps depending on the desired way the species are grouped.

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

def generalized_rps(pairwise_interactions, groups_total, distance, sparcity, grouping_function):
    #There's still something funky about this, but it mostly gives a symmetric interactions matrix, I'll get back to it - Heli
    n = len(pairwise_interactions[0])
    chosen_species = drop_species_from_rps(n, sparcity)
    grouping = grouping_function(chosen_species, groups_total)

    for winner_group in range(groups_total):
        loser_groups = list(range(winner_group+1, winner_group+distance+1))
        for i in range(len(loser_groups)):
            if loser_groups[i] >= groups_total:
                loser_groups[i] = loser_groups[i]-groups_total

        for winner in grouping[winner_group]:
            for group in loser_groups:
                for loser in grouping[group]:
                    pairwise_interactions[winner][loser] = abs(pairwise_interactions[winner][loser])
                    pairwise_interactions[loser][winner] = -pairwise_interactions[winner][loser] 

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
    if n%groups != 0:
        extra_groups = len(grouping)-groups
        extras = []
        for i in range(1, extra_groups+1):
            extras = extras + grouping.pop(-i)
        targets = sample(range(groups), len(extras))
        for i in range(len(extras)):
            grouping[targets[i]].append(extras[i])
    return grouping

def random_groups_for_rps(species, groups):
    grouping = [[] for _ in range(groups)]
    for i in species:
        grouping[sample(range(groups),1)[0]].append(i)
    return grouping