import numpy as np
from scipy.stats import bernoulli
from math import comb
import random
from itertools import permutations
from random import shuffle, sample
#from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions

""" 
generalized_rps: Creates cyclic games in the given matrix of pairwise interactions according to
the given specifications. Input groups is the maximum size of cycles of interaction chains,
distance is the maximum distance along the cycle with which species can interact, sparcity
describes the number or fraction of species left outside of the simulation. grouping_function is
input of either even_groups_for_rps or random_groups_for_rps depending on the desired way the species
are grouped. Between two species the interactions are chosen so that the effect on winner is the absolute
value of the interaction in the original interactions matrix and for the loser the interction is either
this value in negative or if interaction_std != 0 the interaction strenght is drawn from a normal
distribution centered on the winner's interaction. If within_group_interactions is False, no interactions
outside of the rps dynamics exist, if True then original interaction values for interactions outside
of the system are kept.

drop_species_from_rps: Add sparcity to the rps model.Takes input of number of species and desired
sparcity. If given sparcity is lower than 1, it is used as a fraction, otherwise it's assumed to
be number of species dropped.

even_groups_for_rps: Takes input of a list of species to choose from after sparcity function used on
the species group and number of groups, which is the cycle size in the cyclic dynamics. Randomly
divides the species to approximately equal group sizes

random_groups_for_rps: Takes input of a list of species to choose from after sparcity function used on
the species group and number of groups, which is the cycle size in the cyclic dynamics. Randomly assigns
each species to one of the groups.

random_win_lose_system: Function for defining predation interactions in a more random fashion. Number of
win-lose interactions is the reverse fraction of given sparcity value, win-lose direction is defined
pseudorandomly. If sparcity in the interaction map is desired so that some species do not interact with
eachother, sparcity should be added prior to using this function.

"""                                       

def generalized_rps(pairwise_interactions, groups_total, distance, grouping_function, seed_groups,seed_sparcity, sparcity, seed_interactions, interaction_std, within_group_interactions):
    n = len(pairwise_interactions[0])
    chosen_species = choose_species_rps(n, sparcity,seed_sparcity)
    grouping = grouping_function(chosen_species, groups_total, seed_groups)

    np.random.seed(seed_interactions)
    new_interactions = np.empty_like(pairwise_interactions)
    print(new_interactions)
    if within_group_interactions:
        new_interactions[:] = pairwise_interactions
        print(new_interactions)

    for winner_group in range(groups_total):
        loser_groups = list(range(winner_group+1, winner_group+distance+1))
        for i in range(len(loser_groups)):
            if loser_groups[i] >= groups_total:
                loser_groups[i] = loser_groups[i]-groups_total

        for winner in grouping[winner_group]:
            for group in loser_groups:
                for loser in grouping[group]:
                    interaction = abs(pairwise_interactions[winner][loser])
                    new_interactions[winner][loser] = interaction
                    new_interactions[loser][winner] = -abs(np.random.normal(loc=interaction, scale=interaction*interaction_std))
    return new_interactions

def choose_species_rps(n, sparcity, seed_sparcity):
    random.seed(seed_sparcity)
    if sparcity < 1:
        keep = int(n*(1-sparcity))
    else:
        keep = n-sparcity
    return sample(range(n), keep)

def even_groups_for_rps(species, groups, seed_groups):
    random.seed(seed_groups)
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

def random_groups_for_rps(species, groups, seed_groups):
    random.seed(seed_groups)
    grouping = [[] for _ in range(groups)]
    for i in species:
        grouping[sample(range(groups),1)[0]].append(i)
    return grouping

def random_win_lose_system(n, pairwise_interactions, seed_winlose, sparcity, seed_sparcity, interaction_std, seed_interaction):
    np.random.seed(seed_winlose)
    draw = bernoulli(0.5)
    winlose = list(map(bool, draw.rvs(int(((n*n)-n)/2))))
    np.random.seed(seed_sparcity)
    draw = bernoulli(sparcity)
    skip_list = list(map(bool, draw.rvs(int(((n*n)-n)/2))))
    np.random.seed(seed_interaction)
    for i in range(n):
        for j in range(i+1, n):
            skip = skip_list.pop()
            if skip:
                print(i, j)
                continue
            win = winlose.pop()
            interaction_ij = abs(pairwise_interactions[i][j])
            interaction_ji = np.random.normal(loc=interaction_ij, scale=interaction_ij*interaction_std)
            if win:
                pairwise_interactions[i][j] = interaction_ij
                pairwise_interactions[j][i] = -interaction_ji
            else:
                pairwise_interactions[i][j] = -interaction_ij
                pairwise_interactions[j][i] = interaction_ji
    return pairwise_interactions