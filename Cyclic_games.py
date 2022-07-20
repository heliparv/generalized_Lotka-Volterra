import random
import numpy as np
from scipy.stats import bernoulli
from math import comb
import random
from itertools import permutations
#from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions

""" 
even_groups_for_rps: Takes input of number of species and number of groups, which is the cycle size
in the cyclic dynamics and randomly divides the species to approximately equal group sizes

random_groups_for_rps: Takes input of number of species and number of groups, which is the cycle size
in the cyclic dynamics. Randomly assigns each species to one of the groups.

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
         
def even_groups_for_rps(n, groups):
    group_size, extras = divmod(n, groups)
    choice = list(range(n))
    shuffle(choice)
    grouping = list(choice[i*group_size+min(i, extras):(i+1)*group_size+min(i+1, extras)] for i in range(groups))
    group_dict = {}
    for i in range(groups):
        for j in grouping[i]:
            group_dict[j] = i
    return grouping

def random_groups_for_rps(n, groups):
    group_dict = {}
    for i in range(n):
        group_dict[i] = sample(range(groups),1)
    return group_dict
