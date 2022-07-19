import random
import numpy as np
from scipy.stats import bernoulli
from math import comb
import random
from itertools import permutations
#from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions


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
         
