import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_abundances, add_sparcity, adjust_selfinteractions_for_carrying_capacity, calculate_carrying_capacities
from simulations_simple_gLV import simple_gLV_with_extinction
from sim_gLV_RK import gLV_RK
from Cyclic_games import random_win_lose_system
from graphics import abundances_line_chart

"""Index for performing simulations based on generalized Lotka-Volterra dynamics. Only pairwise
interactions. Calls appropriate functions for generating interaction values, starting abundances
and for generating trajectories. To modify simulations, modify input values for functions. Standard
deviation input for generating interactions should be in the desired range, but for generating
intrinsic growth rates and pairwise interactions the standard deviation input should be as a fraction
of the given mean. For example if mean is 100 CFU and desired std is 10 CFU the input should be 0.1
"""

#Number of species
n = 20
#Maximum simulation time
maxtime = 30
time_increment = 0.1

ri = generate_growth_rates(n, mean=0.4, seed_growth=12, std=0.1)

starting_abundances = generate_abundances(n,seed=12,mean= 100, std=0.1)

pairwise_interactions = generate_interactions(n, 1,seed_interactions=12,mean= 0,std= 0.0004)

#pairwise_interactions = add_sparcity(pairwise_interactions, sparcity=0.3, seed_sparcity=12)

pairwise_interactions = random_win_lose_system(n, pairwise_interactions, seed_winlose=12, sparcity=0, seed_sparcity=0, interaction_std=0.1, seed_interaction=33, off_target_interactions=False)

carrying_capacities = generate_abundances(n, seed=45, mean=2000, std=0.3)

pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)

abundances = simple_gLV_with_extinction(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)

if type(abundances) == int:
    print("Encountered error")
else:
    print(abundances[-1])
    abundances_line_chart(n, time_increment, abundances)

abundances = gLV_RK(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)

if type(abundances) == int:
    print("Encountered error")
else:
    print(abundances[-1])
    abundances_line_chart(n, time_increment, abundances)