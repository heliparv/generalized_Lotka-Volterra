import numpy as np
import pandas as pd
from parameters import add_sparcity, generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, calculate_carrying_capacities
from Cyclic_games import generalized_rps, even_groups_for_rps, random_groups_for_rps
from simulations_gLVwK import only_viable_gLVwK, gLVwK_with_extinction, test_gLVwK
from graphics import abundances_line_chart, interactions_heatmap

"""Index for performing simulations based on generalized Lotka-Volterra dynamics with carrying
capacity caps and cyclic games in the interactions. Calls appropriate functions for generating
interaction values, starting abundances and for generating trajectories.

To modify simulations, modify input values for functions. Standard deviation input for generating
interactions should be in the desired range, but for generating intrinsic growth rates and pairwise
interactions the standard deviation input should be as a fraction of the given mean. For example if
mean is 100 CFU and desired std is 10 CFU the input should be 0.1

For cyclic games with function generalized_rps choose input function for the way species are grouped
as either even_groups_for_rps or random_groups_for_rps
"""

#Number of species
n = 10
#Maximum simulation time
maxtime = 100
total_carrying_capacity = 20000000

ri = generate_growth_rates(n, 0.5,seed_growth=12, std=0.1)

starting_abundances = generate_starting_abundances(n, seed_abundance=12, mean=100, std=0.1)

pairwise_interactions = generate_interactions(n, 1,seed_interactions=12,seed_sparcity=12, mean=0, std=0.001)

pairwise_interactions = add_sparcity(pairwise_interactions, 0.3,seed_sparcity=12)

pairwise_interactions = generalized_rps(pairwise_interactions, 6, 2, even_groups_for_rps,seed_groups=12,seed_sparcity=12, sparcity=0.2, interaction_std=0.1)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions,seed_selfinter=12, mean=-0.001, std=0.1)

carrying_capacities = calculate_carrying_capacities(ri, pairwise_interactions)

abundances = gLVwK_with_extinction(n, maxtime, pairwise_interactions, ri, carrying_capacities, starting_abundances, total_carrying_capacity)

if type(abundances) == int:
    if abundances == -1:
        print("Nonviable")
    else:
        print("Encountered error")
else:
    print(abundances[-1])
    abundances_line_chart(n, maxtime, abundances)
    #interactions_heatmap(n, pairwise_interactions)
