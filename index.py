import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_pairwise_interactions, generate_starting_abundances, adjust_selfinteractions
from trajectories import only_viable_trajectories, trajectories_with_extinction, test_trajectories
from graphics import abundances_line_chart, interactions_heatmap

"""Index for performing simulations based on generalized Lotka-Volterra dynamics. Calls appropriate
functions for generating interaction values, starting abundances and for generating trajectories.
To modify simulations, modify input values for functions. Standard deviation input for generating
interactions should be in the desired range, but for generating intrinsic growth rates and pairwise
interactions the standard deviation input should be as a fraction of the given mean. For example if
mean is 100 CFU and desired std is 10 CFU the input should be 0.1
"""

#Number of species
n = 10
#Maximum simulation time
maxtime = 9000

ri = generate_growth_rates(n, 0.001, 0.1)

starting_abundances = generate_starting_abundances(n, 100, 0.1)

pairwise_interactions = generate_pairwise_interactions(n, 0, 0.001, 0.3)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.001, 0.1)

abundances = test_trajectories(n, maxtime, pairwise_interactions, ri, starting_abundances)

print(abundances[-1])
abundances_line_chart(n, maxtime, abundances)
#interactions_heatmap(n, pairwise_interactions)