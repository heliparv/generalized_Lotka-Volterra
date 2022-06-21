import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_pairwise_interactions, generate_starting_abundances
from trajectories import only_viable_trajectories

"""Index for performing simulations based on generalized Lotka-Volterra dynamics. Calls appropriate
functions for generating interaction values, starting abundances and for generating trajectories.
To modify simulations, modify starting values for functions.
"""

#Number of species
n = 5

pairwise_interactions = generate_pairwise_interactions(n)

ri = generate_growth_rates(n)

starting_abundances = generate_starting_abundances(n)

trajectories = only_viable_trajectories(n, 15, pairwise_interactions, ri, starting_abundances)

print(trajectories)