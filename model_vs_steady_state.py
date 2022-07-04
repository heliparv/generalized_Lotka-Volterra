import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_pairwise_interactions, generate_starting_abundances, adjust_selfinteractions
from trajectories import only_viable_trajectories, trajectories_with_extinction, test_trajectories
from graphics import abundances_line_chart, interactions_heatmap
from GLV_steady_state_calculator import steady_state_glv

"""Modelling with generalized Lotka-Volterra is compared to solving steady state abundances
from the set of linear equations. Steady state abundances solved from equations are used
as starting values for the model
"""

#Number of species
n = 20
#Maximum simulation time
maxtime = 7000

ri = generate_growth_rates(n, 0.001, 0.1)

pairwise_interactions = generate_pairwise_interactions(n, 0, 0.00001, 0.3)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.00002, 0.1)

steady_state_abundances = steady_state_glv(n, ri, pairwise_interactions)

#abundances = test_trajectories(n, maxtime, pairwise_interactions, ri, steady_state_abundances)

print(steady_state_abundances)
#print(abundances[-1])
#abundances_line_chart(n, maxtime, abundances)
#interactions_heatmap(n, pairwise_interactions)