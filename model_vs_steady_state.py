import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions
from simulations_simple_gLV import test_simple_gLV, only_viable_simple_gLV, simple_gLV_with_extinction
from graphics import abundances_line_chart, interactions_heatmap
from GLV_steady_state_calculator import steady_state_glv

"""Modelling with generalized Lotka-Volterra is compared to solving steady state abundances
from the set of linear equations. Steady state abundances solved from equations are used
as starting values for the model
"""

#Number of species
n = 3
#Maximum simulation time
maxtime = 100

ri = generate_growth_rates(n, 0.05, 0.1)

pairwise_interactions = generate_interactions(n, 1, 0.01, 0.3)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.2, 0.1)

steady_state_abundances = steady_state_glv(n, ri, pairwise_interactions)

abundances = test_simple_gLV(n,maxtime, pairwise_interactions, ri, steady_state_abundances)

print(steady_state_abundances)
print(abundances[-1])
abundances_line_chart(n, maxtime, abundances)
interactions_heatmap(n, pairwise_interactions)