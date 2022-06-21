import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_pairwise_interactions, generate_starting_abundances
from trajectories import only_viable_trajectories
from graphics import abundances_line_chart

"""Index for performing simulations based on generalized Lotka-Volterra dynamics. Calls appropriate
functions for generating interaction values, starting abundances and for generating trajectories.
To modify simulations, modify input values for functions.
"""

#Number of species
n = 5
#Maximum simulation time
maxtime = 15

pairwise_interactions = generate_pairwise_interactions(n)

ri = generate_growth_rates(n)

starting_abundances = generate_starting_abundances(n)

abundances = only_viable_trajectories(n, maxtime, pairwise_interactions, ri, starting_abundances)

if abundances == -1:
    print("Community not viable")
else:
    abundances_line_chart(n, maxtime, abundances)