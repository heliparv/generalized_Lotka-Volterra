import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, generate_abundance_call
from simulations_simple_gLV_HOI import only_viable_simple_gLV_HOI, simple_gLV_with_extinction_HOI, test_simple_gLV_HOI
from graphics import abundances_line_chart, interactions_heatmap

"""Index for performing simulations based on generalized Lotka-Volterra dynamics. Pairwise and 
higher-order interactions. Here modelled up to quaternary interactions. Calls appropriate functions
for generating interaction values, starting abundances and for generating trajectories. To modify
simulations, modify input values for functions. Standard deviation input for generating interactions
should be in the desired range, but for generating intrinsic growth rates and pairwise interactions
the standard deviation input should be as a fraction of the given mean. For example if mean is
100 CFU and desired std is 10 CFU the input should be 0.1
"""

#Number of species
n = 4
#Maximum simulation time
maxtime = 20

ri = generate_growth_rates(n, 0.3,seed_growth=12,std= 0.1)

starting_abundances = generate_starting_abundances(n,seed_abundance=12,mean= 100, std=0.1)

pairwise_interactions = generate_interactions(n, 1, 0, 0.00004, 0.3)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.008, 0.1)

tertiary_interactions = generate_interactions(n, 2, -0.000004, 0.1, 0.5)

#quaternary_interactions = generate_interactions(n, 3, 0, 0.0000004, 0.5)

#HOInteractions = np.concatenate((tertiary_interactions, quaternary_interactions), axis=1)

#for total competition dynamics, all negative interactions:
#pairwise_interactions = -1*abs(pairwise_interactions)
#HOInteractions = -1*abs(HOInteractions)

call_vector = generate_abundance_call([], list(range(n)), 2)
#+ generate_abundance_call([], list(range(n)), 3)

abundances = test_simple_gLV_HOI(n, maxtime, pairwise_interactions, ri, starting_abundances, tertiary_interactions, call_vector)

if type(abundances) == int:
    if abundances == -1:
        print("Nonviable")
    else:
        print("Encountered error")
else:
    print(abundances[-1])
    abundances_line_chart(n, maxtime, abundances)
    #interactions_heatmap(n, pairwise_interactions)
