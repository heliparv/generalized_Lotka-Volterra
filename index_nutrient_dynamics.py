import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, generate_positive_vector
from Cyclic_games import generalized_rps, even_groups_for_rps
from simulations_nutrient_dynamics import only_viable_nd, nd_with_extinction
from graphics import abundances_line_chart, interactions_heatmap, abundances_and_nutrient_chart

#Number of species
n = 5
#Maximum simulation time
maxtime = 4

ri = generate_growth_rates(n, 0.2,seed_growth=14, std=0.1)

Ki = generate_positive_vector(n, 12, 20, 0.3)

starting_nutrient = 30

gamma = generate_positive_vector(n, 12, 0.001, 0.2)

nutrient_influx = 5

starting_abundances = generate_starting_abundances(n,seed_abundance=12,mean= 100, std=0.1)

pairwise_interactions = generate_interactions(n, 1,seed_interactions=12,seed_sparcity=12,mean= 0,std= 0.001,sparcity= 0)

pairwise_interactions = generalized_rps(pairwise_interactions, 4, 1, even_groups_for_rps, 12, 12, 0, 0.01)

pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions,seed_selfinter=12, mean=-0.008,std= 0.1)

abundances, nutrient = only_viable_nd(n, maxtime, pairwise_interactions, ri, Ki, starting_abundances, starting_nutrient, gamma, nutrient_influx)
print(abundances)
print(nutrient)


if type(abundances) == int:
    if abundances == -1:
        print("Nonviable")
    else:
        print("Encountered error")
else:
    print(abundances[-1])
    abundances_and_nutrient_chart(n, maxtime, abundances, nutrient)
    #interactions_heatmap(n, pairwise_interactions)
