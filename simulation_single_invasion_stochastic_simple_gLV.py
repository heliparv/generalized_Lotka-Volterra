### LIBRARIES ###
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, generate_sigma
from graphics import abundances_line_chart
from invasion import invasion

### DESCRIPTION ###

"""
Index for performing simulation of a single mutant invasion event based on stochastic version of generalized Lotka-Volterra dynamics. Only pairwise
interactions. Calls appropriate functions for generating interaction values, starting abundances and for generating trajectories. 
To modify simulations, modify input values for functions. 
Standard deviation input for generating interactions should be in the desired range, but for generating
intrinsic growth rates, sigma values, and pairwise interactions the standard deviation input should be as a fraction
of the given mean. For example if mean is 100 CFU and desired std is 10 CFU the input should be 0.1
"""

### INITIAL PARAMETERS ###

#Number of species
n = 10
#Maximum simulation time
maxtime = 400
#Choose mutant species (usually we assume species 1 as the mutant)
mutant = 1

ri = generate_growth_rates(n, 0.4,seed_growth=12,std=0.1)
print("The growth rates: ", ri)
print()

starting_abundances = generate_starting_abundances(n,seed_abundance=12, mean=200, std=0.1)
print("The initial abundances: ", starting_abundances)
print()

pairwise_interactions = generate_interactions(n, 1, 0, 0.0004, 0.3)


pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.0008, 0.1)
print("The pairwise interaction grid: \n", pairwise_interactions)
print()

#for total competition dynamics, all negative interactions:
#pairwise_interactions = -1*abs(pairwise_interactions)

sigma = generate_sigma(n,seed_sigma=12, mean=0.1, std=0.1)
print("The sigma terms of the noise are: \n", sigma)
print()

### INVASION FUNCTION ###

final_abundances, ri = invasion(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma, mutant)

### RESULTS ###

if type(final_abundances) == int:
    if final_abundances == -1:
        print("Nonviable")
    elif final_abundances == -3:
        print("Could not form leave-one-out equilibrium.")
    else:
        print("Encountered error")
else:
    print(f"Final abundances after invasion by mutant {mutant} \n", final_abundances[-1])
    print()
    if final_abundances[-1][mutant-1] >= 0.0005:
        print(f"Mutant {mutant} can invade")
        print()
        if all(final_abundances[-1]) > 0:
            print(f"Coexistence possible when mutant {mutant} invades")
            print()
        else:
            print(f"Coexistence impossible when mutant {mutant} invades")
            print()
    else:
        print(f"Mutant {mutant} cannot invade")
        print()
    abundances_line_chart(n, maxtime, final_abundances)
    #interactions_heatmap(n, pairwise_interactions)
