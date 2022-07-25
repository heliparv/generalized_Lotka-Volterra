### LIBRARIES ###

import numpy as np
from simulations_simple_gLV import stochastic_simple_gLV_with_extinction
from graphics import abundances_line_chart, abundances_line_chart_with_ET

### DESCRIPTION ###

"""
Functions used in the simulations of mutant invasions.

Functions:
r_mutant_calculator: this function calculates the growth rate of the mutant after it is inserted in a leave-one-out 
community at equilibrium. The equilibrium abundances should include the invasion abundance of the mutant.

invasion: returns the abundances of all species (including the mutant) after a single invasion of a mutant (if the 
leave_one_out equilibrium can be obtained). It also returns the growth rates of the species, which are the same as the initial ones, 
except for the mutant (which now obeys the result calculated by the r_mutant_calculator).

ET_calculator: calculates the establishment threshold (N*sigma_i)/((2*N*r_mutant)-mean(sigma)). 
It is used in the multiple single species invasion function.

multiple_different_species_invasions: cycles the species checking if they can invade or not the remaining species 
at equilibrium (with invasion abundances being defined by the given initial abundances), 
according to the invasion function above. A plot is printed for each successful try and the function returns the grid with the 
final abundances after the invasion event (the rows specify the species that attempted the invasion).

multiple_single_species_invasions: given an array of different invasion abundances, the function checks if a single mutant
can invade the leave-one-out community at equilibrium with the different abundances. ETs are calculated and returned according
to the ET_calculator above, together with the final abundances after each invasion attempt.

"""

### INVASION FUNCTIONS ###

def r_mutant_calculator(loo_interactions, loo_equilibrium_abundances, mutant_position):
    r_mutant = sum((-1/loo_equilibrium_abundances[mutant_position])*(loo_interactions[mutant_position]*loo_equilibrium_abundances))
    return r_mutant

def invasion(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma, mutant):
    mutant = mutant - 1
    # Storing starting abundance of the mutant
    N_invasion = starting_abundances[mutant]
    # Deleting starting abundance of the mutant (so we can find leave_one_out equilibrium)
    starting_abundances[mutant] = 0
    # Storing sigma of the mutant to use it later (equating it to zero eliminates the noise of the mutant, which is needed in the leave-one-out equilibrium).
    save_sigma_mutant = sigma[mutant]
    sigma[mutant] = 0
    # Notice "loo" means leave-one-out in what follows.
    loo_equilibrium_abundances = np.zeros(n)
    t = 0
    while t<maxtime:
        loo_equilibrium_abundances = stochastic_simple_gLV_with_extinction(n, maxtime, pairwise_interactions, ri, starting_abundances, sigma)
        if type(loo_equilibrium_abundances) != int:
            loo_equilibrium_abundances = loo_equilibrium_abundances[-1]
            loo_equilibrium_abundances[mutant] = N_invasion
            if all(loo_equilibrium_abundances) > 0:
                break
            else:
                t+=1
        else:
            return loo_equilibrium_abundances, ri
    if t>=maxtime:
        return -3, ri
    # Find the new growth rate of the mutant (based on the loo equilibrium)
    r_mutant = r_mutant_calculator(pairwise_interactions, loo_equilibrium_abundances, mutant)
    ri[mutant] = r_mutant
    # Reassign stored variables
    sigma[mutant] = save_sigma_mutant
    starting_abundances[mutant] = N_invasion
    final_abundances = stochastic_simple_gLV_with_extinction(n, maxtime, pairwise_interactions, ri, starting_abundances, sigma)  
    return final_abundances, ri

def ET_calculator(sigma, ri, mutant_position, starting_abundances):
    denominator = (2*sum(starting_abundances)*ri[mutant_position])-np.mean(sigma)
    ET = (sum(starting_abundances)*sigma[mutant_position])/denominator
    return ET

def multiple_different_species_invasions(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma):
    # Here the invasion_abundances refer to the n different species
    final_abundances_after_invasions = np.zeros((n,n))
    for i in range(n):
        final_abundances, r_mutant = invasion(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma, i)
        if type(final_abundances) == int:
            if final_abundances == -1:
                print("Nonviable")
                print()
            elif final_abundances == -3:
                print("Could not form leave-one-out equilibriumm in maxtime")
                print()
            else:
                print("Encountered error")
                print()
        else:
            if final_abundances[-1][i] >= 0.0005:
                print(f"Mutant {i+1} can invade")
                print()
                if all(final_abundances[-1]) > 0:
                    print(f"Coexistence possible when mutant {i+1} invades")
                    print()
                else:
                    print(f"Coexistence impossible when mutant {i+1} invades")
                    print()
            else:
                print(f"Mutant {i+1} cannot invade")
                print()
            final_abundances_after_invasions[i] = final_abundances[-1]
            abundances_line_chart(n, maxtime, final_abundances)
            #interactions_heatmap(n, pairwise_interactions)
    return final_abundances_after_invasions

def multiple_single_species_invasions(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma, mutant, invasion_abundances):
    mutant = mutant - 1
    # Here the invasion_abundances refer to the same mutant species
    k = np.size(invasion_abundances)
    final_abundances_after_invasions = np.zeros((k,n))
    ETs = np.zeros(k)
    for i in range(k):
        starting_abundances[mutant] = invasion_abundances[i]
        final_abundances, r_mutant = invasion(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma, mutant)
        if type(final_abundances) == int:
            if final_abundances == -1:
                print("Nonviable")
                print()
            elif final_abundances == -3:
                print("Could not form leave-one-out equilibriumm in maxtime")
                print()
            else:
                print("Encountered error")
                print()
        else:
            if final_abundances[-1][mutant] >= 0.0005:
                print(f"Mutant {mutant+1} can invade with invasion abundance {invasion_abundances[i]}")
                print()
                if all(final_abundances[-1]) > 0:
                    print(f"Coexistence possible when mutant {mutant+1} invades with invasion abundance {invasion_abundances[i]}")
                    print()
                else:
                    print(f"Coexistence impossible when mutant {mutant+1} invades with invasion abundance {invasion_abundances[i]}")
                    print()
            else:
                print(f"Mutant {mutant+1} cannot invade with initial abundance {invasion_abundances[i]}")
                print()
            final_abundances_after_invasions[i] = final_abundances[-1]
            ET = ET_calculator(sigma, r_mutant, mutant, starting_abundances)
            ETs[i] = ET
            abundances_line_chart_with_ET(n, maxtime, final_abundances, ET)
    return final_abundances_after_invasions, np.around(ETs, decimals=4)
