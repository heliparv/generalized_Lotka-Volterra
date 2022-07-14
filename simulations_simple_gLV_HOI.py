import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.
Simulations include both pairwise and higher-order interactions. All functions test for achieving steady
state and if there's no change to any of the abundances during a simulation round the simulation is
interrupted and steady state abundances are filled in for the remainder of the simulation time.

functions:
only_viable_simple_gLV_HOI: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation halts and returns "-1". If overflow encountered returns "-2"
in stead of abundances.

simple_gLV_with_extinction_HOI: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time. If overflow encountered returns "-2"
in stead of abundances.

test_simple_gLV_HOI: No control for negative or zero values, useful for testing how different starting
values affect the calculated abundances. When simulation encounters error it prints out all values
that were involved in the calculation and returns abundances so abundances before error can be plotted.

calculate_abundance_product_vector: calculates the product of abundances according to the call vector
created for modelling higher order interactions in the model

"""

def only_viable_simple_gLV_HOI(n, maxtime, pw_interactions, ri, starting_abundances, HOInteractions, call_vector):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    nonviable = False
    while time < maxtime+1:
        steady = True
        abund_vector = calculate_abundance_product_vector(call_vector, abundances[time-1])
        if type(abund_vector) == int:
            return -2
        for species in range(0, n):
            try:
                change_per_capita_pairwise = sum(abundances[time-1]*pw_interactions[species])
                change_per_capita_HOI = sum(abund_vector*HOInteractions[species])
                change_per_capita = ri[species] + change_per_capita_pairwise + change_per_capita_HOI
                if change_per_capita != 0 and steady:
                    steady = False
                new_abundance = abundances[time-1][species] + int(abundances[time-1][species]*change_per_capita)
            except:
                return -2
            if new_abundance <= 0:
                nonviable = True
                break
            else:
                abundances[time][species] = new_abundance
        if nonviable:
            return -1
        elif steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def simple_gLV_with_extinction_HOI(n, maxtime, pw_interactions, ri, starting_abundances, HOInteractions, call_vector):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        steady = True
        abund_vector = calculate_abundance_product_vector(call_vector, abundances[time-1])
        if type(abund_vector) == int:
            return -2
        while i < len(livespecies):
            try:
                change_per_capita_pairwise = sum(abundances[time-1]*pw_interactions[livespecies[i]])
                change_per_capita_HOI = sum(abund_vector*HOInteractions[livespecies[i]])
                change_per_capita = ri[livespecies[i]] + change_per_capita_pairwise + change_per_capita_HOI
                if change_per_capita != 0 and steady:
                    steady = False
                new_abundance = abundances[time-1][livespecies[i]] + int(abundances[time-1][livespecies[i]]*change_per_capita)
            except:
                return -2
            if new_abundance <= 0:
                del livespecies[i]
            else:
                abundances[time][livespecies[i]] = new_abundance
                i+=1
        if steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def test_simple_gLV_HOI(n, maxtime, pw_interactions, ri, starting_abundances, HOInteractions, call_vector):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < maxtime+1:
        steady = True
        abund_vector = calculate_abundance_product_vector(call_vector, abundances[time-1])
        if type(abund_vector) == int:
            print("Calculating abundance vector")
            print(f"time: {time}")
            print(f"Abundances: {abundances[time-1]}")
            return abundances
        for species in range(0, n):
            try:
                change_per_capita_pairwise = sum(abundances[time-1]*pw_interactions[species])
            except:
                print("Calculating pairwise change per capita")
                print(f"Time: {time}")
                print(f"species: {species}")
                print(f"abundances t-1: {abundances[time-1]}")
                print(f"interactions: {pw_interactions[species]}")
                return abundances
            #try:
            change_per_capita_HOI = sum(np.around((abund_vector*HOInteractions[species]), decimals=4))
            #except:
             #   print("Calculating change per capita HOI")
              #  print(f"Time: {time}")
               # print(f"species: {species}")
                #print(f"abundance products: {abund_vector}")
               # print(f"interactions: {HOInteractions[species]}")
                #return abundances
            try:
                change_per_capita = ri[species] + change_per_capita_pairwise + change_per_capita_HOI
            except:
                print("Calculating total change per capita")
                print(f"Time: {time}")
                print(f"ri: {ri[species]}")
                print(f"pairwise change per capita: {change_per_capita_pairwise}")
                print(f"HOI change per capita: {change_per_capita_pairwise}")
                return abundances
            if change_per_capita != 0 and steady:    
                steady = False
            try:
                new_abundance = abundances[time-1][species] + int(abundances[time-1][species]*change_per_capita)
                abundances[time][species] = new_abundance
            except:
                print("Calculating new abundance")
                print(f"Time: {time}")
                print(f"Species: {species}")
                print(f"Abundance t-1: {abundances[time-1][species]}")
                print(f"Change per capita: {change_per_capita}")
                return abundances
        if steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def calculate_abundance_product_vector(call, abundances):
    product = np.ones(len(call))
    for i in range(0, len(call)):
        for j in call[i]:
            try:
                product[i] = product[i]*abundances[j]
            except:
                return -1
    return product