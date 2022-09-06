import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.
Only pairwise interactions. All functions test for achieving steady state and if there's no change
to any of the abundances during a simulation round the simulation is interrupted and steady state
abundances are filled in for the remainder of the simulation time.

functions:
only_viable_simple_gLV: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation halts and returns "-1". If overflow encountered returns "-2"
in stead of abundances.

simple_gLV_with_extinction: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time. If overflow encountered returns "-2"
in stead of abundances.

test_simple_gLV: No control for negative or zero values, useful for testing how different starting
values affect the calculated abundances. If overflow encountered, prints out values causing overflow
and returns abundances calculated so far.

stochastic_simple_gLV_with_extinction: does the same thing as the simple_gLV_with_extinction above, except
that now it includes fluctuations that are dependent on the sigma vector and dW (note dW is sampled as is, 
because the time step in the simulations is dt=1, otherwise it would have a different value). These fluctuations are
added while updating abundances.

simple_gLV_return_steady_state_abundances: Returns steady state abundances or abundances at given maximum time.
Also returns a variable that informs if steady state was reached.

"""

def only_viable_simple_gLV(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    nonviable = False
    while time < max_increments+1:
        steady = True
        for species in range(0, n):
            try:
                change_per_capita = (ri[species] + sum(abundances[time-1]*interactions[species]))*time_increment
                if abs(change_per_capita) > 0.0001 and steady:
                    steady = False
                new_abundance = np.around((abundances[time-1][species] + abundances[time-1][species]*change_per_capita), decimals=6)
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
    if steady and time < max_increments:
        return add_steady_state_tail(abundances, time, max_increments)
    return abundances

def simple_gLV_with_extinction(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < max_increments+1:
        i = 0
        steady = True
        while i < len(livespecies):
            try:
                change_per_capita = (ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]]))*time_increment
                if abs(change_per_capita)>0.0001  and steady:
                    steady = False
                new_abundance = np.around((abundances[time-1][livespecies[i]] + abundances[time-1][livespecies[i]]*change_per_capita), decimals=6)
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
    if steady and time < max_increments:
        return add_steady_state_tail(abundances, time, max_increments)
    return abundances

def test_simple_gLV(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < max_increments+1:
        steady = True
        for species in range(0, n):
            try:
                change_per_capita = (ri[species] + sum(abundances[time-1]*interactions[species]))*time_increment
            except:
                print(f"Time: {time}")
                print(f"species: {species}")
                print(f"abundances t-1: {abundances[time-1]}")
                print(f"interactions: {interactions[species]}")
                return -2
            if abs(change_per_capita) > 0.0001 and steady:    
                steady = False
            try:
                new_abundance = np.around((abundances[time-1][species] + abundances[time-1][species]*change_per_capita), decimals=6)
                abundances[time][species] = new_abundance
            except:
                print(f"Time: {time}")
                print(f"Species: {species}")
                print(f"Abundance t-1: {abundances[time-1][species]}")
                print(f"Change per capita: {change_per_capita}")
                return -2
        if steady:
            print("täälä")
            break
        else:
            time += 1
    if steady and time != max_increments:
        return add_steady_state_tail(abundances, time, max_increments)
    return abundances

def stochastic_simple_gLV_with_extinction(n, maxtime, time_increment, interactions, ri, starting_abundances, sigma):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < max_increments+1:
        i = 0
        steady = True
        while i < len(livespecies):
            try:
                change_per_capita = (ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]]))*time_increment
                if abs(change_per_capita) > 0.0001 and steady:
                    steady = False
                dW = np.random.random() - np.random.random()
                new_abundance = np.around((abundances[time-1][livespecies[i]] + abundances[time-1][livespecies[i]]*change_per_capita) + np.sqrt(abundances[time-1][livespecies[i]]*sigma[livespecies[i]])*dW, decimals=6)
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
    if steady and time < max_increments:
        return add_steady_state_tail(abundances, time, max_increments)
    return abundances

def leave_one_out_simple_gLV_with_extinction(n, maxtime, time_increment, loo_interactions, ri, loo_starting_abundances):
    abundances = []
    for i in range(0,n):
        abundances.append(simple_gLV_with_extinction(n-1, maxtime, time_increment, loo_interactions[i], ri, loo_starting_abundances[i]))
    return np.array(abundances)

def simple_gLV_return_steady_state_abundances(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((2,n))
    abundances[0] = starting_abundances
    index = 0
    time = 1
    livespecies = list(range(0,n))
    while time < max_increments+1:
        i = 0
        steady = True
        while i < len(livespecies):
            try:
                change_per_capita = (ri[livespecies[i]] + sum(abundances[index]*interactions[livespecies[i]]))*time_increment
                if abs(change_per_capita) > 0.0001 and steady:
                    steady = False
                new_abundance = np.around((abundances[index][livespecies[i]] + abundances[index][livespecies[i]]*change_per_capita), decimals=6)
            except:
                return -2, -2
            if new_abundance <= 0:
                del livespecies[i]
                abundances[0][livespecies[i]] = 0
                abundances[1][livespecies[i]] = 0
            else:
                abundances[abs(index-1)][livespecies[i]] = new_abundance
                i+=1
        if steady:
            break
        else:
            time += 1
            index = abs(index-1)
    return steady, abundances[index]

def add_steady_state_tail(abundances, time, max_increments):
    steady_abundances = np.repeat(np.array([abundances[time]]), max_increments-time, axis=0)
    abundances[time+1:] = steady_abundances
    return abundances
