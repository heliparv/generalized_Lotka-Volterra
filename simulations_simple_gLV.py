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

"""

def only_viable_simple_gLV(n, maxtime, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    nonviable = False
    while time < maxtime+1:
        steady = True
        for species in range(0, n):
            try:
                change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
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

def simple_gLV_with_extinction(n, maxtime, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        steady = True
        while i < len(livespecies):
            try:
                change_per_capita = ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]])
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

def test_simple_gLV(n, maxtime, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < maxtime+1:
        steady = True
        for species in range(0, n):
            try:
                change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            except:
                print(f"Time: {time}")
                print(f"species: {species}")
                print(f"abundances t-1: {abundances[time-1]}")
                print(f"interactions: {interactions[species]}")
                return abundances
            if change_per_capita != 0 and steady:    
                steady = False
            try:
                new_abundance = abundances[time-1][species] + int(abundances[time-1][species]*change_per_capita)
                abundances[time][species] = new_abundance
            except:
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

def stochastic_simple_gLV_with_extinction(n, maxtime, interactions, ri, starting_abundances, sigma):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        steady = True
        while i < len(livespecies):
            try:
                change_per_capita = ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]])
                if change_per_capita != 0 and steady:
                    steady = False
                dW = np.random.random() - np.random.random()
                new_abundance = abundances[time-1][livespecies[i]] + int(abundances[time-1][livespecies[i]]*change_per_capita) + int(np.sqrt(abundances[time-1][livespecies[i]]*sigma[livespecies[i]])*dW)
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
