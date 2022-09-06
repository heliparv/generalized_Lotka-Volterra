import numpy as np

''' Functions for performing simulations with generalized Lotka-Volterra dynamics
using the fourth order Runge-Kutta method to find the approximate solutions.
All of the functions model extinction by dropping bacteria species from the simulation
when their abundance reaches 0 or below on one round of simulation. In case of error
each simulation returns "-1"

find_change: A function that finds the magnitude of change for given starting abundances
for a time increment. Used for calculating k1-k4 for the Runge-Kutta method

add_steady_state_tail: Function adds steady state abundances from given timepoint in
deterministic functions.

gLV_RK: function simulates generalized Lotka-Volterra dynamics with a fourth order
Runge-Kutta method. Function is deterministic.

stochastic_gLV_RK: Same as previous function, but adds a random term for each timepoint
so that there's stochasticity in the changerate for the bacteria.

leave_one_out_gLV_RK: Function takes a set of bacteria species and returns a set of
leave-one-out experiments for the species set.

gLV_RK_return_steady_state: Only returns the steady state of a deterministic generalized
Lotka Volterra simulation, or the last simulated values if steady state hasn't been found.
Also returns a variable that states if steady state was found.
'''

def find_change(n, abundances, interactions, ri, livespecies, time_increment):
        return_change = np.zeros(n)
        i = 0
        while i < len(livespecies):
            try:
                change_per_capita = (ri[livespecies[i]] + sum(abundances*interactions[livespecies[i]]))*time_increment
                change = abundances[livespecies[i]]*change_per_capita
            except:
                return -1
            return_change[livespecies[i]] = change
            i+=1
        return return_change

def add_steady_state_tail(abundances, time, max_increments):
    steady_abundances = np.repeat(np.array([abundances[time]]), max_increments-time, axis=0)
    abundances[time+1:] = steady_abundances
    return abundances

def gLV_RK(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    if sum(starting_abundances==0)>0:
        i = 0
        while i < len(livespecies):
            if starting_abundances[livespecies[i]] == 0:
                del livespecies[i]
            i +=1
    steady = False
    while time < max_increments+1:
        k1 = find_change(n, abundances[time-1], interactions, ri, livespecies, time_increment)
        if type(k1) == int:
            return -1
        k2 = find_change(n, abundances[time-1]+(k1/2), interactions, ri, livespecies, time_increment)
        if type(k2) == int:
            return -1
        k3 = find_change(n, abundances[time-1]+(k2/2), interactions, ri, livespecies, time_increment)
        if type(k3) == int:
            return -1
        k4 = find_change(n, abundances[time-1]+k3, interactions, ri, livespecies, time_increment)
        if type(k4) == int:
            return -1
        change = np.around((1/6)*(k1+2*k2+2*k3+k4), decimals=6)
        if sum(change) == 0:
            steady = True
            break
        abundances[time] = abundances[time-1]+change
        i = 0
        while i < len(livespecies):
            if abundances[time][livespecies[i]] <= 0:
                abundances[time][livespecies[i]] == 0
                del livespecies[i]
            i +=1
        time += 1
    if steady and time < max_increments:
        return add_steady_state_tail(abundances, time, max_increments)
    return abundances

def stochastic_gLV_RK(n, maxtime, time_increment, interactions, ri, starting_abundances, sigma):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    if sum(starting_abundances==0)>0:
        i = 0
        while i < len(livespecies):
            if starting_abundances[livespecies[i]] == 0:
                del livespecies[i]
            i +=1
    while time < max_increments+1:
        k1 = find_change(n, abundances[time-1], interactions, ri, livespecies, time_increment)
        if type(k1) == int:
            return -1
        k2 = find_change(n, abundances[time-1]+(k1/2), interactions, ri, livespecies, time_increment)
        if type(k2) == int:
            return -1
        k3 = find_change(n, abundances[time-1]+(k2/2), interactions, ri, livespecies, time_increment)
        if type(k3) == int:
            return -1
        k4 = find_change(n, abundances[time-1]+k3, interactions, ri, livespecies, time_increment)
        if type(k4) == int:
            return -1
        change = np.around((1/6)*(k1+2*k2+2*k3+k4), decimals=6)
        stochasticity = np.zeros(n)
        i = 0
        while i < len(livespecies):
            dW = np.random.random() - np.random.random()
            stochasticity[livespecies[i]] = np.sqrt(abundances[time-1][livespecies[i]]*sigma[livespecies[i]])*dW
            i += 1
        abundances[time] = abundances[time-1]+change+stochasticity
        i = 0
        while i < len(livespecies):
            if abundances[time][livespecies[i]] <= 0:
                abundances[time][livespecies[i]] == 0
                del livespecies[i]
            i +=1
        time += 1
    return abundances

def gLV_RK_return_steady_state(n, maxtime, time_increment, interactions, ri, starting_abundances):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((2, n))
    abundances[0] = starting_abundances
    index = 0
    time = 1
    livespecies = list(range(0,n))
    if sum(starting_abundances==0)>0:
        i = 0
        while i < len(livespecies):
            if starting_abundances[livespecies[i]] == 0:
                del livespecies[i]
            i +=1
    steady = False
    while time < max_increments+1:
        k1 = find_change(n, abundances[index], interactions, ri, livespecies, time_increment)
        if type(k1) == int:
            return -1
        k2 = find_change(n, abundances[index]+(k1/2), interactions, ri, livespecies, time_increment)
        if type(k2) == int:
            return -1
        k3 = find_change(n, abundances[index]+(k2/2), interactions, ri, livespecies, time_increment)
        if type(k3) == int:
            return -1
        k4 = find_change(n, abundances[index]+k3, interactions, ri, livespecies, time_increment)
        if type(k4) == int:
            return -1
        change = np.around((1/6)*(k1+2*k2+2*k3+k4), decimals=6)
        if sum(change) == 0:
            steady = True
            break
        abundances[abs(index-1)] = abundances[index]+change
        i = 0
        while i < len(livespecies):
            if abundances[time][livespecies[i]] <= 0:
                abundances[time][livespecies[i]] == 0
                del livespecies[i]
            i +=1
        if steady:
            break
        else:
            index = abs(index-1)
            time += 1
    return steady, abundances[index]
