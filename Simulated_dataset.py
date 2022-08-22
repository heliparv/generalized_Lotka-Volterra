from code import interact
from timeit import repeat
import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_abundances, adjust_selfinteractions_for_carrying_capacity, add_sparcity
from Cyclic_games import generalized_rps, even_groups_for_rps
from simulations_simple_gLV import simple_gLV_with_extinction
from graphics import abundances_line_chart
from data_modification import relative_abundances


def simulated_datasets_abundances(repetition=20,n=23,maxtime=250,time_increment=0.05,ri_mean=0.1,ri_std=0.3,carrying_capacities_mean=2000,carrying_capacities_std=0.5,abundances_mean_range_low=50,abundances_mean_range_up=500,abundances_std=0.1,interactions_mean=-0.00003,interactions_std=0.000002,interactions_sparcity=0,rps_groups=7,rps_distance=2,rps_sparcity=0.4,rps_interactions_std=0.2,off_target_interactions = False):
    rep=0
#Here we want to change the abundance mean for every repetition so we can have different values
    interact_seed=np.random.randint(0, 1000)
    community_seed=np.random.randint(0, 1000)

    np.random.seed(community_seed)
    ri_seed = np.random.randint(0, 2**31)
    carrying_capacities_seed = np.random.randint(0, 2**31)
   
    ri = generate_growth_rates(n, ri_mean, ri_seed, ri_std)
    carrying_capacities = generate_abundances(n, carrying_capacities_seed, carrying_capacities_mean, carrying_capacities_std)

    np.random.seed(interact_seed)      
    interactions_seed = np.random.randint(0, 2**31)
    interactions_sparcity_seed = np.random.randint(0, 2**31)
    rps_group_seed = np.random.randint(0, 2**31)
    rps_sparcity_seed = np.random.randint(0, 2**31)
    rps_interactions_seed = np.random.randint(0, 2**31)

    pairwise_interactions = generate_interactions(n, 1, interactions_seed, interactions_mean, interactions_std)
    pairwise_interactions = add_sparcity(pairwise_interactions, interactions_sparcity, interactions_sparcity_seed)
    pairwise_interactions = generalized_rps(pairwise_interactions, rps_groups, rps_distance, even_groups_for_rps, rps_group_seed, rps_sparcity_seed, rps_sparcity, rps_interactions_seed, rps_interactions_std, off_target_interactions)
    pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)

    interact_columns = ["species", "ri", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23"]
    interact_array = np.zeros((23,25))
    interact_array[:,0] = list(range(1,24))
    interact_array[:,1] = ri
    interact_array[:,2:] = pairwise_interactions
    interact_frame = pd.DataFrame(interact_array, columns=interact_columns)
    
    #community_seeds=np.random.sample(range(0, 1000), repetition)
    #interact_seed=np.random.sample(range(0, 1000), repetition)
    while True:
        save = input("save? y/n")
        if save == "y":
            
#Create a dataframe with as column the experiment number and the species ID
            species_columns = ["Experiment"]+list(range(1, n+1))
            data = pd.DataFrame(columns=species_columns)
            abundance_seeds=[]
            abundance_means=[]
            while rep<repetition:
                np.random.randint(0,1000)
                abundances_seed=np.random.randint(0, 2**31)
                abundance_seeds.append(abundances_seed)
                abundances_mean=np.random.randint(abundances_mean_range_low,abundances_mean_range_up)
                abundance_means.append(abundances_mean)
                starting_abundances = generate_abundances(n,abundances_seed,abundances_mean,abundances_std)
                abundances = simple_gLV_with_extinction(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)
                
                
                if type(abundances) == int:
                    print("Encountered error")
                else:
                    if sum(abundances[-1]<10)>0:
                        print("nope")

                    expframe = pd.DataFrame(np.reshape(([rep]+list(abundances[-1])), (1,n+1)), columns=species_columns)
                    data=pd.concat([data, expframe], ignore_index=True, axis=0)

                    relative_data = relative_abundances(np.array(data.iloc[:,1:]))
                    relative_data = np.insert(relative_data,0, data.iloc[:,0], axis=1)
                    relative = pd.DataFrame(relative_data, columns=species_columns)

                    rep+=1
            parameters = np.reshape([n, maxtime, time_increment, community_seed, abundance_seeds, interact_seed, ri_mean, ri_std, carrying_capacities_mean, carrying_capacities_std, abundance_means, abundances_std, interactions_mean, interactions_std, interactions_sparcity, rps_groups, rps_distance, rps_sparcity, rps_interactions_std, off_target_interactions], (1, 20))
            param_columns = ["n", "maxtime", "time_increment", "community_seed", "abundances_seeds", "interactions_seed", "ri_mean", "ri_std", "carrying_capacities_mean", "carrying_capacities_std", "abundances_mean", "abundances_std", "interactions_mean", "interactions_std", "interactions_sparcity","rps_groups", "rps_distance", "rps_sparcity", "rps_interactions_std", "off_target_interactions"]
            param_frame = pd.DataFrame(parameters, columns=param_columns)
            
            with pd.ExcelWriter("results_abundances.xlsx") as writer:
                param_frame.to_excel(writer, index=False, sheet_name="Parameters")
                interact_frame.to_excel(writer, index=False, sheet_name="Interactions")
                data.to_excel(writer, index=False, sheet_name="Experiments abundances")
                relative.to_excel(writer, index=False, sheet_name="Relative abundances")
            break
        elif save == "n":
            break


def simulated_datasets_ri(repetition=20,n=23,maxtime=250,time_increment=0.05,ri_mean_low=0.0, ri_mean_up=0.5,ri_std=0.3,carrying_capacities_mean=2000,carrying_capacities_std=0.5,abundances_seed=55,abundances_mean=100,abundances_std=0.1,interactions_mean=-0.00003,interactions_std=0.000002,interactions_sparcity=0,rps_groups=7,rps_distance=2,rps_sparcity=0.4,rps_interactions_std=0.2,off_target_interactions = False):
    rep=0
#Here we want to change the ri mean for every repetition so we can have different values
    interact_seed=np.random.randint(0, 1000)
    community_seed=np.random.randint(0, 1000)

    np.random.seed(community_seed)
    carrying_capacities_seed = np.random.randint(0, 2**31)
   

    carrying_capacities = generate_abundances(n, carrying_capacities_seed, carrying_capacities_mean, carrying_capacities_std)

    np.random.seed(interact_seed)      
    interactions_seed = np.random.randint(0, 2**31)
    interactions_sparcity_seed = np.random.randint(0, 2**31)
    rps_group_seed = np.random.randint(0, 2**31)
    rps_sparcity_seed = np.random.randint(0, 2**31)
    rps_interactions_seed = np.random.randint(0, 2**31)

    pairwise_interactions = generate_interactions(n, 1, interactions_seed, interactions_mean, interactions_std)
    pairwise_interactions = add_sparcity(pairwise_interactions, interactions_sparcity, interactions_sparcity_seed)
    pairwise_interactions = generalized_rps(pairwise_interactions, rps_groups, rps_distance, even_groups_for_rps, rps_group_seed, rps_sparcity_seed, rps_sparcity, rps_interactions_seed, rps_interactions_std, off_target_interactions)

    while True:
        save = input("save? y/n")
        if save == "y":
            
#Create a dataframe with as column the experiment number and the species ID
            species_columns = ["Experiment"]+list(range(1, n+1))
            data = pd.DataFrame(columns=species_columns)
            ri_seeds=[]
            ri_means=[]
            interactions=[]
            names=[]
            while rep<repetition:
                np.random.randint(0,1000)
                ri_seed = np.random.randint(0, 2**31)

                ri_seeds.append(ri_seed)
                ri_mean=np.random.uniform(ri_mean_low,ri_mean_up)
                ri_means.append(ri_mean)
                ri = generate_growth_rates(n, ri_mean, ri_seed, ri_std)
                pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)
                interact_columns = ["species", "ri", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23"]
                interact_array = np.zeros((23,25))
                interact_array[:,0] = list(range(1,24))
                interact_array[:,1] = ri
                interact_array[:,2:] = pairwise_interactions
                interact_frame = pd.DataFrame(interact_array, columns=interact_columns)
                interactions.append(interact_frame)
                starting_abundances = generate_abundances(n,abundances_seed,abundances_mean,abundances_std)
                abundances = simple_gLV_with_extinction(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)
    
                
                if type(abundances) == int:
                    print("Encountered error")
                else:
                    if sum(abundances[-1]<10)>0:
                        print("nope")

                    expframe = pd.DataFrame(np.reshape(([rep]+list(abundances[-1])), (1,n+1)), columns=species_columns)
                    data=pd.concat([data, expframe], ignore_index=True, axis=0)

                    relative_data = relative_abundances(np.array(data.iloc[:,1:]))
                    relative_data = np.insert(relative_data,0, data.iloc[:,0], axis=1)
                    relative = pd.DataFrame(relative_data, columns=species_columns)

                    names.append("Interactions_exp. "+str(rep))
                    rep+=1
                    
            parameters = np.reshape([n, maxtime, time_increment, community_seed, abundances_seed,ri_seeds, interact_seed, ri_means, ri_std, carrying_capacities_mean, carrying_capacities_std,abundances_mean, abundances_std, interactions_mean, interactions_std, interactions_sparcity, rps_groups, rps_distance, rps_sparcity, rps_interactions_std, off_target_interactions], (1, 21))
            param_columns = ["n", "maxtime", "time_increment", "community_seed", "abundances_seed","ri_seeds", "interactions_seed", "ri_means", "ri_std", "carrying_capacities_mean", "carrying_capacities_std","abundance_mean", "abundances_std", "interactions_mean", "interactions_std", "interactions_sparcity","rps_groups", "rps_distance", "rps_sparcity", "rps_interactions_std", "off_target_interactions"]
            param_frame = pd.DataFrame(parameters, columns=param_columns)
            
            
            with pd.ExcelWriter("results_ri.xlsx") as writer:
                param_frame.to_excel(writer, index=False, sheet_name="Parameters")
                #interact_frame.to_excel(writer, index=False, sheet_name="Interactions")
                data.to_excel(writer, index=False, sheet_name="Experiments abundances")
                relative.to_excel(writer, index=False, sheet_name="Relative abundances")
                [A.to_excel(writer,sheet_name="{0}".format(names[i])) for i, A in enumerate(interactions)]
            break
        elif save == "n":
            break


def simulated_datasets_ri_abundances(repetition=20,n=23,maxtime=250,time_increment=0.05,ri_mean_low=0.0, ri_mean_up=0.5,ri_std=0.3,carrying_capacities_mean=2000,carrying_capacities_std=0.5,abundances_mean_range_low=50,abundances_mean_range_up=500,abundances_std=0.1,interactions_mean=-0.00003,interactions_std=0.000002,interactions_sparcity=0,rps_groups=7,rps_distance=2,rps_sparcity=0.4,rps_interactions_std=0.2,off_target_interactions = False):
    rep=0
#Here we want to change the ri and abundance mean for every repetition so we can have different values
    interact_seed=np.random.randint(0, 1000)
    community_seed=np.random.randint(0, 1000)

    np.random.seed(community_seed)
    carrying_capacities_seed = np.random.randint(0, 2**31)
   

    carrying_capacities = generate_abundances(n, carrying_capacities_seed, carrying_capacities_mean, carrying_capacities_std)

    np.random.seed(interact_seed)      
    interactions_seed = np.random.randint(0, 2**31)
    interactions_sparcity_seed = np.random.randint(0, 2**31)
    rps_group_seed = np.random.randint(0, 2**31)
    rps_sparcity_seed = np.random.randint(0, 2**31)
    rps_interactions_seed = np.random.randint(0, 2**31)

    pairwise_interactions = generate_interactions(n, 1, interactions_seed, interactions_mean, interactions_std)
    pairwise_interactions = add_sparcity(pairwise_interactions, interactions_sparcity, interactions_sparcity_seed)
    pairwise_interactions = generalized_rps(pairwise_interactions, rps_groups, rps_distance, even_groups_for_rps, rps_group_seed, rps_sparcity_seed, rps_sparcity, rps_interactions_seed, rps_interactions_std, off_target_interactions)

    while True:
        save = input("save? y/n")
        if save == "y":
            
#Create a dataframe with as column the experiment number and the species ID
            species_columns = ["Experiment"]+list(range(1, n+1))
            data = pd.DataFrame(columns=species_columns)
            ri_seeds=[]
            ri_means=[]
            abundance_seeds=[]
            abundance_means=[]
            interactions=[]
            names=[]
            while rep<repetition:
                np.random.randint(0,1000)
                abundances_seed=np.random.randint(0, 2**31)
                abundance_seeds.append(abundances_seed)
                abundances_mean=np.random.randint(abundances_mean_range_low,abundances_mean_range_up)
                abundance_means.append(abundances_mean)

                np.random.randint(0,1000)
                ri_seed = np.random.randint(0, 2**31)
                ri_seeds.append(ri_seed)
                ri_mean=np.random.uniform(ri_mean_low,ri_mean_up)
                ri_means.append(ri_mean)
                ri = generate_growth_rates(n, ri_mean, ri_seed, ri_std)
                pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)
                interact_columns = ["species", "ri", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23"]
                interact_array = np.zeros((23,25))
                interact_array[:,0] = list(range(1,24))
                interact_array[:,1] = ri
                interact_array[:,2:] = pairwise_interactions
                interact_frame = pd.DataFrame(interact_array, columns=interact_columns)
                interactions.append(interact_frame)
                starting_abundances = generate_abundances(n,abundances_seed,abundances_mean,abundances_std)
                abundances = simple_gLV_with_extinction(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)
    
                
                if type(abundances) == int:
                    print("Encountered error")
                else:
                    if sum(abundances[-1]<10)>0:
                        print("nope")

                    expframe = pd.DataFrame(np.reshape(([rep]+list(abundances[-1])), (1,n+1)), columns=species_columns)
                    data=pd.concat([data, expframe], ignore_index=True, axis=0)

                    relative_data = relative_abundances(np.array(data.iloc[:,1:]))
                    relative_data = np.insert(relative_data,0, data.iloc[:,0], axis=1)
                    relative = pd.DataFrame(relative_data, columns=species_columns)

                    names.append("Interactions_exp. "+str(rep))
                    rep+=1
                    
            parameters = np.reshape([n, maxtime, time_increment, community_seed, abundance_seeds, abundance_means, ri_seeds,ri_means, ri_std, interact_seed,carrying_capacities_mean, carrying_capacities_std,abundances_mean, abundances_std, interactions_mean, interactions_std, interactions_sparcity, rps_groups, rps_distance, rps_sparcity, rps_interactions_std, off_target_interactions], (1, 22))
            param_columns = ["n", "maxtime", "time_increment", "community_seed", "abundances_seed","abundance_means","ri_seeds",  "ri_means", "ri_std", "interactions_seed","carrying_capacities_mean", "carrying_capacities_std","abundance_mean", "abundances_std", "interactions_mean", "interactions_std", "interactions_sparcity","rps_groups", "rps_distance", "rps_sparcity", "rps_interactions_std", "off_target_interactions"]
            param_frame = pd.DataFrame(parameters, columns=param_columns)
            
            
            with pd.ExcelWriter("results_ri_abundances.xlsx") as writer:
                param_frame.to_excel(writer, index=False, sheet_name="Parameters")
                #interact_frame.to_excel(writer, index=False, sheet_name="Interactions")
                data.to_excel(writer, index=False, sheet_name="Experiments abundances")
                relative.to_excel(writer, index=False, sheet_name="Relative abundances")
                [A.to_excel(writer,sheet_name="{0}".format(names[i])) for i, A in enumerate(interactions)]
            break
        elif save == "n":
            break


simulated_datasets_abundances()
simulated_datasets_ri()
simulated_datasets_ri_abundances()
