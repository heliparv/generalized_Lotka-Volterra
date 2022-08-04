import matplotlib.pyplot as plt
import numpy as np

""" Tools for graphic representations of the simulated data.

functions:
abundances_line_chart: Draws a line chart of the given abundance data.

interactions heatmap: Draws a heatmap based on the given interactions between species.
"""

def abundances_line_chart(n, maxtime, time_increment, abundances):
    x = np.arange(0, maxtime+time_increment, time_increment)
    plt.figure()
    plt.plot(x, abundances)
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.legend(range(1,n+1), title="Species", loc="center left", fontsize='small', bbox_to_anchor=(1, 0.5))
    return plt.show()

def abundances_and_nutrient_chart(n, maxtime, time_increment, abundances, nutrient):
    x = np.arange(0, maxtime+time_increment, time_increment)
    fig, ax1 = plt.subplots() 
    ax1.set_xlabel('') 
    ax1.set_ylabel('Species abundance') 
    plot_1 = ax1.plot(x, abundances) 
    ax1.tick_params(axis ='y')
    ax2 = ax1.twinx()   
    ax2.set_ylabel("Nutrient")
    plot_2 = ax2.plot(x, nutrient) 
    ax2.tick_params(axis ='y')
    return plt.show()

def interactions_heatmap(n, interactions):
    plt.imshow(interactions, cmap="hot", interpolation="nearest")
    plt.ylabel("Species affected")
    plt.yticks(range(0,n), range(1, n+1))
    plt.xlabel("Effector species")
    plt.xticks(range(0,n), range(1, n+1))
    return plt.show()
