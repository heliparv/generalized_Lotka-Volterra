import matplotlib.pyplot as plt
import numpy as np

""" Tools for graphic representations of the simulated data.

functions:
abundances_line_chart: Draws a line chart of the given abundance data.

interactions heatmap: Draws a heatmap based on the given interactions between species.
"""

def abundances_line_chart(n, maxtime, abundances):
    y = range(0, maxtime+1)
    plt.figure()
    plt.plot(abundances, y)
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.legend(range(1,n+1), title="Species", loc="center left", fontsize='small', bbox_to_anchor=(1, 0.5))
    return plt.show()

def interactions_heatmap(n, interactions):
    plt.imshow(interactions, cmap="hot", interpolation="nearest")
    plt.ylabel("Species affected")
    plt.yticks(range(0,n), range(1, n+1))
    plt.xlabel("Effector species")
    plt.xticks(range(0,n), range(1, n+1))
    return plt.show()