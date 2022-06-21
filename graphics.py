import matplotlib.pyplot as plt
import numpy as np
from sqlalchemy import null

""" Tools for graphic representations of the simulated data.

functions:
abundances_line_chart: Draws a line chart of the given abundance data.
"""

#TODO add function for a heatmap of species interactions

def abundances_line_chart(n, maxtime, abundances):
    y = range(0, maxtime+1)
    plt.figure()
    plt.plot(abundances, y)
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.legend(range(1,n+1), title="Species", loc="center left", fontsize='small', bbox_to_anchor=(1, 0.5))
    return plt.show()
