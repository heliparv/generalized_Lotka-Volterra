'''Methods for creating data output frames and sampling data 

relative_abundances: Input is an array of absolute abundances of different
species at given time points. Returns an array with relative abundances

initiate_df: initiates dataframe for data output. The same function can be used for
both absolute abundance frames and relative abundance frames.

addable_to_frame: takes given data and formats it into a pandas dataframe that can be
added to the total dataframe

output: takes given dataframes and outputs it into an excel file
'''
import numpy as np
import pandas as pd

def relative_abundances(abundances):
    sums = np.sum(abundances, axis=1)
    return abundances/sums[:,None]

def initiate_df(n):
    columns = ["experiment", "time"]+list(range(n))
    frame = pd.DataFrame(columns=columns)
    return columns, frame

def addable_to_frame(columns, sim_no, times, data):
    sim = sim_no*np.ones(len(times))
    data = np.insert(data, 0, times, axis=1)
    data = np.insert(data, 0, sim, axis=1)
    return pd.DataFrame(data, columns=columns)

def output(output_name, param_df, abs_df, rel_df):
    with pd.ExcelWriter(f'{output_name}.xlsx') as writer:
        param_df.to_excel(writer, sheet_name='Parameters', index=False)
        abs_df.to_excel(writer, sheet_name='Absolute abundances', index=False)
        rel_df.to_excel(writer, sheet_name='Relative abundances', index=False)
