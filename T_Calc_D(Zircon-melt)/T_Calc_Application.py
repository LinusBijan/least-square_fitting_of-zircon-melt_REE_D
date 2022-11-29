# -*- coding: utf-8 -*-
"""
Created on Nov 2021
@author: Linus Streicher 
Vrije Universiteit Amsterdam, Faculty of Science, Geology and Geochemistry Cluster
For questions and remarks write to 
linusbijan@yahoo.de 


With this Code you can calculate crystallisation temperatures of zircons from REE zircon-melt/bulk-rock partition coefficents.
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# constants 
constraint   = 6489e-21        # Q1
#constraint   = 7827-21       # Q2
N_a = 6.02214086 * 10**23;     # 1/mol
R   = 8.314472;                # J/(mol*K)
r_0 = 0.94                     # Angström

# import inital data 
initial_data = pd.read_excel("Partition_coefficents_data.xlsx")
# get dataset names
column_names = []
for col in initial_data.columns:
    column_names.append(col)
column_names = column_names[1:]
# create empty table for results
results    = pd.DataFrame(columns=['Sample','Temperature [K]','err_T'])
# Independent variables 
indi_var =  (-4*np.pi*N_a*constraint)/(R)
# model function
def model_constrainted(x, T):
    exponent = ((indi_var/T)*(0.5 * r_0 * (x - r_0) ** 2 + (1 / 3) * (x - r_0) ** 3))/((1.38 + r_0)**3)
    #factor   = np.exp((22420/T)-14.221)    # Rubatto and Hermann (2007)
    factor   = np.exp((13594/T)-7.1266)     # Streicher et al. (2022)
    return factor * np.exp(exponent)

# curve fit to determine T 
for a in np.arange(0,len(column_names)):
    y            = initial_data[column_names[a]].to_numpy()
    x            = initial_data['Radii'].to_numpy()
    
    popt, pcov = curve_fit(model_constrainted, x, y,
          maxfev  = 10000,
          bounds  = ((0),(2000)),
          # if using bounds switch method to 'trf' or 'dogbox'. If not, also method to 'lm' works 
          method='trf');
    
    # Array with errors for D0,r0 and E
    perr = np.sqrt(np.diag(pcov))
    
    res_help = pd.DataFrame([[column_names[a], popt, perr,]],
                        columns=['Sample','Temperature [K]','err_T'])
    results = pd.concat([res_help, results])
