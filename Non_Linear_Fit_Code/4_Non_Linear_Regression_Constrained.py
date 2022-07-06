# -*- coding: utf-8 -*-
"""
Created on Nov 2020

@author: Linus Streicher 

Vrije Universiteit Amsterdam, Faculty of Science, Geology and Geochemistry Cluster


For questions and remarks write to 

linusbijan@yahoo.de 

"""

"""
This python code is made to fit partitioning coefficents of REE and Sc in 
zircon based on a predictive model by Blundy and Wood, 1994 & 1997, using the Lattice Strain Model by Brice,1975  
"""

import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.optimize import curve_fit

plt.close('all')

# Import data from excel file

"""intial_data enter a .xlsx - file with the same structur as the .xlsx in the same folder of this code. It is important that the xlsx - file 
used and the program code are in the same folder   """

initial_data = pd.read_excel("2_Lattice_Strain_Data.xlsx") 

# Select data from Dataframe columns and transform columns to numpy lists
names = initial_data['Element'].to_numpy()
# Get ionic radii
x     = initial_data['ri'].to_numpy()          # Angstroem (10**(-10 m))
# Get partition coefficents and error 
y     = initial_data['Di'].to_numpy()
y_err = initial_data['1s'].to_numpy()
# Get temperature
T     = initial_data['T [C]'].to_numpy()
T     = T[0] + 273.15;                         # K
# Get pressure
P     = initial_data['P [GPa]'].to_numpy()
P     = P[0]                                   # GPa

# Constants
N_a = 6.02214086 * 10**23;                     # 1/mol
R   = 8.314472;                                # J/(K*mol)

### Calculation
# Constant Q / constraint factor Q1 = 6489, Q2 = 7687
Q = 6489
# Q = 7687
constraint_constant = (Q)*1e-21 # J
indi_var          = (-4*np.pi*N_a*constraint_constant)/R * 1/T 

# Define odel function to which the data will be Fit 
# Including constraint
def model(x, D_0, r_0):
    exponent = (indi_var*(0.5 * r_0 * (x - r_0) ** 2 + (1 / 3) * (x - r_0) ** 3))/((1.38 + r_0)**3)
    return D_0 * np.exp(exponent)

# Curve fitting
"""
Additional Information about the fitting methods:

‘trf’   : Trust Region Reflective algorithm, particularly suitable for large sparse problems with bounds. Generally robust method.

‘dogbox’ : dogleg algorithm with rectangular trust regions, typical use case is small problems with bounds. Not recommended for 
           problems with rank-deficient Jacobian.

‘lm’    : Levenberg-Marquardt algorithm as implemented in MINPACK. Doesn’t handle bounds and sparse Jacobians. Usually the most efficient 
           method for small unconstrained problems.

for more Information see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares

"""
             
popt, pcov = curve_fit(model, x, y,
                      sigma=y_err,
                      # set this on "True" if errors are real errors on y_err. If "False" the errors are only used as relative weights
                      absolute_sigma =True,
                      # Number of iterations
                      maxfev  = 50,
                      # if using bounds switch method to 'trf' or 'dogbox'. If not, also method to 'lm' works 
                      bounds  = ((0,0),(2000,2)),
                      method='dogbox');

# Array with errors for D0 and r0
perr = np.sqrt(np.diag(pcov))

### Prepare data for plotting
# Factor including the constants 
D0_out        = round(popt[0],2);
D0_Sig        = round(np.absolute(perr[0]),2);
r0_out        = round(popt[1],3);
r0_Sig        = round(np.absolute(perr[1]),3);
# Calculate E via Q
E_out         = round(((constraint_constant/((1.38
                                            +popt[1])*10**(-10))**3)/1e9)
                      ,1                      );
E_out         = round(((constraint_constant/((1.38
                                            +popt[1])*1e-10)**3)/1e9)   # GPa, *10**(-10) conversion from Angstroem in Meters
                      ,1);  
# Taylor-Polynom 2nd degree to approximate the error of E
E_Sig  = round(np.sqrt((((-3*constraint_constant/((1.38+popt[1])*1e-10)**4)*perr[1]*1e-10 + (6*constraint_constant/((1.38+popt[1])*1e-10)**5)*(perr[1]*1e-10)**2)/1e9)**2),1);
        
### Plotting Section
# helper array to plot curve
x_help   = np.array(np.arange(0.8,1.2,step=0.01));
# Plotting Data and Curve, Figsize in px/100 
fig, ax = plt.subplots(figsize=(14.5,13))

# Plots a Graph in the Range of x_help
plt.plot(x_help, model(x_help, *popt),
# color: 'g' = green,'b' = blue, 'r' = red, 'c'= cyan, 'm' = magenta, 'y' = yellow, 'k' = black, 'w' = white 
          color = 'r')     
# Plots data points with errorbar
plt.errorbar(x, y, yerr=y_err, fmt='o',
             color = 'k')

# Log on y axis in numbers 
plt.yscale("log")
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))

# Labels axis
plt.xlabel('Ionic radius [\u00C5]',
            fontsize = 14,
            )
plt.ylabel('D$_{(Zircon/Melt)}$',
           fontsize = 14
           )

# Limits axis
plt.xlim(left=0.8, right=1.3)
plt.ylim(bottom=0.001, top=10)

# # Label data points with a 0.5 % offset
# for j, txt in enumerate(names):
#     ax.annotate(txt, 
#                 (x[j], y[j]), 
#                 xytext=(x[j]+x[j]*0.005, 
#                         y[j]+y[j]*0.005))
# # Adds grid
# ax.grid(True)

# Text box with results including errors // To change the last digit change the number before "f"
anchored_text_1 = AnchoredText("D$_{0}$"+ " = " 
                              + str('{:.2f}'.format(D0_out)) + " \u00B1 " 
                              + str('{:.2f}'.format(D0_Sig)) + "\n" + "r$_{0}$" + " = " 
                              + str('{:.3f}'.format(r0_out)) + " \u00B1 " 
                              + str('{:.3f}'.format(r0_Sig)) + " \u00C5" + "\n" + "E" + " = " 
                              + str('{:.0f}'.format(E_out)) + " \u00B1 " 
                              + str('{:.0f}'.format(E_Sig)) + " GPa",
                               prop=dict(fontsize = 11),
# Number changes location of text box in graph
                              loc=3)                                    
anchored_text_2 = AnchoredText(str('{:.0f}'.format(P)) + " GPa",
                              prop=dict(fontsize = 11),
                              loc=1)
ax.add_artist(anchored_text_1)
ax.add_artist(anchored_text_2)
# Print Plot 
plt.show()
