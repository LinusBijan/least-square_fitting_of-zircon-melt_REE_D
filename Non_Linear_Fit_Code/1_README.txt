Manual:
0	Save the all folders in the same directory
1. 	Open the excel file "2_Lattice_Strain_Data"
2. 	Fill in the cells highlighted in green (Di = measured partition coefficient of element i, 
	ri = ionic radius of element i, 1 s% to calculate 1s (abs) automatically or just 1s as well 
	as T [°C] and P, under which the temperature in °C and pressure in GPa 
3.	Run the Code with Python (e.g., via Spyder) 

To change the style of the plot, the Code has to be modified (See section ##Plotting at the bottom of the Code)
Keep attention that no spaces are in the cells of "2_Lattice_Strain_Data" ; this could result in trouble. 
The Code will not work if the folder structure will be changed


The code 3_Non_Linear_Regression.py works for all regressions to the lattice strain model just change the ionic radius ri for the elements. The  4_Non_Linear_Regression_constrained.py works for zircon only.

The data in the test file is from Burnham, A.D., Berry, A.J., 2012. An experimental study of trace element partitioning between zircon and melt as a function of oxygen fugacity. Geochimica et Cosmochimica Acta 95, 196-212.

