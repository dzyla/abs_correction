# abs_correction
A small python script for Scattering correction for protein concentration determination using UV absorbance (200 - 350 nm)


it takes a csv file (check sample file) with only two columns and calculates the protein concentration @280 nm corrected for light scattering (baseline is fitted in log-log plot for 350 - 310 nm). Plots the data and baseline on log-log plot. User has to know the pathlenght and extinction coefficient. 

It works for now only with 2 column data and 2 lines of headers. Few packages are required: scipy, numpy and matplotlib. Python 3+
