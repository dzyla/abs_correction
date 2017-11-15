# abs_correction
A small python script to correct absorbance for scattering

it takes a csv file with only two columns and calculates the protein concentration corrected for light scattering (baseline is fitted in log-log plot for 350 - 310 nm).

It works for now only with 2 column data and 2 lines of headers. Few packages are required: scipy, numpy and matplotlib. Python 3+
