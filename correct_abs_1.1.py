'''script written by Dawid Zyla for automatic scattering correction of Absorbance spectrum.
work with files directly exported from Varian Cary as *.csv. Applies linear regression in 
log-log plot between log(350) and log(310).
To run it you need umpy, matplotlib and scipy.
Use is inside the directory with csv files.
example of use in powershell with python 3.5:
C:\Python35\python.exe .\correct_abs.py
v1.1 - files with less or more than 3 columns (x,y and '') will not be taken for calculations (added file_check)
'''
import numpy as np
import matplotlib.pyplot as pyplot
import math
from scipy import stats
import os
import csv
from sys import exit
from time import sleep


corrected_abs = 0
abs280 = 0
protein_conc = 0
r_value = 0
skip = 0

def generate_plot(file):
    data = np.genfromtxt(str(file).replace('.csv','')+'_log.csv', delimiter=',', skip_header=2, names=['x', 'y'])
    #print(data)
    data_baseline = np.genfromtxt('baseline.csv', delimiter=',', names=['x', 'y'])
    pyplot.xlabel ('log(Wavelength)')
    pyplot.ylabel ('log(Absorbance)')
    pyplot.plot(data['x'], data['y'], label='log(Absorbance)')
    pyplot.plot(data_baseline['x'], data_baseline['y'], label='log(Baseline)')
    pyplot.legend(loc=1)  # upper right corner, anticlockwise
    pyplot.suptitle(file, size=15) #titlesize and name
    pyplot.annotate('Protein concentration:  '+str(round(protein_conc,2))+' uM\n$R^{2}$ = '+str(round(-r_value,3)), xy=(math.log10(280), math.log10(abs280)),
                xytext=(math.log10(280)-0.05, math.log10(float(abs280)+0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.05))
    #axes = pyplot.gca() #for axes adjustment
    #axes.set_ylim(ymax=0) #set the y range
    pyplot.show()
    #dodac export wykresu!


def log_log(file): #generation of log-log plot and 350-310 nm for baseline correction
    data = np.genfromtxt(file, delimiter=',', skip_header=2, names=['x', 'y'])
    with open(file, 'r', newline='') as f:
        with open(str(file).replace('.csv','')+'_log.csv', 'w', newline='') as f_out:
            with open('log_fit.tmp', 'w', newline='') as flog_out:
            #writer = csv.writer(f_out,delimiter=',')
                for i, line in enumerate(data):
                    x_log = math.log10(float(line[0]))
                    y_log = math.log10(float(math.sqrt(line[1]**2)))
                    data = (str(x_log)+','+str(y_log))
                    if float(line[0]) <= 350 and float(line[0]) >= 310:
                        print(data, file=flog_out)
                    print(data,file=f_out)
                    #print(data)

def fitting(file):
    data_log = np.genfromtxt("log_fit.tmp", delimiter=',', skip_header=2, names=['x', 'y']) #for baseline
    global r_value
    slope, intercept, r_value, p_value, std_err = stats.linregress(data_log['x'], data_log['y']) #fitting
    if math.sqrt(r_value**2) < 0.8: #for bad baseline fitting
        print('\nR Value is less than 0.8. It means that something went wrong. Do you want to continue? [YES / NO]: \n')
        question = input('\nYES / NO  ').upper()
        if question == 'NO':
            print('Program exited!')
            exit()
        elif question == 'YES':
            print('continuing...')
        else:
            print('write yes or no, nothing else, mate')
    data = np.genfromtxt(file, delimiter=',', skip_header=2, names=['x', 'y'])
    num_lines = sum(1 for line in open(file))
    row1 = data[0]
    value1 = row1[0]
    row_last = data[num_lines-4]
    value_last = row_last[0]
    baseline = np.arange(value_last,value1,0.5) #generate baseline x points
    #baseline = np.arange(220, 350, 0.5)
    with open('baseline.csv', 'w', newline='') as f_baseline:
        for i, number in enumerate(baseline):
            print(str(math.log10(number))+','+str(math.log10(number)*slope+intercept),file=f_baseline) #generates whole range baseline
    for point in data:
        if float(point[0]) == 280: #checks the proper 280 abs value
            global  abs280
            abs280 = point[1]
    data_baseline = np.genfromtxt('baseline.csv', delimiter=',', names=['x', 'y'])
    for basepoint in data_baseline:
        #print(round(float(10**basepoint[0]),2))
        if round(float(10**basepoint[0]),2) == 280.0:
            base280 = 10 ** float(basepoint[1])
            global corrected_abs
            corrected_abs = (abs280 - base280)
        if corrected_abs < 0:
            print('Protein concentration is negative, check the baseline!')
            global skip
            skip = 1
            break
    #print(corrected_abs)

def cuvette_extintion(correctedabs):
    if skip != 1:
        while True:
            try:
                cuvette = float(input('\nWhat was the pathlenght of the cuvette used for the measurement? [in cm] :   '))
                extintion = int(input('\nWhat is the extinction coefficient of the measured protein? :   '))
                global protein_conc
                protein_conc = correctedabs / cuvette / extintion *1000000
                print('Final protein concentration is:  ',round(protein_conc,2),' uM')
                break
            except ValueError:
                print('\ntry again, mate...')

def clean(file):
    if os.path.isfile('log_fit.tmp'):
        os.remove('log_fit.tmp')
        os.remove('baseline.csv')
        os.remove(str(file).replace('.csv','')+'_log.csv')

def file_list(): #use for many csv files in the directory
    pliki_csv = []
    lista = os.listdir('.')
    for i,file in enumerate(lista):
        if file.endswith('.csv'):
            pliki_csv += [file]
    if pliki_csv == []:
        print('No *.csv files found. Exiting...')
        exit()
    for j,file in enumerate(pliki_csv):
            print(j,'\t',file, )
    global file_name
    while True:
        try:
            choise = int(input('\n\nProvide the file name (number):\t',))
            break
        except ValueError:
            print('Try again, mate...')
    file_name = pliki_csv[choise]

def file_check(file):
    with open(file, 'r', newline='') as f:
        for line in csv.reader(f):
            line_num = 0
            for line1 in line:
                line_num+=1
                if line_num > 3:
                    print('\nYour file has more than two columns. I cannot hadle it yet, mate.\t',line)
                    exit()
            if line_num < 2 and line != []:
                print('\nYour file has less than two columns. I cannot hadle it yet, mate.\t',line)
                exit()
print('\n\n'
      '============================================================\n'
      'This script calculates linear regression in log-log plot \n'
      'and corrects Absorbance for scattering. You can use it with \n'
      'exported *.csv files from Cary \n'
      '(change it in Setup of Spectrum measurement)\n'
      'To run script you need: numpy, matplotlib and scipy is required\n'
      '============================================================\n\n')
sleep(1)
file_list()
file_check(file_name)
log_log(file_name)
fitting(file_name)
cuvette_extintion(corrected_abs)
generate_plot(file_name)
clean(file_name)
