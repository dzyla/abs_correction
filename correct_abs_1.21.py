'''script written by Dawid Zyla for automatic scattering correction of Absorbance spectrum.
work with files directly exported from Varian Cary as *.csv. Applies linear regression in
log-log plot between log(350) and log(310).
To run it you need numpy, matplotlib and scipy.
Use is inside the directory with csv files.
example of use in powershell with python 3.5:
C:\Python35\python.exe .\correct_abs.py
v1.1 - files with less or more than 3 columns (x,y and '') will not be taken for calculations (added file_check)
v1.2 - added multiple column splitting,
v.1.21 - revision and cleaning of the code and many changes, changed to 3 digits of protein conc and plot export automaticaly
         low R^2 alert decreased to 0.5
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
splitfile = False


def file_list(): #use for many csv files in the directory
    global file_name
    pliki_csv = []
    file_numbers = []
    lista = os.listdir('.')
    for i,file in enumerate(lista):
        if file.endswith('.csv'):
            pliki_csv += [file]
    if not pliki_csv:
        print('No *.csv files found. Exiting...')
        exit()
    for j,file in enumerate(pliki_csv):
            print(j,'\t',file, )
            file_numbers += [j]
    while True:
        try:
            choice = int(input('\n\nProvide the file name (number):\t', ))
            if choice in file_numbers:
                file_name = pliki_csv[choice]
                break
        except ValueError:
            print('Try again, mate...')


def file_check(file):
    with open(file, 'r', newline='') as f:
        num_lines = sum(1 for line in open(file)) -1
        for i, row in enumerate(csv.reader(f)):
            if i == int(num_lines / 2): #in the middle of data check number of columns
                row_size = int((len(row)))  # number of columns in half of dataset
                if '' in row:
                    row_size -= 1
                if row_size > 2 : #check first 5 rows for number of columns
                    print(file,'\tFile has more than 2 columns, splitting file to new files')
                    global splitfile
                    splitfile = True
                    break

def split_file(file):
    if splitfile == 1:
        headers = []
        with open(file) as csvfile:  # check for headers with error when try to float the name
            reader = csv.reader(csvfile, delimiter=',')
            for i, row in enumerate(reader):
                for element in row:
                    header = True
                    if element != '': #the last element tends to be ''
                        while header:
                            try:
                                float(element) #checks if headers are names
                                break
                            except ValueError:
                                if i not in headers: #if error then add line to headers list
                                    headers += [i]
                                header = False
        for fn in os.listdir('.'):
            if os.path.isfile(fn) and "_col-" in fn:
                os.remove(fn)
        with open(file) as csvfile:  # check for headers with error when try to float the name
            reader = csv.reader(csvfile, delimiter=',')
            for i, row in enumerate(reader):
                if i not in headers and row != []: #for all rows which are not in headers list
                    rows = len(row)
                    if row[rows - 1] == '':  # remove the last empty row from row number
                        rows -= 1
                    for r in range(rows):
                        if r == 0 or r % 2 == 0:  # every second column is x 0,2 etc..
                            newrow = (row[rows - rows + r],
                                      row[rows - rows + r + 1])  # newrow gets every two colums together, starting from [0]
                            file_name = open(file.replace('.csv', '_col-') + str(r) + '.csv', 'a', newline='')
                            with file_name:
                                writer = csv.writer(file_name, delimiter=',')
                                writer.writerow(newrow)

def log_log(file): #generation of log-log plot and 350-310 nm for baseline correction
    if not splitfile:
        data = np.genfromtxt(file, delimiter=',', skip_header=2, names=['x', 'y'])
        with open(str(file).replace('.csv','')+'_log.csv', 'w', newline='') as f_out: #whole dataset
            with open(file.replace('.csv','')+'_log_fit.tmp', 'w', newline='') as flog_out: #only 350-310 nm for baseline fitting
                for i, line in enumerate(data):
                    x_log = math.log10(float(line[0]))
                    y_log = math.log10(float(math.sqrt(line[1]**2)))
                    data = (str(x_log)+','+str(y_log))
                    if 350 >= float(line[0]) >= 310:
                        print(data, file=flog_out)
                    print(data,file=f_out)
    if splitfile:
        print('\n')
        #creates new log-log for all datasets
        for fn in os.listdir('.'):
            if os.path.isfile(fn) and "_col-" in fn:
                print('generating  ',fn)
                data = np.genfromtxt(fn, delimiter=',', names=['x', 'y'])
                with open(str(fn).replace('.csv', '') + '_log.csv', 'w', newline='') as f_out:
                    with open(fn.replace('.csv','')+'_log_fit.tmp', 'w', newline='') as flog_out:
                        for i, line in enumerate(data):
                            x_log = math.log10(float(line[0]))
                            y_log = math.log10(float(math.sqrt(line[1] ** 2)))
                            data = (str(x_log) + ',' + str(y_log))
                            if 350 >= float(line[0]) >= 310:
                                print(data, file=flog_out)
                            print(data, file=f_out)

def fitting(file):
    #global declaration
    global r_value
    global skip
    global abs280
    global corrected_abs
    global file_name
    '''---------------------------------------------------------------------------------------------'''
    if not splitfile:
        data_log = np.genfromtxt(file.replace('.csv','')+'_log_fit.tmp', delimiter=',', skip_header=2, names=['x', 'y']) #np array for the baseline fitting
        slope, intercept, r_value, p_value, std_err = stats.linregress(data_log['x'], data_log['y']) #real fitting part
        if math.sqrt(r_value**2) < 0.8: #checking Rsquared if less than 80%, which would mean that fit is poor
            print('\nR Value is less than 0.5. It means that something went wrong. Do you want to continue? [YES / NO]: \n')
            question = input('\nYES / NO  ').upper()
            if question == 'NO':
                print('Program exited!')
                exit()
            elif question == 'YES':
                print('continuing...')
            else:
                print('write yes or no, nothing else, mate')
        data = np.genfromtxt(file, delimiter=',', skip_header=2, names=['x', 'y']) #taking the raw file to generate x values for full axis baseline
        num_lines = sum(1 for line in open(file))
        row1 = data[0]
        value1 = row1[0]
        row_last = data[num_lines-4]
        value_last = row_last[0]
        baseline = np.arange(value_last,value1,0.5) #generate baseline x points
        #baseline = np.arange(220, 350, 0.5) for fixed baseline creation
        with open(file.replace('.csv','')+'_baseline.csv', 'w', newline='') as f_baseline: #creative the *_baseline.csv file
            for i, number in enumerate(baseline):
                print(str(math.log10(number))+','+str(math.log10(number)*slope+intercept),file=f_baseline) #generates whole range baseline
        for point in data:
            if float(point[0]) == 280: #checks the proper 280 abs value
                abs280 = point[1]
        data_baseline = np.genfromtxt(file.replace('.csv','')+'_baseline.csv', delimiter=',', names=['x', 'y'])
        for basepoint in data_baseline:
            if round(float(10**basepoint[0]),2) == 280.0:
                base280 = 10 ** float(basepoint[1])
                corrected_abs = (abs280 - base280)
            if corrected_abs < 0:
                print('Protein concentration is negative, check the baseline!')
                skip = 1
                break
    if splitfile:
        '''-------------------File choosing of splited file---------------------------------------------------'''
        print('\nYour file was divided to separate files for each measurement. Please chose column:')
        pliki_csv = []
        lista = os.listdir('.')
        print('\n')
        for i, filen in enumerate(lista):
            if '_col-' in filen and not '_log' in filen:
                pliki_csv += [filen]
        for j, filen in enumerate(pliki_csv):
            print(j, '\t', filen, )
        while True:
            try:
                choice = int(input('\n\nProvide the file name (number):\t', ))
                break
            except ValueError:
                print('Try again, mate...')
        file_name = pliki_csv[choice]
        print(file_name)
        '''--------------------------------Fitting part------------------------------------------------------'''
        data_log = np.genfromtxt(file_name.replace('.csv','')+'_log_fit.tmp', delimiter=',', names=['x', 'y'])  # for baseline
        slope, intercept, r_value, p_value, std_err = stats.linregress(data_log['x'], data_log['y'])  # fitting
        if math.sqrt(r_value ** 2) < 0.5:  # for bad baseline fitting
            print(
                '\nR Value is less than 0.5. It means that something went wrong. Do you want to continue? [YES / NO]: \n')
            question = input('\nYES / NO  ').upper()
            if question == 'NO':
                print('Program exited!')
                exit()
            elif question == 'YES':
                print('\ncontinuing...')
            else:
                print('write yes or no, nothing else, mate')

        '''-------------------------Baseline creation and calculating corrected ABS--------------------------'''
        data = np.genfromtxt(file_name, delimiter=',', names=['x', 'y']) #open the chosen file and find first and last x value
        num_lines = sum(1 for line in open(file_name))
        row1 = data[0]
        value1 = row1[0]
        row_last = data[num_lines-1]
        value_last = row_last[0]
        baseline = np.arange(value_last, value1, 0.5)  # generate baseline with all x points every spacing, here 0.5
        # baseline = np.arange(220, 350, 0.5) #alternative baseline creation
        with open(str(file_name).replace('.csv','')+'_baseline.csv', 'w', newline='') as f_baseline:
            for i, number in enumerate(baseline):
                print(str(math.log10(number)) + ',' + str(math.log10(number) * slope + intercept),
                      file=f_baseline)  # generates whole range baseline according to linear fit done before
        for point in data:
            if float(point[0]) == 280:  # checks measured 280 abs value
                abs280 = point[1]
        data_baseline = np.genfromtxt(str(file_name).replace('.csv','')+'_baseline.csv', delimiter=',', names=['x', 'y'])
        for basepoint in data_baseline:
            # print(round(float(10**basepoint[0]),2))
            if round(float(10 ** basepoint[0]), 2) == 280.0:
                base280 = 10 ** float(basepoint[1])
                corrected_abs = (abs280 - base280)
            if corrected_abs < 0:
                print('Protein concentration is negative, check the baseline!')
                skip = 1
                break

def cuvette_extintion(correctedabs):
    global protein_conc
    if skip != 1:
        while True:
            try:
                cuvette = float(input('\nWhat was the pathlenght of the cuvette used for the measurement? [in cm] :   '))
                extintion = int(input('\nWhat is the extinction coefficient of the measured protein? :   '))
                protein_conc = correctedabs / cuvette / extintion *1000000
                print('\nFinal protein concentration is:  ',round(protein_conc,3),' uM') #changed to 3 digits
                break
            except ValueError:
                print('\ntry again, mate...')

def generate_plot(file_name):
    if not splitfile:
        data = np.genfromtxt(str(file_name).replace('.csv','')+'_log.csv', delimiter=',', skip_header=2, names=['x', 'y'])
    if splitfile:
        data = np.genfromtxt(str(file_name).replace('.csv', '') + '_log.csv', delimiter=',', names=['x', 'y'])

    # print(data)
    '''-----------------Graph generation----------------------------------------------------------------------'''
    data_baseline = np.genfromtxt(file_name.replace('.csv','')+'_baseline.csv', delimiter=',', names=['x', 'y'])
    pyplot.xlabel('log(Wavelength)')
    pyplot.ylabel('log(Absorbance)')
    pyplot.plot(data['x'], data['y'], label='log(Absorbance)')
    pyplot.plot(data_baseline['x'], data_baseline['y'], label='log(Baseline)')
    pyplot.legend(loc=1)  # upper right corner, anticlockwise
    pyplot.suptitle(file_name, size=15)  # titlesize and name
    pyplot.annotate(
        'Protein concentration:  ' + str(round(protein_conc, 2)) + ' uM\n$R^{2}$ = ' + str(round(-r_value, 3)),
        xy=(math.log10(280), math.log10(math.sqrt(abs280**2))),
        xytext=(math.log10(280) - 0.05, math.log10(float(abs280) + 0.2)),
        arrowprops=dict(facecolor='black', shrink=0.05))
    # axes = pyplot.gca() #for axes adjustment
    # axes.set_ylim(ymax=0) #set the y range
    pyplot.savefig(file_name.replace('.csv', '' + '.png'), dpi=200) #exports file automatically :)
    pyplot.show()

def clean():
    pliki_csv =[]
    lista = os.listdir('.')
    for i,file in enumerate(lista):
        if '_col-' in file and not 'png' in file:
            pliki_csv += [file]
    for i in pliki_csv:
        os.remove(i)
    pliki_csv = []
    lista = os.listdir('.')
    for i,file in enumerate(lista):
        if '.tmp' in file or '_baseline' in file or '_log' in file:
            pliki_csv += [file]
    for i in pliki_csv:
        os.remove(i)


'''==============='''

print('\n\n'
      '============================================================\n'
      'This script calculates linear regression in log-log plot \n'
      'and corrects Absorbance for scattering. You can use it with \n'
      'exported *.csv files from Cary or any other files win many columns\n'
      'To run script you need: numpy, matplotlib and scipy is required\n'
      '============================================================\n\n')
sleep(1)


file_list()
print(file_name)
file_check(file_name)
split_file(file_name)
log_log(file_name)
fitting(file_name)
cuvette_extintion(corrected_abs)
generate_plot(file_name)
clean()