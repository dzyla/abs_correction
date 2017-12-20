'''script written by Dawid Zyla for automatic scattering correction of Absorbance spectra.
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
v.1.22 - added check for A260/A280 ratio for DNA contaminated samples. To use this pure protein A260/A280 ratio is needed.
v1.3 - code was rewrited, some adjustments, added proper DNA contamination part, now all files are 'splited' and also the
        column name is printed in file name. Code should be more robust to different csv files
'''

def file_list():  # use for many csv files in the directory, return filename
    import os
    pliki_csv = []
    file_numbers = []
    lista = os.listdir('.')
    for i, file in enumerate(lista):
        if file.endswith('.csv'):
            pliki_csv += [file]
    if not pliki_csv:
        print('No *.csv files found. Exiting...')
        exit()
    for j, file in enumerate(pliki_csv):
        print(j, '\t', file, )
        file_numbers += [j]
    while True:
        try:
            choice = int(input('\n\nProvide the file name (number):\t', ))
            if choice in file_numbers:
                file_name = pliki_csv[choice]
                return file_name
        except ValueError:
            print('Try again, mate...')

def split_file(file):  #now splitfile is not required because it works for all kind of files
    import csv
    import os
    headers = []
    splited_files =[]
    with open(file) as csvfile:  # check for headers with error when try to float the name
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            for element in row:
                header = True
                if element != '':  # the last element tends to be ''
                    while header:
                        try:
                            float(element)  # checks if headers are names
                            break
                        except ValueError:
                            if i not in headers:  # if error then add line to headers list
                                headers += [i]
                                if i == 0:
                                    sample_name = row
                            header = False
    for fn in os.listdir('.'):  #remove old column files
        if os.path.isfile(fn) and "_col-" in fn:
            os.remove(fn)
    with open(file) as csvfile:  # check for headers with error when try to float the name
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            if i not in headers and row != []:  # for all rows which are not in headers list
                rows = len(row)
                if row[rows - 1] == '':  # remove the last empty row from row number
                    rows -= 1
                for r in range(rows):
                    if r == 0 or r % 2 == 0:  # every second column is x 0,2 etc..
                        newrow = (row[rows - rows + r],
                                  row[
                                      rows - rows + r + 1])  # newrow gets every two colums together, starting from [0]
                        namerow = str(sample_name[rows - rows + r]) #gets the name from sample names
                        file_name = open(file.replace('.csv', '_col-'+str(namerow)+'.csv'), 'a', newline='')
                        file_name_str = str(file.replace('.csv', '_col-'+str(namerow)+'.csv'))
                        with file_name:
                            writer = csv.writer(file_name, delimiter=',')
                            writer.writerow(newrow)
    return file_name_str

def subselection(file):
    import os

    pliki_csv = []
    file_name = ''
    lista = os.listdir('.')
    print('\n')
    file_body = file.replace('.csv','')
    for i, filen in enumerate(lista):
        if file_body in filen and not '_log' in filen and not str(file) in filen and not '.png' in filen:
            pliki_csv += [filen]
    if len(pliki_csv) > 1:
        print('\nYour file was divided to separate files for each measurement. Please choose column:')
        for j, filen in enumerate(pliki_csv):
            print(j, '\t', str(filen)[len(file)+1:-4:])
        while True:
            try:
                choice = int(input('\n\nProvide the file name (number):\t', ))
                file_name = pliki_csv[choice]
                return file_name
            except ValueError:
                print('Try again, mate...')
    elif len(pliki_csv) == 1:
        file_name = pliki_csv[0]
        return file_name
    else:
        print('No files?')
        exit(1)

def log_log(file_from_subselection):
    import numpy as np
    import math

    '''==========================creates new log-log for all datasets======================================='''

    print('generating  log-log data for ',file_from_subselection)
    data = np.genfromtxt(file_from_subselection, delimiter=',', names=['x', 'y'])  #loads data from text using numpy
    with open(str(file_from_subselection).replace('.csv', '') + '_log.csv', 'w', newline='') as f_out:  #for the loglog
        loglog_file = str(file_from_subselection).replace('.csv', '') + '_log.csv'
        with open(str(file_from_subselection).replace('.csv', '') + '_log_fit.tmp', 'w', newline='') as flog_out:  #for baseline
            loglog_fit = str(file_from_subselection).replace('.csv', '') + '_log_fit.tmp'
            for i, line in enumerate(data):
                x_log = math.log10(float(line[0]))
                y_log = math.log10(float(abs(line[1])))  #gets the absolute value of y
                data = (str(x_log) + ',' + str(y_log))
                if 350 >= float(line[0]) >= 310:
                    print(data, file=flog_out)
                print(data, file=f_out)
    return loglog_file,loglog_fit

def fitting(file,logfile,log_fit):  #file - converted file by split_file, logfile - file in loglog, log_fit - template for baseline
    import numpy as np
    from scipy.stats import stats
    import math
    skip = False

    '''--------------------------------Fitting part------------------------------------------------------'''
    data_log = np.genfromtxt(log_fit, delimiter=',',
                             names=['x', 'y'])  # data generation for baseline from baseline tmp file

    slope, intercept, r_value, p_value, std_err = stats.linregress(data_log['x'], data_log['y'])  # fitting

    if abs(r_value) < 0.5:  # for bad baseline fitting
        print(
            '\nR Value is less than 0.5. It means that something went wrong. Do you want to continue? [YES / NO]: \n')
        question = input('\nYES / NO  ').upper()
        if question == 'NO':
            print('Program exited!')
            exit()
        elif question == 'YES' or question == '':
            print('\ncontinuing...')
        else:
            print('write yes or no, nothing else, mate')
    data = np.genfromtxt(file, delimiter=',',
                             names=['x', 'y'])  # taking the raw file to generate x values for full axis baseline
    #num_lines = sum(1 for line in open(logfile))  #check how many lines are in the input file

    '''--------------------------------------baseline generation-----------------------------------------------------'''

    vmax = 350 #default settings for vmax and min
    vmin = 220
    row1 = data[0] #choose first tuple
    value1 = row1[0]  #choose first value of the tuple, which is wavelength
    row2 = data[1]
    value2 = row2[0]
    row_last = data[-1]  #last row of the file
    value_last = row_last[0]  #last wavelength
    if value1 > value_last:  #check which value is bigger
        vmax = value1
        vmin = value_last
    elif value_last > value1:
        vmax = value_last
        vmin = value1
    spacing = abs(value1-value2) #spacing between two following points
    baseline = np.arange(vmin, vmax, spacing)  # generate baseline x points with spacing
    baseline_file = str(logfile.replace('.csv', '') + '_baseline.csv')
    with open(baseline_file, 'w',
              newline='') as f_baseline:  # creative the *_baseline.csv file
        for i, number in enumerate(baseline):
            print(str(math.log10(number)) + ',' + str(math.log10(number) * slope + intercept),
                  file=f_baseline)  # generates whole range baseline

    '''-----------------------------------checking for A280 if not negative----------------------------------------'''
    for point in data:  #data in raw file
        if float(point[0]) == 280:  # checks the measured 280 abs value
            abs280 = point[1]
    data_baseline = np.genfromtxt(baseline_file, delimiter=',', names=['x', 'y'])  #load created baseline
    corrected_abs = 0
    for basepoint in data_baseline:
        if round(float(10 ** basepoint[0]), 2) == 280.0:
            base280 = 10 ** float(basepoint[1])
            corrected_abs = (abs280 - base280)
            if corrected_abs < 0:
                print('Protein concentration is negative, check the baseline!')
                skip = True
                return baseline_file,skip,0,r_value,abs280
    return baseline_file,skip,corrected_abs, r_value,abs280

def dna_contamination_check(logfile,baseline_file):
    from statistics import stdev
    import numpy as np
    dna_cont = False

    '''----------------checking corrected abs from log-log plot----------------------'''

    data = np.genfromtxt(logfile, delimiter=',', names=['x', 'y'])
    for point in data:
        if round(10 ** (float(point[0])), 2) == 280.0:  # checks measured 280 abs value
            abs280 = 10 ** (float(point[1]))
        if round(10 ** (float(point[0])), 2) == 260.0:  # checks measured 260 abs value
            abs260 = 10 ** (float(point[1]))
    data_baseline = np.genfromtxt(baseline_file, delimiter=',', names=['x', 'y'])
    for basepoint in data_baseline:
        if round(float(10 ** basepoint[0]), 2) == 280.0:
            base280 = 10 ** float(basepoint[1])
        if round(float(10 ** basepoint[0]), 2) == 260.0:
            base260 = 10 ** float(basepoint[1])
    a260 = abs260 - base260
    a280 = abs280 - base280


    '''-----------------------calculation of DNA contamination--------------------------'''
    protein_ratio = 0.67734008  # ratio of your protein A260 to A280 - you need pure protein spectrum to check this!
    dna_ratio = 2  # estimated DNA A260 / A280 ratio. If other value needed, change it here
    ratio_measured = a260 / a280
    proc_protein = 100  #if there is no DNA protein is 100% of the sample

    saved_standev = 1.0
    if a260 > a280:
        dna_cont = True
        calc_dna = str(input(('\nThere is possible DNA contamination in the sample. Do you want to calc DNA contamination? (default - NO)')))
        if calc_dna.upper() == 'YES':
            for i in range(0, 5000, 1):  # itarates through 0.001 till 5.0 as a fraction (from 0.1% till 83%)
                fraction = i / 1000  # divides number i to generate fine slicing of search data (fraction of protein)
                average = (dna_ratio * 1 + (protein_ratio * fraction)) / (
                    1 + fraction)  # calculates the avarage with given fraction
                standdev = stdev([ratio_measured, average])  # checks how close are those two values
                if standdev < saved_standev:  # checks the lowest stdev
                    saved_standev = standdev  # sets new lowest
                    proc_protein = fraction / (
                        fraction + 1) * 100  # calculates the fraction of protein when DNA fraction is 1

            print(round(proc_protein, 2), '% of the sample is protein.\nA260 / A280 ratio is: ', round(ratio_measured, 2))
        elif calc_dna.upper() == 'NO' or calc_dna.upper() == '':
            print('Skipping...')
    return dna_cont,proc_protein

def cuvette_extintion(correctedabs,skip):
    calc_mass = False
    if skip != 1:
        while True:
            try:
                cuvette = float(
                    input('\nWhat was the pathlenght of the cuvette used for the measurement? [in cm] :   '))
                extintion = int(input('\nWhat is the extinction coefficient of the measured protein? :   '))
                mass = input('\nWhat is the size in kDa of your protein? (enter to skip)   ')
                if mass != '':
                    mass = float(mass)
                    calc_mass = True
                protein_conc = correctedabs / cuvette / extintion * 1000000
                if calc_mass:
                    mgml = protein_conc*mass/1000
                print('\nFinal protein concentration is:  ', round(protein_conc, 3), ' uM')  # changed to 3 digits
                if calc_mass:
                    print('\nFinal protein concentration is:  ', round(mgml, 2), ' mg/ml')  # mg/ml
                    return protein_conc,mgml
                if not calc_mass:
                    return protein_conc,0
                break
            except ValueError:
                print('\ntry again, mate...')
    if skip == 1:
        return 0, 0

def generate_plot(file,log_file,baseline_file,protein_conc,r_value,abs280,mgml):
    import matplotlib.pyplot as pyplot
    import numpy as np
    import math


    '''-----------------Graph generation----------------------------------------------------------------------'''

    data = np.genfromtxt(log_file, delimiter=',', names=['x', 'y'])
    data_baseline = np.genfromtxt(baseline_file, delimiter=',', names=['x', 'y'])
    pyplot.xlabel('log(Wavelength)')
    pyplot.ylabel('log(Absorbance)')
    pyplot.plot(data['x'], data['y'], label='log(Absorbance)')
    pyplot.plot(data_baseline['x'], data_baseline['y'], label='log(Baseline)')
    pyplot.legend(loc=1)  # upper right corner, anticlockwise
    pyplot.suptitle(log_file.replace('.csv',''), size=11)  # titlesize and name
    pyplot.rcParams["figure.figsize"] = (20,15)
    if mgml != 0:
        pyplot.annotate(
                'Protein concentration:  ' + str(round(protein_conc, 2)) + ' uM\n                                     '+str(round(mgml,2))+' mg/ml\n$R^{2}$ = ' + str(round(-r_value, 3)),
        xy=(math.log10(280), math.log10(abs(abs280))),
        xytext=(0.5, 0.7), textcoords='axes fraction',
        arrowprops=dict(arrowstyle="-|>",facecolor='black'))
    if mgml == 0:
        pyplot.annotate(
                'Protein concentration:  ' + str(round(protein_conc, 2)) + ' uM\n$R^{2}$ = ' + str(round(-r_value, 3)),
        xy=(math.log10(280), math.log10(abs(abs280))),
        xytext=(0.5, 0.7), textcoords='axes fraction',
        arrowprops=dict(arrowstyle="-|>", facecolor='black'))


    pyplot.savefig(file.replace('.csv', '' + '.png'), dpi=200)  # exports file automatically :)
    pyplot.show()

def clean():
    import os
    pliki_csv = []
    lista = os.listdir('.')
    for i, file in enumerate(lista):
        if '_col-' in file and not 'png' in file:
            pliki_csv += [file]
    for i in pliki_csv:
        os.remove(i)
    pliki_csv = []
    lista = os.listdir('.')
    for i, file in enumerate(lista):
        if '.tmp' in file or '_baseline' in file or '_log' in file:
            pliki_csv += [file]
    for i in pliki_csv:
        os.remove(i)

def check_modules():
    try:
        import scipy
        import numpy
        import matplotlib
        import math
        import os
    except ImportError:
        print('Required modules not found (scipy, numpy, matplotlib). Install them to run the script. Try pip3 install [modulename]')

def measure_again(file):
    import os

    pliki_csv = []
    file_name = ''
    lista = os.listdir('.')
    file_body = file.replace('.csv', '')
    for i, filen in enumerate(lista):
        if file_body in filen and not '_log' in filen and not str(file) in filen and not '.png' in filen:
            pliki_csv += [filen]
    if len(pliki_csv) > 1:
        question = str(input('\nThere are more column in your file. Do you want to check other column? (YES / NO), or press ENTER to exit ').upper())
        if question == '' or question == 'NO':
            return False
        elif question == 'YES':
            return True


'''========================================================================='''

print('\n\n'
      '============================================================\n'
      'This script calculates linear regression in log-log plot \n'
      'and corrects Absorbance for scattering. You can use it with \n'
      'exported *.csv files from Cary or any other files with many columns\n'
      'To run the script, numpy, matplotlib and scipy is required\n'
      '============================================================\n\n')


check_modules()
file = file_list()
converted_file = split_file(file)
again = True
while again:
    selected_file = subselection(file)
    log_file, log_fit = log_log(selected_file)
    baseline_file,skip,corrected_abs,r_value,abs280 = fitting(selected_file,log_file,log_fit)
    dna_cont,proc_prot = dna_contamination_check(log_file,baseline_file)
    protein_conc, mgml = cuvette_extintion(corrected_abs,skip)
    generate_plot(selected_file,log_file,baseline_file,protein_conc,r_value,abs280,mgml)
    again = measure_again(file)
clean()

