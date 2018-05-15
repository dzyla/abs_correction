'''script written by Dawid Zyla for automatic scattering correction of Absorbance spectra.
work with files directly exported from Varian Cary as *.csv. Applies linear regression in
log-log plot between log(350) and log(310).
To run it you need numpy, matplotlib and scipy.
Use is inside the directory with csv files.
example of use in powershell with python 3.5:
C:\Python35\python.exe .\correct_abs.py
v1.1 - files with less or more than 3 columns (x,y and '') will not be taken for calculations (added file_check)
v1.2 - added multiple column splitting,
v.1.21 - revision and cleaning of the code and many changes, changed to 3 digits of protein conc and plot export automatically
         low R^2 alert decreased to 0.5
v.1.22 - added check for A260/A280 ratio for DNA contaminated samples. To use this pure protein A260/A280 ratio is needed.
v1.3 - code was rewritten, some adjustments, added proper DNA contamination part, now all files are 'splited' and also the
        column name is printed in file name. Code should be more robust to different csv files
v1.41 - basically everything was written again, no more temporary files, implemented gui elements, fitting the baseline to
        wanted range, saving settings for known protein and reading them from file, improved speed of loading file.
'''

import os
import sys

import easygui
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import stats

'''general use functions'''


def is_number(s):
    try:
        float(s)  # for int, long and float
    except ValueError:
        try:
            complex(s)  # for complex
        except ValueError:
            return False
        except TypeError:
            return False
    except TypeError:
        return False
    return True


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


'''-----------------------------------------'''


def choose_file():
    csv_file = easygui.fileopenbox()
    # print(csv_file)
    if csv_file:
        return csv_file
    else:
        print('No file provided')
        quit()


def open_file(file):
    import pandas as pd
    data = pd.read_csv(file, header=None)
    # get the names from file

    name_list = np.array(data[:2][0])
    data = data.dropna(axis=1, how='all')

    headers = []
    header = False

    for row, element in enumerate(name_list):
        if not is_number(element):
            headers.append(row)
            header = True

    if header:
        sample_row = data.iloc[[0]].as_matrix()
        sample_row = sample_row[~pd.isnull(sample_row)]

        data_names = []
        samples = []
        for sample_name in sample_row:
            data_names.append('Wavelength')
            data_names.append(sample_name)
            samples.append(sample_name)
        data = data[2:]
        # change data to np.array
        data = pd.DataFrame.as_matrix(data)

        print('Opening: ' + file + '\n')
        print('Data columns found in file: ' + ''.join(samples) + '\n')
        return data_names, data, samples, file

    elif not header:
        print('Opening: ' + file + '\n')
        print('No data description found in ' + file + '\n')
        # remove empty columns

        # change data to np.array
        data = pd.DataFrame.as_matrix(data)
        data_names = 0
        samples = None
        return data_names, data, samples, file


def choose_column(data, data_names, samples):
    msg = "Multiple data columns found. Which one should I analyze?"
    title = "Choose the data column"
    choices = samples
    choice = easygui.choicebox(msg, title, choices)
    if choice != None:
        data_column = data_names.index(choice)
        return data[:, data_column - 1:data_column + 1]
    else:
        print('No data file chosen!')
        quit()

def get_x_y(data):
    # change the data to float arrays
    x = np.asfarray(data[:, 0], float)
    y = np.asfarray(data[:, 1], float)
    return x, y

def log_log_plot(data, settings):
    # change the data to float arrays
    x = np.asfarray(data[:, 0], float)
    y = np.asfarray(data[:, 1], float)

    # log-log  - abs to remove negative values
    log_x = np.log10(np.abs(x))
    log_y = np.log10(np.abs(y))
    if not settings:
        msg = "Fit the baseline between 310 / 350 nm?"
        choices = ["Yes", "No"]
        reply = easygui.buttonbox(msg, choices=choices)

        # here might be added a window to change the range
        if reply == 'Yes':
            fitting_region_low = int(np.where(x == find_nearest(x, 310))[0])
            fitting_region_high = int(np.where(x == find_nearest(x, 350))[0])

            if fitting_region_low < fitting_region_high:
                temp = fitting_region_low
                fitting_region_low = fitting_region_high
                fitting_region_high = temp

            Wavelength = x[np.where(x == find_nearest(x, 280))[0]]

        elif reply == 'No':
            msg = "Enter the fitting region for the baseline correction"
            title = "Baseline fitting range"
            fieldNames = ["Lower boundary (310):", "Higher boundary (350):", "Peak Wavelength"]
            fieldValues = []  # we start with blanks for the values
            fieldValues = easygui.multenterbox(msg, title, fieldNames)

            while 1:
                errmsg = ""
                for i in range(len(fieldNames)):
                    if fieldValues[i].strip() == "":
                        if i == "":
                            if i != 2:
                                errmsg = errmsg + ('"%s" is a required field.\n\n' % fieldNames[i])
                    if fieldValues[2].strip() == "":
                        Wavelength = None
                if float(fieldValues[0]) < np.min(x):
                    errmsg = errmsg + "lower boundary is not present in dataset!"
                elif float(fieldValues[1]) > np.max(x):
                    errmsg = errmsg + "higher boundary is not present in dataset!"
                if errmsg == "": break  # no problems found
                fieldValues = easygui.multenterbox(errmsg, title, fieldNames, fieldValues)

            x1 = float(fieldValues[1])
            x2 = float(fieldValues[0])
            fitting_region_low = int(np.where(x == find_nearest(x, x1))[0])
            fitting_region_high = int(np.where(x == find_nearest(x, x2))[0])
            Wavelength = float(fieldValues[2])
            Wavelength = x[np.where(x == find_nearest(x, Wavelength))[0]]
            if fitting_region_low < fitting_region_high:
                temp = fitting_region_low
                fitting_region_low = fitting_region_high
                fitting_region_high = temp

        fitting_region_x = log_x[fitting_region_high:fitting_region_low]
        fitting_region_y = log_y[fitting_region_high:fitting_region_low]
        return log_x, log_y, fitting_region_x, fitting_region_y, x, y, Wavelength

    elif settings:
        fitting_region_low = int(np.where(x == find_nearest(x, 310))[0])
        fitting_region_high = int(np.where(x == find_nearest(x, 350))[0])

        if fitting_region_low < fitting_region_high:
            temp = fitting_region_low
            fitting_region_low = fitting_region_high
            fitting_region_high = temp

        Wavelength = x[np.where(x == find_nearest(x, 280))[0]]
        fitting_region_x = log_x[fitting_region_high:fitting_region_low]
        fitting_region_y = log_y[fitting_region_high:fitting_region_low]
        return log_x, log_y, fitting_region_x, fitting_region_y, x, y, Wavelength


def baseline_correction(log_x, log_y, fitting_region_x, fitting_region_y, x_data, y_data, wavelength):
    slope, intercept, r_value, p_value, std_err = stats.linregress(fitting_region_x, fitting_region_y)
    baseline_y = log_x * slope + intercept

    if wavelength == 280:
        # get the value of A280 from raw file, from baseline and subtract them
        abs280_raw = y_data[int((np.where(x_data == find_nearest(x_data, 280))[0]))]
        log_abs280_baseline = baseline_y[int((np.where(x_data == find_nearest(x_data, 280))[0]))]
        abs280_baseline = 10 ** log_abs280_baseline
        abs280_corrected = abs280_raw - abs280_baseline

        skip_protein_conc = False

        # check for negative A280 difference
        if abs280_corrected < 0:
            print('Protein concentration is negative, check the baseline!')
            skip_protein_conc = True

        # print(abs280_corrected)
        return baseline_y, abs280_corrected, r_value, skip_protein_conc
    elif wavelength != 280:
        wavelength = float(wavelength)
        # get the value of selected wavelength
        abs_raw = y_data[int((np.where(x_data == find_nearest(x_data, wavelength))[0]))]
        log_abs_baseline = baseline_y[int((np.where(x_data == find_nearest(x_data, wavelength))[0]))]
        abs_baseline = 10 ** log_abs_baseline
        abs280_corrected = abs_raw - abs_baseline

        skip_protein_conc = False

        # check for negative A280 difference
        if abs280_corrected < 0:
            print('Protein concentration is negative, check the baseline!')
            skip_protein_conc = True

        # print(abs280_corrected)
        return baseline_y, abs280_corrected, r_value, skip_protein_conc


def cuvette_extintion(correctedabs, skip):
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
                    mgml = protein_conc * mass / 1000
                print('\nFinal protein concentration is:  ', round(protein_conc, 3), ' uM')  # changed to 3 digits
                if calc_mass:
                    print('\nFinal protein concentration is:  ', round(mgml, 2), ' mg/ml')  # mg/ml
                    return protein_conc, mgml
                if not calc_mass:
                    return protein_conc, 0
                break
            except ValueError:
                print('\ntry again, mate...')
    if skip == 1:
        return 0, 0


def dna_contamination_check(log_x, log_y, baseline):
    from statistics import stdev
    import numpy as np
    dna_cont = False

    '''----------------checking corrected abs from log-log plot----------------------'''

    abs280 = log_y[int(np.where(log_x == np.log10(280))[0])]
    abs260 = log_y[int(np.where(log_x == np.log10(260))[0])]

    base280 = baseline[int(np.where(log_x == np.log10(280))[0])]
    base260 = baseline[int(np.where(log_x == np.log10(260))[0])]

    a260 = abs260 - base260
    a280 = abs280 - base280

    '''-----------------------calculation of DNA contamination--------------------------'''
    protein_ratio = 0.67734008  # ratio of your protein A260 to A280 - you need pure protein spectrum to check this!
    dna_ratio = 2  # estimated DNA A260 / A280 ratio. If other value needed, change it here
    ratio_measured = a260 / a280
    proc_protein = 100  # if there is no DNA protein is 100% of the sample

    saved_standev = 1.0
    if a260 > a280:
        dna_cont = True
        msg = "There is possible DNA contamination in your protein sample. Do you want to estimate the contamination?"
        choices = ["Yes", "No"]
        reply = easygui.buttonbox(msg, choices=choices)

        # here might be added a window to change the range
        if reply == 'Yes':
            for i in range(0, 5000, 1):  # itarates through 0.001 till 5.0 as a fraction (from 0.1% till 83%)
                fraction = i / 1000  # divides number i to generate fine slicing of search data (fraction of protein)
                average = (dna_ratio * 1 + (protein_ratio * fraction)) / (
                        1 + fraction)  # calculates the avarage with given fraction
                standdev = stdev([ratio_measured, average])  # checks how close are those two values
                if standdev < saved_standev:  # checks the lowest stdev
                    saved_standev = standdev  # sets new lowest
                    proc_protein = fraction / (
                            fraction + 1) * 100  # calculates the fraction of protein when DNA fraction is 1

            print(round(proc_protein, 2), '% of the sample is protein.\nA260 / A280 ratio is: ',
                  round(ratio_measured, 2))
        elif reply == 'Yes':
            print('Skipping...')
    return dna_cont, proc_protein


def cuvette_extintion_gui(correctedabs, skip):
    msg = "Enter the protein information"
    title = "Protein information"
    fieldNames = ["Pathlenght", "Extinction coefficient", "Protein mass in kDa"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = easygui.multenterbox(msg, title, fieldNames)

    # make sure that none of the fields was left blank

    check_mgml = True

    while 1:
        errmsg = ""
        for i in range(len(fieldNames)):
            if fieldValues[i].strip() == "":
                if i != 2:
                    errmsg = errmsg + ('"%s" is a required field.\n\n' % fieldNames[i])
                if fieldValues[2].strip() == "":
                    mgml = 0
                    check_mgml = False
        if errmsg == "": break  # no problems found
        fieldValues = easygui.multenterbox(errmsg, title, fieldNames, fieldValues)
    cuvette = fieldValues[0]
    extintion = fieldValues[1]
    protein_conc = correctedabs / float(cuvette) / float(extintion) * 1000000

    if check_mgml:
        mgml = fieldValues[2]
        protein_mgml = protein_conc * float(mgml) / 1000
    else:
        protein_mgml = 0
    return float(protein_conc), float(protein_mgml), cuvette, extintion, mgml


def calc_protein_conc_from_settings(correctedabs, cuvette, extinction, mgml):
    protein_conc = correctedabs / float(cuvette) / float(extinction) * 1000000
    protein_mgml = protein_conc * float(mgml) / 1000
    return protein_conc, protein_mgml


def generate_plot(x_log, y_log, baseline, protein_conc, mgml, r_value, abs280, filename, wavelength):
    import os

    filename = os.path.split(filename)[1]

    plt.xlabel('log(Wavelength)')
    plt.ylabel('log(Absorbance)')
    plt.plot(x_log, y_log, label='log(Absorbance)')
    plt.plot(x_log, baseline, label='log(Baseline)')
    plt.legend(loc=1)  # upper right corner, anticlockwise
    plt.suptitle(filename, size=11)  # titlesize and name
    # plt.rcParams["figure.figsize"] = (15, 15)
    if mgml != 0:
        plt.annotate(
            'Protein concentration:  ' + str(
                round(protein_conc, 2)) + ' uM\n                                     ' + str(
                round(mgml, 2)) + ' mg/ml\n'
                                  '$R^{2}$ = ' + str(round(-r_value, 3)),
            xy=(np.log10(wavelength), log_y[int((np.where(log_x == np.log10(wavelength))[0]))]),
            xytext=(0.5, 0.7), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="-|>", facecolor='black'))

    if mgml == 0:
        plt.annotate(
            'Protein concentration:  ' + str(round(protein_conc, 2)) + ' uM\n$R^{2}$ = ' + str(round(-r_value, 3)),
            xy=(np.log10(wavelength), log_y[int((np.where(log_x == np.log10(wavelength))[0]))]),
            xytext=(0.5, 0.7), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="-|>", facecolor='black'))

    filename = filename.replace('.csv', '') + str(round(abs280, 3)).replace('.', '') + '.png'
    plt.savefig(filename, dpi=200)  # exports file automatically :)
    plt.show()
    plt.gcf().clear()
    plt.clf()
    plt.cla()
    plt.close('all')


def check_modules():
    try:
        import scipy
        import numpy
        import matplotlib
        import math
        import os
        import pandas
        import easygui
    except ImportError:
        print(
            'Required modules not found (scipy, numpy, matplotlib, pandas and easygui). Install them to run the script. Try pip3 install [modulename]')

def show_graph(data_x, data_y, filename):

    filename = os.path.split(filename)[1]
    plt.xlabel('log(Wavelength)')
    plt.ylabel('log(Absorbance)')
    plt.plot(data_x, data_y, label='log(Absorbance)')
    plt.legend(loc=1)  # upper right corner, anticlockwise
    plt.suptitle(filename, size=11)  # titlesize and name
    # plt.rcParams["figure.figsize"] = (15, 15)
    plt.show()
    plt.gcf().clear()
    plt.clf()
    plt.cla()
    plt.close('all')

def write_settings(cuvette, extintion, mgml):
    with open('correct_abs.conf', 'a') as fh:
        for i in ['cuvette pathlenght', 'extinction coefficient', 'Molar Mass']:
            print(i, file=fh)
        for i in [cuvette,extintion,mgml]:
            print(i,file=fh)

def check_arg():
    normal_run = True
    conc = sys.argv[1:]
    if len(conc) == 0:
        normal_run = True
    if '-g' in conc:
        try:
            normal_run = False
        except IndexError:
            exit(0)
        return normal_run

def read_settings():
    settings = []
    with open('correct_abs.conf', 'r') as fh:
        for i in fh:
            settings.append(i.rstrip('\n'))
    cuvette = float(settings[3])
    extintion = float(settings[4])
    mgml = float(settings[5])
    return cuvette, extintion, mgml


'''Program starts here'''

print('\n\n'
      '============================================================================\n'
      '\tThis script calculates linear regression in log-log plot \n'
      '\tand corrects Absorbance for scattering. You can use it with \n'
      '\texported *.csv files from Cary or any other files with many columns\n\n'
      'To run the script, numpy, matplotlib, easygui, pandas and scipy is required\n'
      '============================================================================\n\n')


check_modules()
normal_run = check_arg()

if normal_run == None:

    # check for the old config file
    if os.path.isfile('correct_abs.conf'):
        msg = 'Do you want to load previous protein settings? (cuvette, extinction coefficient, Molecular mass)'
        title = 'Load previous settings?'
        if easygui.boolbox(msg, title):
            cuvette, extintion, mgml = read_settings()
            old_settings = True
        else:
            old_settings = False
    else:
        old_settings = False

    # pre-settings
    run_again = True
    mutltiple_column = False

    # get the data from file and sample names
    data_names, data, samples, filename = open_file(choose_file())
    data_shape = data.shape[1]

    if data_shape > 2:
        data_column = choose_column(data, data_names, samples)
        mutltiple_column = True
    else:
        data_column = data

    while run_again:
        # #For multiple data column file
        # the log-log data and part for fitting
        log_x, log_y, fitting_region_x, fitting_region_y, x_data, y_data, wavelength = log_log_plot(data_column,
                                                                                                    old_settings)

        # do the baseline
        baseline, abs280, r_value, skip = baseline_correction(log_x, log_y, fitting_region_x, fitting_region_y, x_data,
                                                              y_data, wavelength)

        # calculate concentration
        if not old_settings:
            protein_uM, protein_mgml, cuvette, extintion, mgml = cuvette_extintion_gui(abs280, skip)


        elif old_settings:
            print('\nUsing loaded settings: Pathlength: %f, Extinction coefficient: %d, Molar Mass: %f \n' % (
                round(cuvette, 1), extintion, round(mgml, 2)))
            protein_uM, protein_mgml = calc_protein_conc_from_settings(abs280, cuvette, extintion, mgml)
        # dna_contamination_check(log_x,log_y,baseline) #uncomment to use

        if protein_mgml != 0:
            print('\nProtein concentration: ' + str(round(protein_uM, 2)) + ' uM [or ' + str(
                round(protein_mgml, 2)) + ' mg/ml]. Abs@' + str(wavelength) + ' = ' + str(
                round(abs280, 4)) + ' with pathlength: ' + str(round(float(cuvette), 2)) + ' cm')
        else:
            print('\nProtein concentration: ' + str(round(protein_uM, 2)) + ' uM. Abs@' + str(wavelength) + ' = ' + str(
                round(abs280, 4)) + ' with pathlength: ' + str(round(float(cuvette), 2)) + ' cm')

        # plot the graph. Abs@'+str(wavelength)+' = '+str(abs280)
        generate_plot(log_x, log_y, baseline, protein_uM, protein_mgml, r_value, abs280, filename, wavelength)

        if mutltiple_column:
            data_column = choose_column(data, data_names, samples)
        else:
            run_again = False

    # saving settings
    if not os.path.isfile('correct_abs.conf'):
        msg = 'Do you want to save current protein settings? (cuvette, extinction coefficient, Molecular mass)'
        title = 'Save current settings?'
        if easygui.boolbox(msg, title):
            write_settings(cuvette, extintion, mgml)
            print('\nSaving setting to correct_abs.conf')

#for showing the plot only
elif not normal_run:
    # pre-settings
    run_again = True
    mutltiple_column = False

    # get the data from file and sample names
    data_names, data, samples, filename = open_file(choose_file())
    data_shape = data.shape[1]

    if data_shape > 2:
        data_column = choose_column(data, data_names, samples)
        mutltiple_column = True
    else:
        data_column = data
    while run_again:
        # #For multiple data column file
        # the log-log data and part for fitting
        old_settings = False
        x_data, y_data = get_x_y(data_column)

        show_graph(x_data,y_data, filename)
