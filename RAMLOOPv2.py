# -*- coding: utf-8 -*-
"""
RAMLOOP Version 2.1.2

NOTE: REQUIRES PYTHON 2.7! WILL NOT RUN ON PYTHON 3.x!

Authors:
George Belsten
Alexanger Hage
"""

from os import chdir
import sys
import csv
import subprocess
import shutil
import signal
import os
import multiprocessing

# Verify version number; requires version 2.7.x
v = sys.version_info
version = [v[0], v[1]]
if not version == [2, 7]:
    raise Exception("Wrong Python version installed; RAMLOOP requires 2.7.x")

# Attempt to set SIGPIPE to SIG_DFL to prevent broken pipe errors; this will
# not work on Windows, however it may function by being run in an IDE
try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except AttributeError:
    pass


def working_directory(ram_dir=os.getcwd()):

    """
    Checks the given directory (interpreter current directory by default) for
    the ram.exe and, if it does not contain it, prompts the user to specify
    a valid directory.

    :param ram_dir: directory path
    :type ram_dir: str
    :return: directory path
    """

    # Verify that ram.exe is in the specified directory; if not, require the
    # user to update the directory
    while not os.access(ram_dir + ".\\ram.exe", os.R_OK):
        print "Directory not valid; does not contain ram.exe"
        ram_dir = raw_input("Specify valid directory: ")

    return ram_dir


def input_parameters(edit=True):

    """
    Contains the default RAM input parameters, and allows for adjustment in
    the interpreter if edit = True.

    :param edit: if True, allows for editing the default RAM parameters in
        the interpreter
    :type edit: bool
    :return: list of RAM input parameters
    """

    # Default list of parameters; see __print_params for a description of each
    # parameter in the list
    params = [
        5.0, 20.0, 10000.0, 10.0, 1, 8000.0, 1.0, 1, 200.0, 1600.0, 8,
        1.0, 0.0,
        [0.0, 25000.0, 50000.0], [2000.0, 2000.0, 2000.0],
        [0.0], [1700.0],
        [0.0], [1.5],
        [0.0, 900, 1000], [0.5, 0.5, 10.0]
    ]

    def __print_params():
        # Print the list of parameters to the interpreter
        print "-------------------------------------------------------"
        print "0 - Source Depth: %dm" % params[0]
        print "1 - Receiver Depth: %dm" % params[1]
        print "2 - Simulation Maximum Range: %dm" % params[2]
        print "3 - Simulation Range Step: %dm" % params[3]
        print "4 - Range Decimation Factor: %d" % params[4]
        print "5 - Simulation Depth: %dm" % params[5]
        print "6 - Simulation Depth Step: %dm" % params[6]
        print "7 - Depth Decimation Factor: %d" % params[7]
        print "8 - tl.grid Depth Maximum Output: %dm" % params[8]
        print "9 - Reference Sound Speed: %dm/s" % params[9]
        print "10 - Number of rational approximation terms: %d" % params[10]
        print "11 - Number of stability constraints (1 or 2): %d" % params[11]
        print "12 - Radius of stability constraints: %d" % params[12]
        print "13 - Bathymetry profile (range list): %sm" % \
              "m, ".join(str(x) for x in params[13])
        print "14 - Bathymetry profile (depth list): %sm" % \
              "m, ".join(str(x) for x in params[14])
        print "15 - Sediment sound profile (depth list): %sm" % \
              "m, ".join(str(x) for x in params[15])
        print "16 - Sediment sound profile (sound speed list): %sm/s" % \
              "m/s, ".join(str(x) for x in params[16])
        print "17 - Sediment density profile (depth list): %sm" % \
              "m, ".join(str(x) for x in params[17])
        print "18 - Sediment density profile (density value list): %s" % \
              ", ".join(str(x) for x in params[18])
        print "19 - Sediment attenuation profile (depth list): %sm" % \
              "m, ".join(str(x) for x in params[19])
        print "20 - Sediment attenuation profile (attenuation value): %s" % \
              ", ".join(str(x) for x in params[20])

    def __edit_params():
        # Show interface in interpreter for user to edit parameters
        while True:
            print "-------------------------------------------------------"
            edit = raw_input("Edit parameter values? (y/n): ")
            if edit == "n":
                break
            elif edit == "y":
                edit_select = raw_input("Enter index of parameter to edit: ")
            else:
                print "Invalid input"
                continue

            try:
                edit_select = int(edit_select)
            except:
                print "Invalid index input format"
                continue

            # Values 0-12 are floats, 13-20 are lists; handle differently
            if 0 <= edit_select <= 12:
                edit_value = raw_input("Enter new value for parameter: ")
                try:
                    edit_value = float(edit_value)
                except:
                    print "Invalid parameter input format"
                    continue

            elif 13 <= edit_select <= 20:
                edit_value = raw_input("Enter new list for parameter, "
                    "delimited by commas (without brackets): ")
                edit_value = edit_value.split(",")
                for i, a in enumerate(edit_value):
                    try:
                        edit_value[i] = float(a)
                    except:
                        print "Invalid parameter input format"
                        continue

            else:
                print "Invalid parameter index"
                continue

            # Update params list and print new list to console
            params[edit_select] = edit_value
            __print_params()

    if edit:
        # Load edit interface if edit parameter is left as True
        print "_______________________________________________________"
        print "RAM INPUT PARAMETERS"
        print "Verify and adjust parameters as necessary"
        print "Consult RAM documentation for description of parameters"
        __print_params()
        __edit_params()

    # Verify that the connected lists are of equal length, or RAM will not run
    if not len(params[13]) == len(params[14]):
        print "Bathymetry profile range/depth lists must be equal length"
        __edit_params()
    if not len(params[15]) == len(params[16]):
        print "Sediment profile depth/sound speed lists must be equal length"
        __edit_params()
    if not len(params[17]) == len(params[18]):
        print "Sediment density profile depth/density lists must be equal " \
              "length"
        __edit_params()
    if not len(params[19]) == len(params[20]):
        print "Sediment attenuation depth/att lists must be equal length"
        __edit_params()

    return params


def run_parameters(edit=True):

    """
    Contains the default RAM running parameters, and allows for adjustment in
    the interpreter if edit = True.

    :param edit: if True, allows for editing the default parameters in
        the interpreter
    :type edit: bool
    :return: list of running parameters
    """

    # Default parameters: day range, frequency range, grid output, thread count
    params = [
        range(1827),
        [50, 101, 1],
        False,
        multiprocessing.cpu_count() - 1
    ]

    def __print_params():
        # Print the list of parameters to the interpreter
        print "-------------------------------------------------------"
        if len(params[0]) > 10:
            print("0 - Day List: [%d, %d, ..., %d, %d]" %
                  (params[0][0], params[0][1], params[0][-2], params[0][-1]))
        else:
            print("0 - Day List: ", params[0])
        print ("1 - Frequencies: range(%d, %d, %d)" %
               (params[1][0], params[1][1], params[1][2]))
        print "2 - Output Grid: ", params[2]
        print "3 - Thread Count: ", params[3]

    def __edit_params():
        # Show interface in interpreter for user to edit parameters
        while True:
            print "-------------------------------------------------------"
            edit = raw_input("Edit parameter values? (y/n): ")
            if edit == "n":
                break
            elif edit == "y":
                edit_select = raw_input("Enter index of parameter to edit: ")
            else:
                print "Invalid input"
                continue

            try:
                edit_select = int(edit_select)
            except:
                print "Invalid index input format"
                continue

            # 0 and 1 are range parameters (requiring 3 values); 2 is bool;
            # 3 is an integer; handle inputs differently
            if edit_select == 0:
                day_list = raw_input("Enter list of days separated by "
                                     "commas: ")
                try:
                    day_list = day_list.split(",")
                    for i, x in enumerate(day_list):
                        day_list[i] = int(x)
                except:
                    print "Invalid input format"
                    continue
                day_list.sort()
                params[0] = day_list

            elif edit_select == 1:
                range_params = raw_input("Enter new range start/stop/step "
                    "parameters, separated by commas (all 3 are needed): ")
                try:
                    range_params = range_params.split(",")
                    for i, x in enumerate(range_params):
                        range_params[i] = int(x)
                except:
                    print "Invalid range input parameter format"
                    continue
                if not len(range_params) == 3:
                    print "Three parameters are required (start/stop/step)"
                    continue
                params[1] = range_params

            elif edit_select == 2:
                g_o = raw_input("Enter True/False for grid output: ")
                if g_o == "True":
                    params[2] = True
                elif g_o == "False":
                    params[2] = False
                else:
                    print "Invalid input; must be True/False"
                    continue

            elif edit_select == 3:
                th = raw_input("Enter number of processing threads to use: ")
                try:
                    th = int(th)
                except:
                    print "Invalid input parameter; must be positive integer"
                    continue
                if th > 0:
                    params[3] = th
                else:
                    print "Invalid input; must be positive integer"

            else:
                print "Invalid parameter index"
                continue

            __print_params()

    if edit:
        # Load edit interface if edit parameter is left as True
        print "_______________________________________________________"
        print "SIMULATION PARAMETERS"
        print "Verify and adjust parameters as necessary"
        __print_params()
        __edit_params()

    # Turn parameter lists into range objects
    params[1] = range(params[1][0], params[1][1], params[1][2])

    return params


def input_file_names(edit=True):

    """
    Contains the default input file names, and allows for adjustment of these
    input files if edit = True.

    :param edit: if True, allows for editing the input files. (Editing also
        enabled if default files are not found.)
    :type edit: bool
    :return: list of file names
    """

    # Default input files; depth, salinity, temperature
    input_files = [
        ".\\input_data\\NS_HYCOM_depth.txt",
        ".\\input_data\\NS_HYCOM_salt.txt",
        ".\\input_data\\NS_HYCOM_temp.txt"
    ]

    def __print_input_files():
        # Print the list of input file paths into the interpreter
        print "-------------------------------------------------------"
        print "0 - Depth data: %s" % input_files[0]
        print "1 - Salinity data: %s" % input_files[1]
        print "2 - Temperature data: %s" % input_files[2]

    def __edit_input_files():
        # Show interface in interpreter for user to edit file paths
        while True:
            print "-------------------------------------------------------"
            edit = raw_input("Edit input file paths? (y/n): ")
            if edit == "n":
                break
            elif edit == "y":
                edit_select = raw_input("Enter index of file name to edit: ")
            else:
                print "Invalid input"
                continue

            try:
                edit_select = int(edit_select)
            except:
                print "Invalid index input format"
                continue

            if 0 <= edit_select <= 2:
                edit_value = raw_input(
                    "Enter new input file name: ")
                try:
                    input_files[edit_select] = edit_value
                except:
                    print "Invalid name input format"
                    continue

            else:
                print "Invalid parameter index"
                continue

    def __verify_input_files(filenames):
        # Verify that the output file paths lead to readable files
        v = 0
        if os.access(filenames[0], os.R_OK):
            if os.access(filenames[1], os.R_OK):
                if os.access(filenames[2], os.R_OK):
                    v = 1

        return v

    if not __verify_input_files(input_files):
        # Check default file paths; if they are wrong, require edit
        print "Default input files could not be verified; please specify " \
              "correct input file names"
        edit = True

    if edit:
        # Load edit interface; verify paths after changing them
        print "_______________________________________________________"
        print "INPUT DATA PATHS"
        print "Verify and adjust file paths as necessary"
        __print_input_files()
        __edit_input_files()
        while not __verify_input_files(input_files):
            print("Could not verify input files; please review and edit")
            __edit_input_files()

    return input_files


def standard_data(depth_list, salinity_list, temperature_list):

    """
    Reads data from the specified input files, returns lists of depths,
    salinities, temperatures and sound speeds (estimated via the Coppens
    equation). Input files are expected in .txt or .csv form with values
    delimited by commas, and data for each day on a separate row. (Depth file
    should be a single row and apply to all days.)

    :param depth_list: specified depth list file
    :type depth_list: str
    :param salinity_list: specified salinity list file
    :type salinity_list: str
    :param temperature_list: specified temperature list file
    :type temperature_list: str
    :return: list of depths, dict of salinities, dict of temperatures, dict of
        sound speeds
    """

    depths = []
    reader = csv.reader(open(depth_list), delimiter=',',
                        quoting=csv.QUOTE_NONE)
    for row in reader:
        for x in row:
            depths.append(float(x))

    salinities = {}
    reader = csv.reader(open(salinity_list), delimiter=',',
                        quoting=csv.QUOTE_NONE)

    for i, row in enumerate(reader):
        s = []
        for x in row:
            s.append(float(x))
        salinities[i] = s

    temperatures = {}
    reader = csv.reader(open(temperature_list), delimiter=',',
                        quoting=csv.QUOTE_NONE)

    for i, row in enumerate(reader):
        t = []
        for x in row:
            t.append(float(x))
        temperatures[i] = t

    sound_speeds = {}

    for x, y in zip(salinities, temperatures):
        speeds = []
        for D, S, T in zip(depths, salinities[x], temperatures[y]):
            t = T/10
            d = D/1000
            v = (1449.05 + (45.7*t) - (5.21 * t**2) + (0.23 * t**3) +
                 (S-35)*(1.333 - (0.126*t) + (0.009 * t**2)) +
                 ((16.23 + (0.253*t))*d) + ((0.213-(0.1*t)) * d**2) +
                 ((0.016 + 0.0002*(S-35))*(S-35)*t*d))
            speeds.append(v)
        sound_speeds[x] = speeds

    return depths, salinities, temperatures, sound_speeds


def refresh_working_directories(threads, frequencies, delete_old_files=True,
                                output_grid=False):

    """
    Prepares the working directory for a new RAM processing run.

    :param threads: number of processing threads (for the RAM workers)
    :type threads: int
    :param delete_old_files: Removes old output file trees
    :type delete_old_files: bool
    :param output_grid: Creates grid output folders
    :type output_grid: bool
    :return: none
    """

    if delete_old_files:
        try:
            shutil.rmtree(".\\ram_output")
        except:
            pass
        os.mkdir(".\\ram_output")
        os.mkdir(".\\ram_output\\ram_ins")
        os.mkdir(".\\ram_output\\line_outputs")
        if output_grid:
            os.mkdir(".\\ram_output\\grid_outputs")

    for f in frequencies:
        os.mkdir(".\\ram_output\\line_outputs\\%d" % f)
        os.mkdir(".\\ram_output\\ram_ins\\%d" % f)
        if output_grid:
            os.mkdir(".\\ram_output\\grid_outputs\\%d" % f)

    if os.path.isdir(".\\worker_threads"):
        shutil.rmtree(".\\worker_threads")

    os.mkdir(".\\worker_threads")
    for i in range(threads):

        os.mkdir(".\\worker_threads\\ram_worker_%d" % i)
        shutil.copyfile(".\\ram.exe",
                        ".\\worker_threads\\ram_worker_%d\\ram.exe" % i)


def __ram_worker_thread(depths, sound_speeds, current_freq, frequencies,
                        current_day, days, params, output_grid, thread_index,
                        shutdown):

    while not shutdown.value:

        with current_day.get_lock():
            with current_freq.get_lock():
                day = current_day.value
                f = current_freq.value
                if current_day.value == days[-1]:
                    current_day.value = days[0]
                    i = frequencies.index(current_freq.value)
                    try:
                        current_freq.value = frequencies[i+1]
                    except IndexError:
                        shutdown.value = 1
                else:
                    i = days.index(current_day.value)
                    current_day.value = days[i+1]

        with open('.\\worker_threads\\ram_worker_%d\\ram.in' %
                  thread_index, 'wb') as raminput:
            ramwriter = csv.writer(raminput, delimiter=' ',
                                   quoting=csv.QUOTE_MINIMAL)
            ramwriter.writerow(
                ['Auto generated ram.in; day %d, frequency %dHz' %
                 (day, f)])
            ramwriter.writerow(['%f' % f, '%f' % params[0],
                                '%f' % params[1]])
            ramwriter.writerow(['%f' % params[2], '%f' % params[3],
                                '%d' % params[4]])
            ramwriter.writerow(['%f' % params[5], '%f' % params[6],
                                '%d' % params[7], '%f' % params[8]])
            ramwriter.writerow(['%f' % params[9], '%d' % params[10],
                                '%d' % params[11], '%f' % params[12]])
            for r, d in zip(params[13], params[14]):
                ramwriter.writerow(['%f' % r, '%f' % d])
            ramwriter.writerow(['-1', '-1'])
            for d, v in zip(depths, sound_speeds[day]):
                ramwriter.writerow(['%f' % d, '%f' % v])
            ramwriter.writerow(['-1', '-1'])
            for d, s in zip(params[15], params[16]):
                ramwriter.writerow(['%f' % d, '%f' % s])
            ramwriter.writerow(['-1', '-1'])
            for d, x in zip(params[17], params[18]):
                ramwriter.writerow(['%f' % d, '%f' % x])
            ramwriter.writerow(['-1', '-1'])
            for d, a in zip(params[19], params[20]):
                ramwriter.writerow(['%f' % d, '%f' % a])
            ramwriter.writerow(['-1', '-1'])

        main_dir = os.getcwd()
        chdir(".\\worker_threads\\ram_worker_%d" % thread_index)
        subprocess.call(".\\ram.exe")
        chdir(main_dir)
        shutil.copyfile(".\\worker_threads\\ram_worker_%d\\tl.line" %
                        thread_index,
                        ".\\ram_output\\line_outputs\\%i\\tl%d.line" %
                        (f, day))
        shutil.copyfile(".\\worker_threads\\ram_worker_%d\\ram.in" %
                        thread_index,
                        ".\\ram_output\\ram_ins\\%i\\ram%d.in" %
                        (f, day))
        if output_grid:
            shutil.copyfile(".\\worker_threads\\ram_worker_%d\\tl.grid" %
                            thread_index,
                            ".\\ram_output\\grid_outputs\\%i\\tl%d.grid" %
                            (f, day))


def run_ram(depths, sound_speeds, frequencies, days, thread_count, params,
            output_grid):

    """
    Runs RAM with the specified parameters and data.

    :param depths: list of depths
    :type depths: list
    :param sound_speeds: dictionary containing sound speed lists for each day
    :type sound_speeds: dict
    :param frequencies: list of frequencies
    :type frequencies: list
    :param days: list of days
    :type days: list
    :param thread_count: number of RAM threads to run
    :type thread_count: int
    :param params: RAM input parameters
    :type params: list
    :param output_grid: Specify True to output tl.grid files
    :type output_grid: bool
    :return: none
    """

    current_freq = multiprocessing.Value('i', frequencies[0])
    current_day = multiprocessing.Value('i', days[0])
    shutdown = multiprocessing.Value('i', 0)
    threads = []
    for i in range(thread_count):
        p = multiprocessing.Process(target=__ram_worker_thread, args=(
            depths, sound_speeds, current_freq, frequencies, current_day, days,
            params, output_grid, i, shutdown))
        threads.append(p)

    print "Starting worker threads"
    for x in threads:
        x.start()

    for x in threads:
        x.join()

    print "Finished!"
    shutil.rmtree(".\\worker_threads")


def copy_output(ram_dir=os.getcwd()):

    """
    Asks the user if they want to copy the RAM output files to a different
    directory, and if so, to specify the new directory.

    :param ram_dir: working directory for RAM
    :type ram_dir: str
    :return: none
    """

    while True:
        c = raw_input("Copy output files to new directory? (y/n)")
        if c == "y":
            output_dir = raw_input("Enter save directory")
            try:
                os.mkdir(output_dir)
            except:
                print "Invalid save directory"
                continue
            try:
                shutil.copytree(ram_dir + "\\ram_output\\ram_ins",
                                output_dir + "\\ram_ins")
                shutil.copytree(ram_dir + "\\ram_output\\line_outputs",
                                output_dir + "\\line_outputs")
                if os.path.isdir(ram_dir + "\\ram_output\\grid_outputs"):
                    shutil.copytree(ram_dir + "\\ram_output\\grid_outputs",
                                    output_dir + "\\grid_outputs")
            except:
                print "Cannot copy files to specified directory"
                continue
            print "Output data successfully copied to ", output_dir
            break
        elif c == "n":
            break
        else:
            print "Invalid input"
            continue


if __name__ == '__main__':
    print "RAMLOOP V2.1.1"
    ram_dir = working_directory()
    chdir(ram_dir)
    days, frequencies, output_grid, thread_count = run_parameters(True)
    print "_______________________________________________________"
    print "WARNING! All data in old ram_ins, line_outputs and grid_outputs " \
        "folders will be deleted! Copy this data elsewhere if you want to " \
        "keep it!"
    raw_input("Press Enter to continue")
    file_names = input_file_names(True)
    refresh_working_directories(thread_count, frequencies, True, output_grid)
    ram_params = input_parameters(True)
    d, S, t, v = standard_data(file_names[0], file_names[1], file_names[2])
    run_ram(d, v, frequencies, days, thread_count, ram_params, output_grid)
    copy_output(ram_dir)
