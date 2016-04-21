# -*- coding: utf-8 -*-
"""
RAM Data Handler v2.7

Authors:
George Belsten
Alexander Hage

PLEASE NOTE: Script is designed to be imported into a Python interpreter for
use. It is designed to work with output data from the RAMLOOP script.
"""

import os
from os import chdir
import csv
import sys
import ast
import re
import itertools
import subprocess
import shutil
import multiprocessing
import matplotlib
import matplotlib.pyplot as plt
import RAMLOOPv2 as ramloop
import numpy as np
from scipy.interpolate import spline
from matplotlib import gridspec
from datetime import date, timedelta
from statsmodels.nonparametric.smoothers_lowess import lowess
matplotlib.use('Agg')


class ProcessLines:

    def __init__(self,
                 data_dir="D:\college\Year_4\Project\Output\Nova_Scotia",
                 frequencies=range(50, 100, 1),
                 days=range(5113),
                 start_date=[1999, 1, 1],
                 smoothing=0,
                 refresh_averages=False):

        """
        Sets up the data processor.

        :param data_dir: directory containing the line_outputs and ram_ins
            folders from a RAMLOOP session output.
        :type data_dir: str
        :param frequencies: list of frequencies for input data
        :type frequencies: list (int elements)
        :param days: list of days for input data
        :type days: list (int elements)
        :param start_date: date corresponding to first day of data, in list
            format [year, month, day]
        :type start_date: list (int elements)
        :param smoothing: smoothing parameter for line averaging;
            0 = no smoothing
            1 = fine smoothing
            2 = broad smoothing
        :type smoothing: int
        :param refresh_averages: if True, re-calculate average line files
        :type refresh_averages: bool
        """

        # Check directory format is correct
        try:
            chdir(data_dir)
        except:
            raise Exception("Not a valid working directory!")
        self.data_dir = data_dir

        # Check start date format
        if not type(start_date) == list or len(start_date) != 3:
            raise Exception("Invalid start date input type")

        # Generate lists of frequencies and days for internal use while
        # checking that the specified range is valid
        self.frequencies = []
        for i in frequencies:
            if os.path.isdir('line_outputs\\%i' % i):
                self.frequencies.append(i)
            else:
                raise Exception("Frequency range is invalid!")

        # Create (and check) list of days
        self.days = []
        for day in days:
            if os.path.isfile('line_outputs\\%i\\tl%d.line' %
                              (self.frequencies[0], day)):
                self.days.append(day)
            else:
                raise Exception("Day range is invalid!")

        # Generate Ranges list
        self.ranges = []
        line = csv.reader(open('line_outputs\\%i\\tl%d.line' %
                               (self.frequencies[0], self.days[0])),
                          delimiter=' ')
        for row in line:
            # Range given by 3rd value of each row
            self.ranges.append(float(row[3]))

        # Create a dictionary of dates corresponding to each day index
        self.date_dict = {index: date(start_date[0], start_date[1],
                                      start_date[2])+timedelta(days=index)
                          for index in self.days}

        # Generate folder of frequency-averaged tl.line files if they don't
        # already exist
        for day in self.days:
            if refresh_averages:
                print "Refreshing average tl.line files"
                self.gen_avg_lines(smoothing)
                break
            elif not os.path.isfile(".\\line_outputs\\avg\\tl%d.line" %
                                    day):
                print "Producing average tl.line files for data processing"
                self.gen_avg_lines(smoothing)
                break

    def read_line(self, file_name):

        """
        Reads a given tl.line file and outputs the lists of ranges and TLs.

        :param file_name: tl.line file to read
        :type file_name: str
        :return: list of ranges, list of TLs
        """

        ranges = []
        tl = []
        line = csv.reader(open(file_name), delimiter=' ')
        for row in line:
            # Range and TL given by 3rd and 10th value of each row respectively
            ranges.append(float(row[3]))
            tl.append(float(row[10]))

        return ranges, tl

    def read_ram_in(self, file_name):

        """
        Reads a given ram.in file and outputs the list of parameters extracted
        from the file.

        :param file_name: ram.in to read
        :type file_name: str
        :return: list of RAM params, list of ranges, list of sound speeds
        """

        def __append_lists(a, b):
            # Appends a, b to ram_params and returns blank lists to
            # reset them
            ram_params.append(a)
            ram_params.append(b)
            return [], []

        line = csv.reader(open(file_name), delimiter=' ')

        # Append rows to a list
        rows = []
        for row in line:
            rows.append(row)

        # Remove first row (string)
        rows.remove(rows[0])

        # Find locations of -1 dividers
        divider_points = []
        for i, row in enumerate(rows):
            if int(float(row[0])) == -1:
                divider_points.append(i)

        # If list is empty, raise exception
        if not divider_points:
            raise Exception("ram.in file fault")

        # Complex structure required to regenerate original parameter list
        ram_params = []
        ranges = []
        sound_speeds = []

        # First parameter value is frequency which we don't want, so discard
        for x in rows[0][1:]:
            ram_params.append(float(x))

        for row in rows[1: 4]:
            for x in row:
                # Integers must be recovered as ints, but int('float') returns
                # a ValueError; use this to our advantage
                try:
                    ram_params.append(int(x))
                except ValueError:
                    ram_params.append(float(x))
        # Set up lists for appending the list parameters to ram_params
        list0, list1 = [], []

        # Rows 4 to first divider are bathymetry profile, append as lists
        for row in rows[4: divider_points[0]]:
            list0.append(float(row[0]))
            list1.append(float(row[1]))
        list0, list1 = __append_lists(list0, list1)

        # Values between first two dividers are ranges and sound speeds; append
        # to separate lists and keep separate from parameters
        for row in rows[divider_points[0] + 1: divider_points[1]]:
            ranges.append(float(row[0]))
            sound_speeds.append(float(row[1]))

        # Append remaining 3 sets of parameters as lists
        for i in range(3):
            for row in rows[divider_points[1+i] + 1: divider_points[2+i]]:
                list0.append(float(row[0]))
                list1.append(float(row[1]))
            list0, list1 = __append_lists(list0, list1)

        return ram_params, ranges, sound_speeds

    def write_line(self, file_name, line_in, ranges=None):

        """
        Writes a new tl.line file with the specified name and list of TLs.

        :param file_name: new tl.line to write
        :type file_name: str
        :param line_in: list of TLs for line (will use class default ranges)
        :type line_in: list
        :param ranges: list of ranges (defaults to default range list)
        :type ranges: list
        :return: none
        """

        if ranges is None:
            ranges = self.ranges
        linewriter = csv.writer(open(file_name, 'wb'), delimiter=' ',
                                quoting=csv.QUOTE_NONE)

        for r, tl in zip(ranges, line_in):
            # Maintain same formatting as standard tl.line files; convert value
            # to string, add 10 trailing 0's to decimal then limit to 10 char
            r = str(r).ljust(10, '0')[:10]
            tl = str(tl).ljust(10, '0')[:10]
            linewriter.writerow(
                ['', '', '', r, '', '', '', '', '', '', tl, '', '', '', ''])

    def gen_avg_lines(self, smoothing=0, days=None, working_dir=0):

        """
        Reads the tl.line files across all frequencies and produces average
        tl.line files for each day.

        :param days: input day list; if left as default (None), uses standard
            list of days
        :type days: list
        :param working_dir: allows for specific directory to be specified for
            function; else, uses current directory
        :type working_dir: str
        :return: none
        """

        if not working_dir:
            working_dir = os.getcwd()

        # If list of days is not specified, use standard list of days
        if days is None:
            days = self.days

        # Check for line_outputs folder and set base directory
        if os.path.exists(working_dir + "\\line_outputs"):
            base_dir = working_dir + "\\line_outputs\\"
        else:
            raise Exception("Directory does not contain line outputs folder")

        # Read first day file to get line length
        ranges = self.read_line(
            base_dir + "%d\\tl%d.line" % (self.frequencies[0], days[0])
            )[0]
        line_length = len(ranges)
        # Check average folder exists; if not, create it
        if not os.path.exists(base_dir + "\\avg"):
            os.mkdir(base_dir + "\\avg")

        for day in days:
            # Build sums of TLs for each range then divide to take average
            avg_line = [0] * line_length
            for f in self.frequencies:
                line = self.read_line(base_dir + "\\%d\\tl%d.line"
                                      % (f, day)
                                      )[1]
                # Append line values to avg_line list
                for i, value in enumerate(line):
                    avg_line[i] += value
            # Take average
            for i in range(line_length):
                avg_line[i] /= len(self.frequencies)

            if smoothing is not 0:
                if smoothing == 1:
                    fr = 0.01
                elif smoothing == 2:
                    fr = 0.02
                else:
                    raise Exception("Invalid smoothing parameter")

                smooth_input = np.array(avg_line)
                r = np.array(ranges)
                filtered = lowess(
                    smooth_input, r, is_sorted=True, frac=fr, it=0)
                avg_line = filtered[:,1].tolist()

            # Write average line file
            self.write_line(base_dir + "\\avg\\tl%d.line" % day,
                            avg_line, ranges)

    def line_range(self, range_list, line_tl, target_tl):

        """
        Takes an input line (read from a tl.line file) and returns the range
        at which it drops below target TL.

        :param range_list: list of ranges corresponding to line_tl
        :type range_list: list
        :param line_tl: list read from tl.line file
        :type line_tl: list
        :param target_tl: specified TL
        :type target_tl: int
        :return: range
        """

        # Reverse lines that are not already reversed
        if line_tl[0] < line_tl[1]:
            line_tl = reversed(line_tl)

        range_interval = range_list[1] - range_list[0]
        min_range = range_list[0]
        max_range = range_list[-1]
        r = 0

        # Record longest range at which TL is still higher than target
        for i, x in enumerate(line_tl):
            if float(x) < target_tl:
                if not i == len(range_list):
                    r = min((range_list[-i-1]+range_interval), max_range)
                else:
                    pass
                break

        # If target_tl is always exceeded, return minimum range
        if r == 0:
            r = min_range

        return r

    def find_day_tl_range(self, day, target_tl):

        """
        Finds the range at which the TL (averaged over all frequencies)
        reaches a given transmission loss, for a given day.

        :param day: day index for input data
        :type day: int
        :param target_tl: transmission loss to find data for
        :type target_tl: int
        :return: average range
        """

        r_tl, line_tl = self.read_line(".\\line_outputs\\avg\\tl%d.line" % day)
        r = self.line_range(r_tl, line_tl, target_tl)

        return r

    def range_by_day(self, days, target_tl):

        """
        Produces tl_range for a given tl for a given set of days.

        :param days: given list of days
        :type days: list
        :param target_tl: transmission loss to find data for
        :type target_tl: int
        :return: tl_range for each day
        """

        tl = []
        for day in days:
            if type(day) == int:
                r = self.find_day_tl_range(day, target_tl)
                tl.append(r)
            else:
                raise Exception("Day input not int; requires int days in list")

        return tl

    def median(self, input_list):

        """
        Returns the median value of a given list. If the list is of even
        length, it returns the value just above the median (instead of
        interpolating).

        :param input_list: list to find median value
        :type input_list: list
        :return: selected value at median
        """

        input_list = sorted(input_list)
        if len(input_list) == 0:
            raise Exception("Empty list")
        elif len(input_list) == 1:
            return input_list[0]
        # Treat even and odd differently; use remainder function
        elif len(input_list) % 2 == 1:
            index = (((len(input_list)+1)/2)-1)
            return input_list[index]
        elif len(input_list) % 2 == 0:
            # Do not interpolate if list length is even; take next value
            index = ((len(input_list)+1)/2)
            return input_list[index]

    def list_percentile(self, input_list, percentile):

        """
        Returns the value at the given percentile of a given list. (Index for
        selection is rounded to the nearest value.)

        :param input_list: list to find percentile value
        :type input_list: list
        :param percentile: selection percentile (from 0 to 99)
        :type percentile: int
        :return: selected value at percentile
        """

        input_list = sorted(input_list)
        if len(input_list) == 0:
            raise Exception("Empty list")
        elif len(input_list) == 1:
            return input_list[0]
        elif len(input_list) > 1:
            index = int(round(len(input_list)*(percentile/100.)))
            return input_list[index]

    def identify_percentile(self, input_list, value):

        """
        Returns the percentile within the given list for the given value.

        :param input_list: list to find percentile
        :type input_list: list
        :param value: value to select percentile
        :type value: int or float
        :return: percentile of selected value
        """

        input_list = sorted(input_list)
        if len(input_list) == 0:
            raise Exception("Empty list")
        elif len(input_list) == 1:
            return 100
        elif len(input_list) > 1:
            index = input_list.index(value) + 1
            percentile = int(round(index*100./len(input_list)))
            return percentile

    def generate_pdf(self, days, target_tl, ranges=None):

        """
        Returns a dictionary containing a list of days that return the range
        for a given target_tl.

        :param days: list of days
        :type days: list
        :param target_tl: transmission loss to find data for
        :type target_tl: int
        :param ranges: overwrites default list of ranges to use
        :type ranges: list
        :return: dictionary containing list of days for each range
        """

        if ranges is None:
            ranges = self.ranges

        tl_ranges = self.range_by_day(days, target_tl)

        data_dict = {int(r): [] for r in ranges}
        for day, r in zip(days, tl_ranges):
            data_dict[r].append(day)

        return data_dict

    def flatten_pdf(self, pdf):

        """
        Flattens a pdf data dict from generate_pdf to produce a full list of
        ranges (without days) for median/percentile selection.

        :param pdf: data_dict from generate_pdf
        :type pdf: dict
        :return: list of ranges (incl. duplicates)
        """

        range_list = []
        for r in pdf:
            range_list.extend([r for i in range(len(pdf[r]))])

        range_list = sorted(range_list)
        return range_list

    def mean_range(self, days, target_tl):

        """
        Returns the mean range at which a given transmission loss is reached,
        for a specified list of days.

        :param days: specified list of days
        :type days: list
        :param target_tl: transmission loss to find data for
        :type target_tl: int
        :return: mean range
        """

        data_dict = self.generate_pdf(days, target_tl)
        avg = 0
        # Sum each range multiplied by number of occurrences of that range
        for r in data_dict.keys():
            avg += r*len(data_dict[r])
        avg /= len(days)
        return avg

    def percentile(self, days, target_tl, percentile, ranges=None):

        """
        Returns the range at the given percentile for the data given by the
        given transmission loss, the list of days for which this range was
        calculated, and a list of days above this percentile.

        :param days: specified list of days
        :type days: list
        :param target_tl: transmission loss to find data for
        :type target_tl: int
        :param percentile: percentile to find
        :type percentile: int

        :return: range, days for range, days above range
        """

        if ranges is None:
            ranges = self.ranges

        pdf = self.generate_pdf(days, target_tl, ranges)
        range_list = self.flatten_pdf(pdf)
        percentile_range = self.list_percentile(range_list, percentile)

        days = pdf[percentile_range]
        outliers = []
        # Check that percentile range is less than maximum range
        if percentile_range < ranges[-1]:
            # Append lists of days for all greater ranges
            for r in ranges[(ranges.index(percentile_range)+1):]:
                outliers.extend(pdf[r])

        return percentile_range, days, outliers

    def vertical_binning(self, bin_days, bin_range):

        """
        Takes a bin of days and returns the day with the lowest TL at a given
        range.

        :param bin_days: list of bin days
        :type bin_days: list
        :param bin_range: range of the bin, in m
        :type bin_range: float
        :return: day of lowest TL, lowest TL
        """

        min_tl = 500.
        min_tl_day = -1
        # Increment max TL and the day if a day is higher than the current
        # values
        for day in bin_days:
            ranges, tls = self.read_line(
                ".\\line_outputs\\avg\\tl%d.line" % day)
            try:
                i = ranges.index(bin_range)
            except ValueError:
                raise Exception("Invalid bin range; day %d, bin_range %d" %
                                (day, bin_range))

            tl = tls[i]
            if tl < min_tl:
                min_tl = tl
                min_tl_day = day

        if min_tl_day == -1:
            raise Exception("Function error")

        return min_tl_day, min_tl

    def day_sound_profile(self, days):

        """
        Returns average sound speed profile for a given day (or days).

        :param days: given list of days
        :type days: list
        :return: list of depths, list of sound speeds
        """

        # Convert days to list if int
        if type(days) != list:
            try:
                days = [int(days)]
            except:
                raise Exception("Invalid input value")

        ranges = []
        sound_speeds = []
        for day in days:
            file_path = ".\\ram_ins\\%d\\ram%d.in" % (self.frequencies[0], day)
            null, r, s = self.read_ram_in(file_path)
            # Set ranges and initial sound speed list on first loop
            if not ranges:
                ranges = r
                sound_speeds = s
            # Add new speeds to old list for each further loop
            else:
                sound_speeds = [x+y for x, y in zip(sound_speeds, s)]

        sound_speeds = [i/len(days) for i in sound_speeds]

        return ranges, sound_speeds

    def plot_sound_profile(self, days):

        """
        Function to easily generate sound speed plots for a given day range.

        :param days: given day range
        :type days: int list
        :param graph_title: title for graph
        :type graph_title: str
        :return: none
        """

        ranges, speeds = self.day_sound_profile(days)

        y = np.linspace(ranges[0], 500)
        x = spline(ranges, speeds, y)

        plt.plot(x, y)

        plt.title('Sound Speed Profile Day %i' % days[0])
        plt.gca().invert_yaxis()
        plt.xlabel('Sound speed (m/s)')
        plt.ylabel('Depth (m)')
        plt.grid()
        plt.show()

    def find_outliers(self, min_range, target_tl):

        """
        Finds all days with a range above given range for target_tl.

        :param min_range: given minimum outlier range in m
        :type min_range: int or float
        :param target_tl: target TL in dB
        :type target_tl: int or float
        :return: list of days
        """

        days = []
        for day in self.days:
            day_range = self.find_day_tl_range(day, target_tl)
            if day_range >= min_range:
                days.append(day)

        return days

    def plot_percentile_sound_profile(self, target_tl, percentile):

        """
        Plots the sound profile for all days above a given percentile.

        :param target_tl: target TL in dB
        :type target_tl: int or float
        :param percentile: percentile to use
        :type percentile: int
        :return: no return; gives a plot
        """

        days = self.percentile(
            self.days, target_tl, percentile)[2]
        self.plot_sound_profile(days,
                                "Sound profile for days beyond %d percentile" %
                                percentile)

    def sound_peak_size(self, day):

        """
        Returns size of initial low-depth peak in sound speed, if peak exists.

        :param day: specified day
        :type day: int
        :return: size of sound speed peak (m/s), depth of turning point(m)
        """

        ranges, sound_speeds = self.day_sound_profile(list(day))

        for i, s in enumerate(sound_speeds[1:]):
            if s < sound_speeds[i]:
                if s+3 <= sound_speeds[0]:
                    raise Exception("No initial sound profile peak!")
                elif s <= sound_speeds[0]:
                    pass
                else:
                    return sound_speeds[i] - sound_speeds[0], ranges[i]

    def monthly_bins(self, target_tl, percentile):

        """
        Finds the binned and outlier days for a given percentile for each
        month.

        :param target_tl: target tl
        :type target_tl: int
        :param percentile: target percentile
        :type percentile: int
        :return: dictionary of [[binned], [outlier]] days for each month
        """

        # Sort all of the days into months
        sort_days = {i: [] for i in range(1, 13)}
        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        output = {i: [[], []] for i in range(1, 13)}
        # Find binned days for each month
        for b in sort_days:
            input_days = sort_days[b]
            null, bin_days, outlier_days = self.percentile(
                input_days, target_tl, percentile)
            output[b] = [bin_days, outlier_days]

        return output

    def ram_averages(self, ram_dir, target_tl, percentile,
                     ram_params=None, orig_dir=os.getcwd()):

        """
        Runs RAM for monthly average sound profiles taken from the bin around
        the given percentile in the probability density function.

        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param target_tl: specified TL
        :type target_tl: int
        :param percentile: specified percentile
        :type percentile: int
        :param ram_params: specified input parameters for RAM; if left as
            default (None), runs ramloop method to specify params
        :type ram_params: list
        :param orig_dir: allows for specifying original data directory;
            defaults to original directory
        :type orig_dir: str
        :return: none
        """

        day_dict = self.monthly_bins(target_tl, percentile)
        # Produce average sound speed profiles from binned days
        month_avg_sound = {i: [] for i in range(1, 13)}
        for x in day_dict:
            bin_days = day_dict[x][0]
            sound_speeds = self.day_sound_profile(bin_days)[1]
            month_avg_sound[x] = sound_speeds
        depths = self.day_sound_profile(day_dict[1][0])[0]
        # Default run parameters for RAM
        run_params = [range(1, 13), self.frequencies, False,
                      multiprocessing.cpu_count() - 1]
        if ram_params is None:
            ram_params = self.read_ram_in(
                ".\\ram_ins\\%d\\ram%d.in" %
                (self.frequencies[0], self.days[0])
            )[0]

        chdir(ram_dir)
        print "_______________________________________________________"
        print(
            "WARNING! All data in old ram_ins, line_outputs and "
            "grid_outputs folders will be deleted! Copy this data "
            "elsewhere if you want to keep it!")
        raw_input("Press Enter to continue")
        ramloop.refresh_working_directories(run_params[3], self.frequencies)
        # Run RAM using average sound profile data
        ramloop.run_ram(depths, month_avg_sound, self.frequencies,
                        range(1, 13), run_params[3], ram_params, False)
        # Copy the resulting data to a new folder and name appropriately
        shutil.copytree(".\\ram_output", orig_dir + "\\month_avg_%ddB_%dth" %
                        (target_tl, percentile))
        chdir(orig_dir)

    def ram_avg_sweep(self, ram_dir, tls, percentiles):

        """
        Runs ram_averages for a range of TLs and percentiles.

        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param tls: list of specified TLs
        :type tls: list
        :param percentiles: list of specified percentiles
        :type percentiles: list
        :return: none
        """
        orig_dir = os.getcwd()
        cores = multiprocessing.cpu_count() - 1
        chdir(ram_dir)
        ram_params = ramloop.input_parameters()
        print "_______________________________________________________"
        print(
            "WARNING! All data in old ram_ins, line_outputs and grid_outputs "
            "folders will be deleted! Copy this data elsewhere if you want to "
            "keep it!")
        raw_input("Press Enter to continue")
        ramloop.refresh_working_directories(cores,
                                            self.frequencies)
        chdir(orig_dir)
        for tl in tls:
            for per in percentiles:
                self.ram_averages(ram_dir, tl, per, ram_params,
                                  orig_dir)
                chdir(ram_dir)
                ramloop.refresh_working_directories(cores, self.frequencies)
                chdir(orig_dir)

    def save_day_sound_plots(self, save_dir):

        """
        Plots the sound speed/depth and TL/range profiles for each day.

        :param save_dir: Directory for saving plots
        :type save_dir: str
        """

        # Create output folder if it doesn't exist
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        output_path = save_dir + '\\daily_plots'
        if not os.path.isdir(output_path):
            os.mkdir(output_path)

        # Plot sound profile for each day
        for day in self.days:

            prof = self.day_sound_profile(day)

            fig, ax = plt.subplots(1, 2)
            fig.set_size_inches(19, 10)
            ax[0].plot(prof[1], prof[0])
            ax[0].set_title('Sound profile, day %d' % day)
            ax[0].set_xlabel('Sound speed (m/s)')
            ax[0].set_ylabel('Depth (m)')
            # Limit axes to range that we're interested in
            ax[0].set_xlim(1460, 1520)
            ax[0].set_ylim(500, 0)
            ax[0].grid()

            ax[1].plot(self.ranges, self.read_line(
                ".\\line_outputs\\avg\\tl%d.line" % day
                        )[1]
                       )
            ax[1].set_title('Averaged TL over distance')
            ax[1].set_xlabel('Distance (m)')
            ax[1].set_ylabel('TL (dB)')
            ax[1].set_xlim(0, 10000)
            ax[1].set_ylim(100, 20)
            ax[1].grid()

            plt.savefig(output_path + "\\day_%d.pdf" % day,
                        dpi=100)
            plt.tight_layout()
            plt.clf()

    def save_percentile_sound_plots(self, save_dir, tls):

        """
        Plots the sound speed/depth and TL/range for the median and 90th
        percentile days for each month and specified TL.

        :param save_dir: Directory for saving plots
        :type save_dir: str
        :param tls: list of TLs
        :type tls: list
        """

        # Create output folder if it doesn't exist
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        output_path = save_dir + '\\percentile_plots'
        if not os.path.isdir(output_path):
            os.mkdir(output_path)

        # Sort days into months
        sort_days = {i: [] for i in range(1, 13)}

        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        for tl in tls:
            try:
                os.mkdir(output_path + '\\%ddB' % tl)
            except:
                pass
            for month in range(1, 13):
                m_r, m_bin, null = self.percentile(
                        sort_days[month], tl, 50)
                p_r, p_bin, null = self.percentile(
                    sort_days[month], tl, 90)
                day50 = self.vertical_binning(
                    m_bin, m_r
                )[0]
                day90 = self.vertical_binning(
                        p_bin, p_r
                )[0]
                prof50 = self.day_sound_profile(day50)
                prof90 = self.day_sound_profile(day90)

                fig, ax = plt.subplots(1, 2)
                fig.set_size_inches(19, 10)
                ax[0].plot(prof50[1], prof50[0], 'b-')
                ax[0].plot(prof90[1], prof90[0], 'g-')
                ax[0].set_title('Sound profile')
                ax[0].set_xlabel('Sound speed (m/s)')
                ax[0].set_ylabel('Depth (m)')
                # Limit axes to range that we're interested in
                ax[0].set_xlim(1460, 1520)
                ax[0].set_ylim(500, 0)
                ax[0].grid(which='major')
                ax[0].grid(which='minor')

                ax[1].plot(self.ranges, self.read_line(
                    ".\\line_outputs\\avg\\tl%d.line" % day50
                            )[1],
                           'b-', label='median')
                ax[1].plot(self.ranges, self.read_line(
                    ".\\line_outputs\\avg\\tl%d.line" % day90
                            )[1],
                           'g-', label='90th percentile')
                ax[1].set_title('Averaged TL over distance')
                ax[1].set_xlabel('Distance (m)')
                ax[1].set_ylabel('TL (dB)')
                ax[1].set_xlim(0, 10000)
                ax[1].set_ylim(100, 20)
                ax[1].grid(which='major')
                ax[1].grid(which='minor')
                handles, labels = ax[1].get_legend_handles_labels()
                ax[1].legend(handles, labels)
                plt.tight_layout()

                plt.savefig(output_path + "\\%ddB\\Month %d.pdf" % (tl, month),
                            dpi=100)
                plt.clf()

    def extract_grid_data(self, grid_file):

        """
        Extract data from specified tl.grid file.

        :return: grid data
        """

        # Open file and extract floats and ints
        f = open(grid_file, 'rb')
        grid_data = np.fromfile(f, dtype='float32', count=-1)
        f = open(grid_file, 'rb')
        grid_data_ints = np.fromfile(f, dtype='int', count=-1)

        # Splice specific int values into grid_data
        grid_data[6] = grid_data_ints[6] # range printing freq
        grid_data[9] = grid_data_ints[9] # depth printing freq

        return grid_data

    def process_grid_data(self, grid_data):

        """
        Converts grid data for plotting, using data from extract_grid_data.

        :return: grid data for plot, params
        """

        # params are freq, src_d, rec_d, max_r, r_step, r_dec, max_d, d_step,
        # d_dec, grid_d
        params = grid_data[1:11]

        ranges = int((params[3]/params[4]) / params[5])
        depths = int((params[9]/params[7]) / params[8])

        grid_out = np.zeros((ranges, depths))
        # Format data from grid_data
        for r in range(ranges):
            for d in range(depths):
                try:
                    grid_out[r, d] = grid_data[
                        18 + (r*(depths+2)) + d
                    ]
                except IndexError:
                    raise Exception("Index error; r=%d, d=%d" % (r, d))

        # Transpose for plotting
        grid_out = np.transpose(grid_out)

        return grid_out, params

    def average_grid_data(self, day, frequencies, data_dir=None):

        """
        Produces an average grid_out from a selection of frequencies.
        """

        # If a new directory is specified, change directory
        orig_dir = os.getcwd()
        if data_dir is not None:
            try:
                chdir(data_dir)
            except:
                raise Exception("Invalid directory")

        # Produce initial grid_out and set params
        grid_data = self.extract_grid_data(
            ".\\grid_outputs\\%d\\tl%d.grid" % (frequencies[0], day)
        )
        grid_out, params = self.process_grid_data(grid_data)

        # Add additional grid data from frequencies
        if len(frequencies) > 1:
            for f in frequencies[1:]:
                grid_data = self.extract_grid_data(
                    ".\\grid_outputs\\%d\\tl%d.grid" % (f, day)
                )
                grid_out += self.process_grid_data(grid_data)[0]

        # Average
        grid_out /= len(frequencies)
        chdir(orig_dir)

        return grid_out, params

    def grid_plot_tool(self, tl_grid_data, params, sound_speed_data,
                       state='show', filename='placeholder'):

        """
        Function used to create 2D plots from tl.grid files. Takes grid data
        and params from grid_data function for a specific frequency or
        averaged across all frequencies.
        """

        max_range = params[3]
        max_grid_depth = params[9]
        aspect = (max_range/max_grid_depth)

        fig = plt.figure()
        fig.set_size_inches(19, 10)
        ax1 = plt.subplot2grid((1, 3), (0, 0))
        ax1.plot(sound_speed_data[1], sound_speed_data[0])
        ax1.set_ylim(max_grid_depth, 0)
        ax1.set_xlim(1460, 1520)
        ax1.grid()
        ax1.set_ylabel('Depth (m)')
        ax1.set_xlabel('Sound Speed (m/s)')
        ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=2)
        im = ax2.matshow(
            tl_grid_data,
            interpolation='none',
            extent=[0, max_range, -max_grid_depth, 0],
            aspect=aspect,
            vmin=40,
            vmax=120
        )

        fig.colorbar(
            im,
            orientation='vertical',
            label='Transmission loss (dB re 1m)',
        )
        ax2.set_xlabel('Range (m)')
        ax2.set_ylim(-max_grid_depth, 0)
        plt.tight_layout()

        if state == 'show':
            plt.show()
        elif state == 'save':
            if not os.path.isdir(".\\grid_plots"):
                os.mkdir(".\\grid_plots")
            plt.savefig(".\\grid_plots\\%s.pdf" % filename, dpi=100)
            plt.clf()
        else:
            raise Exception("Invalid show/save state entered")

    def plot_grid(self, day, frequencies, state='show', data_dir=None,
                  filename='placeholder'):

        """
        Extracts data from tl.grid files and plots.
        """

        # If single int frequency entered, turn into single item list
        if type(frequencies) == int:
            frequencies = [frequencies]

        tl_grid_data, params = self.average_grid_data(day, frequencies,
                                                      data_dir)
        sound_speed_data = self.day_sound_profile(day)

        self.grid_plot_tool(tl_grid_data, params, sound_speed_data, state,
                            filename)

    def interesting_grid_plots(self, frequencies):

        """
        Tool for automatically plotting relevant grid plots
        """

        days = [2946, 2966, 2999, 3028, 3063, 3095, 3125, 3139, 3172, 3200,
                3243, 3256, 2969, 2998, 3041, 3047, 3087, 3118, 3147, 3184,
                3214, 3244, 3264]

        for day in days:
            print "Plotting grid for day %d" % day
            self.plot_grid(day, frequencies, state='save',
                           filename="Grid plot, day %d, frequencies %d-%d" %
                                    (day, frequencies[0], frequencies[-1]))

    def save_monthly_data(self, save_dir, tls):

        sort_days = {i: [] for i in range(1, 13)}

        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        save_file = open(save_dir + '\\tl_data.csv', 'wb')

        if not os.path.isdir(save_dir + "\\monthly_data_graphs"):
            os.mkdir(save_dir + "\\monthly_data_graphs")

        with save_file as data_input:
            data_writer = csv.writer(data_input, delimiter=',',
                                     quoting=csv.QUOTE_MINIMAL)
            colours = ['b', 'c', 'g', 'y', 'r', 'm', 'w', 'k']

            # Create and flatten line to write
            line = [["%idB" % tl] + [" "]*(len(tls)+2) for tl in tls]
            line = [x for y in line for x in y]
            data_writer.writerow(line)

            line = [["Month", "Day"] + ["%idB" % tl for tl in tls] + [" "]
                    for tl in tls]
            line = [x for y in line for x in y]
            data_writer.writerow(line)

            for percentile in [50, 90]:
                data_writer.writerow(["%dth Percentile" % percentile])
                plot_data = {tl: {month: [] for month in range(1, 13)}
                             for tl in tls}
                for month in range(1, 13):
                    # Get binned days for each tl
                    tl_day = {}
                    for tl in tls:
                        r, day_bin, null = self.percentile(
                            sort_days[month], tl, percentile
                        )
                        select_day = self.vertical_binning(
                            day_bin, r
                        )[0]
                        tl_day[tl] = select_day
                    # Generate full list of ranges for each tl
                    ranges = [[self.find_day_tl_range(tl_day[tl], db)
                              for db in tls]
                              for tl in tls]
                    # Produce line for writing
                    line = [[month, tl_day[tl]] + range_row + [" "]
                            for tl, range_row in zip(tls, ranges)]
                    # Flatten list
                    line = [x for y in line for x in y]
                    data_writer.writerow(line)

                    # Repeat with percentiles
                    percentiles = [
                        [self.identify_percentile(
                            self.flatten_pdf(
                                self.generate_pdf(sort_days[month], db)
                            ),
                            self.find_day_tl_range(tl_day[tl], db)
                        ) for db in tls
                        ] for tl in tls]

                    line2 = [
                        [" ", "Percentiles"] + row + [" "]
                        for row in percentiles
                    ]
                    line2 = [x for y in line2 for x in y]
                    data_writer.writerow(line2)

                    # Save ranges to plot_data for bar graph plotting
                    for i, row in enumerate(ranges):
                        plot_data[tls[i]][month] = row

                for tl_select_index, tl in enumerate(tls):
                    bars = {month: [] for month in range(1, 13)}
                    plot_list = [[] for db in tls]
                    for month in range(1, 13):
                        for i, db in enumerate(tls):
                            value = plot_data[db][month][tl_select_index]
                            bars[month].append(value)
                            plot_list[i].append(value)

                    ind = np.arange(12)
                    width = 0.8/len(tls)
                    fig, ax = plt.subplots()
                    fig.set_size_inches(19, 10)
                    rects = []
                    for i, bars in enumerate(plot_list):
                        r = ax.bar(ind + i*width, bars, width,
                                   color=colours[i])
                        rects.append(r)

                    ax.set_title("Ranges for binned days at %ddB" % tl)
                    ax.set_ylabel("Range")
                    ax.set_xticks(ind + 0.4)
                    ax.set_xticklabels(
                        ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                         "Aug", "Sep", "Oct", "Nov", "Dec"])
                    legend = ["%ddB selected bin" % db for db in tls]
                    ax.legend([r[0] for r in rects], legend)
                    ax.grid()
                    output_path = \
                        save_dir + \
                        "\\monthly_data_graphs\\%dth_percentile_%ddB.pdf" % \
                        (percentile, tl)
                    plt.tight_layout()
                    plt.savefig(output_path, dpi=100)
                    plt.clf()


##############################################################################
                ### Below this is Alex' code for various tests###
##############################################################################

    def true_percentiles(self, ram_dir, tls, orig_dir=0):

        """
        Finds the days corresponding to the p90 and p50 of the wide-spaced data,
        then runs those days through RAMLOOP with finer spacing.

        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param tls: list of TLs for which to run RAMLOOP
        :type tls: list
        :param orig_dir: allows for specifying original data directory;
            defaults to original directory
        :type orig_dir: str
        :return: none
        """

        # Putting this in the parameters meant it was called after the directory
        # was changed to ram_dir, breaking the loop.
        if orig_dir == 0:
            orig_dir=os.getcwd()

        sort_days = {i: [] for i in range(1, 13)}

        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        ram_params = ramloop.input_parameters(edit=False)
        ram_params[3] = 1
        cores = multiprocessing.cpu_count() - 1

        for tl in tls:
            for month in range(1,13):
                for i in range(50, 100, 40):
                    # gets bin of days for percentile
                    none, days, none2 = self.percentile(sort_days[month], tl,
                                                        i)
                    # creates dictionary of sound speeds per day
                    sound_profiles ={i: [] for i in range(len(days))}
                    for x in range(len(days)):
                        sound_speeds = self.day_sound_profile(days[x])[1]
                        sound_profiles[x] = sound_speeds
                    depths = self.day_sound_profile(days[0])[0]

                    chdir(ram_dir)
                    ramloop.refresh_working_directories(cores,
                                                        self.frequencies)
                    # Run RAM using binned days sound profile data
                    ramloop.run_ram(depths, sound_profiles, self.frequencies,
                                    range(len(days)), cores, ram_params,
                                    False)
                    # Copy the resulting data to a new folder
                    shutil.copytree(".\\ram_output", orig_dir +
                                    "\\fine\\month%d_%ddB_%dth" %(month, tl, i))
                    chdir(orig_dir)

    def attenuation_test(self, constants, days, ram_dir, orig_dir=0):
        """
        Runs RAMLOOP for several different attenuation constants and saves the
        data in subfolders.
        :param constants: list of attenuation constants for which to run RAMLOOP
        :type constants: list
        :param days: days for which to run RAMLOOP
        :type days: list
        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param orig_dir: allows for specifying original data directory;
            defaults to original directory
        :type orig_dir: str
        :return: none
        """

        if orig_dir == 0:
            orig_dir=os.getcwd()

        ram_params = ramloop.input_parameters(edit=False)
        cores = multiprocessing.cpu_count() - 4

        for constant in constants:
            ram_params[20] = [constant, constant, 10]
            sound_profiles ={i: [] for i in range(len(days))}
            for x in days:
                sound_speeds = self.day_sound_profile(x)[1]
                sound_profiles[days.index(x)] = sound_speeds
            depths = self.day_sound_profile(days[0])[0]
            chdir(ram_dir)
            ramloop.refresh_working_directories(cores, self.frequencies)
            # Run RAM
            ramloop.run_ram(depths, sound_profiles, self.frequencies,
                            range(len(days)), cores, ram_params, False)
            # Copy the resulting data to a new folder
            shutil.copytree(".\\ram_output", orig_dir +
                            "\\attenuation\\constant%d"
                            %constants.index(constant))
            chdir(orig_dir)

        for i in range(len(constants)):
            chdir(orig_dir +"\\attenuation\\constant%d" %i)
            ranges, tls = self.average_tl(range(len(days)))
            chdir(orig_dir)
            plt.plot(ranges, tls, label=constants[i])

        plt.ylabel('TL /dB')
        plt.xlabel('Range /m')
        plt.legend()
        plt.savefig('attenuation.pdf')

    def density_test(self, densities, days, ram_dir, orig_dir=0):
        """
        Runs RAMLOOP for several different densities and saves the data in
        subfolders.
        :param densities: list of densities for which to run RAMLOOP
        :type densities: list
        :param days: days for which to run RAMLOOP
        :type days: list
        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param orig_dir: allows for specifying original data directory;
            defaults to original directory
        :type orig_dir: str
        :return: none
        """

        if orig_dir == 0:
            orig_dir=os.getcwd()

        ram_params = ramloop.input_parameters(edit=False)
        cores = multiprocessing.cpu_count() - 4

        for density in densities:
            ram_params[18] = [density]
            sound_profiles ={i: [] for i in range(len(days))}
            for x in days:
                sound_speeds = self.day_sound_profile(x)[1]
                sound_profiles[days.index(x)] = sound_speeds
            depths = self.day_sound_profile(days[0])[0]
            chdir(ram_dir)
            ramloop.refresh_working_directories(cores, self.frequencies)
            # Run RAM
            ramloop.run_ram(depths, sound_profiles, self.frequencies,
                            range(len(days)), cores, ram_params, False)
            # Copy the resulting data to a new folder
            shutil.copytree(".\\ram_output", orig_dir +
                            "\\density\\density%d"
                            %densities.index(density))
            chdir(orig_dir)

        for i in range(len(densities)):
            chdir(orig_dir +"\\density\\density%d" %i)
            ranges, tls = self.average_tl(range(len(days)))
            chdir(orig_dir)
            plt.plot(ranges, tls, label=densities[i])

        plt.ylabel('TL /dB')
        plt.xlabel('Range /m')
        plt.legend()
        plt.savefig('density.pdf')


    def source_test(self, sources, days, ram_dir, orig_dir=0):
        """
        Runs RAMLOOP for several different source depths and saves the data in
        subfolders.
        :param sources: list of source depths for which to run RAMLOOP
        :type sources: list
        :param days: days for which to run RAMLOOP
        :type days: list
        :param ram_dir: directory of the RAM program files
        :type ram_dir: str
        :param orig_dir: allows for specifying original data directory;
            defaults to original directory
        :type orig_dir: str
        :return: none
        """

        if orig_dir == 0:
            orig_dir=os.getcwd()

        ram_params = ramloop.input_parameters(edit=False)
        cores = multiprocessing.cpu_count() - 4

        for source in sources:
            ram_params[0] = source
            sound_profiles ={i: [] for i in range(len(days))}
            for x in days:
                sound_speeds = self.day_sound_profile(x)[1]
                sound_profiles[days.index(x)] = sound_speeds
            depths = self.day_sound_profile(days[0])[0]
            chdir(ram_dir)
            ramloop.refresh_working_directories(cores, self.frequencies)
            # Run RAM
            ramloop.run_ram(depths, sound_profiles, self.frequencies,
                            range(len(days)), cores, ram_params, False)
            # Copy the resulting data to a new folder
            shutil.copytree(".\\ram_output", orig_dir +
                            "\\source_depths\\depth_%d"
                            %source)
            chdir(orig_dir)

        for i in range(len(sources)):
            chdir(orig_dir + "\\source_depths\\depth_%d" %sources[i])
            self.monthly_range('source_depth_%d' %i)
            chdir(orig_dir)

    def average_tl(self, days):
        """
        Find the average TL for each range for a specified set of days.

        :param days: specified list of days
        :type days: list (int elements)
        :return: list of ranges, list of TLs
        """
        transmission_losses = [0] * len(self.ranges)
        for f in self.frequencies:
            # Add together the TL for each day, for each range
            for day in days:
                line = csv.reader(open('line_outputs\\%i\\tl%d.line' %
                                       (f, day)), delimiter=' ')
                r = 0
                for row in line:
                    # TL given by 10th row of tl.line file (with space
                    # delimiter); value is a string and so must be floated
                    transmission_losses[r] += float(row[10])
                    r += 1

        # Take the average of the TL for each day (and invert to produce a
        # negative number)
        for i in range(len(transmission_losses)):
            transmission_losses[i] = (-transmission_losses[i] /
                (len(days)*len(self.frequencies)))

        return self.ranges, transmission_losses

    def monthly_range(self, tls, save='y',percentiles=[50, 90],
                      value=ramloop.input_parameters(edit=False)[0]):
        """
        Creates plots of p90/p50 range vs month for different TL's.
        :param tls: list of transmission losses to plot for
        :type tls: list
        :param save: determines if the plot is saved via this function
        :type save: string 'y' to save, anything else to not.
        :param percentiles: percentiles for which to plot
        :type percentiles: list
        :param value: depth of source or receiver
        :type value: float or int
        :return: none
        """
        sort_days = {i: [] for i in range(1, 13)}
        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)
        for i in tls:
            self.__monthly_range_iterate(i, value, sort_days, percentiles)
            plt.xlabel('month')
            plt.ylabel('range /m')
            plt.title('Monthly ranges for target %d dB' %i)
            plt.legend()
            if save == 'y':
                plt.savefig('monthly_comparison_tl%d' %i)

    def __monthly_range_iterate(self, tl, value, days, percentiles):
        """
        :param name: source or receiver
        :type name: string
        :param tl: transmission loss
        :type tl: int
        :param value: depth of source or receiver
        :type value: float or int
        :param days: dictionary of days sorted into months
        :type days: dict
        :param percentiles: percentiles for which to plot
        :type percentiles: list
        :return:
        """
        ranges = {i: [] for i in range(len(percentiles))}
        for month in range(1,13):
            if days[month] == []:
                for i in range(len(percentiles)):
                    ranges[i].append(0)
            else:
                for i in range(len(percentiles)):
                    ranges[i].append(self.percentile(days[month], tl,
                                                     percentiles[i])[0])
        plt.clf()
        for i in range(len(percentiles)):
            plt.plot(range(1,13), ranges[i], label='p%d'
                     % (percentiles[i]))

    def channel_check(self, days, depth=100.):
        """
        Checks if a sound channel exists in the days specified by checking the
        slope of the sound profile for shallow depths.
        :param days: list of days for which to check
        :type days: list
        :param depth: depth up to which to check
        :type depth: float
        :return: fraction of days with sound channels, list of days with sound
        channels, and list of days without.
        """
        channel_days = []
        no_channel_days = []
        for day in days:
            file_path = ".\\ram_ins\\%d\\ram%d.in" % (self.frequencies[0], day)
            null, r, s = self.read_ram_in(file_path)
            depths = [0., 2., 4., 6., 8., 10., 12., 15., 20., 25., 30., 35.,
                      40., 45., 50., 60., 70., 80., 90., 100., 125., 150.,
                      200., 250., 300., 350., 400., 500., 600., 700., 800.,
                      900., 1000., 1250., 1500., 2000.]
            index = depths.index(depth)
            avg_slope = (s[index] - s[0])/depths[index]
            if avg_slope > 0:
                channel_days.append(day)
            else:
                no_channel_days.append(day)

        fraction = float(len(channel_days))/float(len(days))

        return fraction, channel_days, no_channel_days

    def channel_check_percentile(self, target_tl, percentile, depth=100.):
        """
        Checks if the days above the specified percentile contain a sound
        channel.
        :param target_tl: transmission loss for which to calculate percentiles
        :type target_tl: int or float
        :param percentile: percentile above which we check for channels
        :type percentile: int
        :param depth: depth to which we check the sound profile
        :type depth: float
        :return: dictionaries with monthly values for fractions of days with
        sound channels, lists of days with channels, and lists of days without.
        """
        monthly_fractions = {i: 0 for i in range(1, 13)}
        monthly_channel_days = {i: [] for i in range(1, 13)}
        monthly_no_channel_days = {i: [] for i in range(1, 13)}
        sort_days = {i: [] for i in range(1, 13)}
        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        for b in sort_days:
            input_days = sort_days[b]
            null, bin_days, outlier_days = self.percentile(
                input_days, target_tl, percentile)

            total_days = bin_days + outlier_days

            fraction, channel_days, no_channel_days = \
                self.channel_check(total_days, depth)

            monthly_fractions[b] = fraction
            monthly_channel_days[b] = channel_days
            monthly_no_channel_days[b] = no_channel_days

        return monthly_fractions, monthly_channel_days, monthly_no_channel_days

    def consistency_check(self, target_tl, percentile):
        """
        Returns several parameters useful for checking the consistency of the
        percentile sound profile with the days in a higher percentile.
        :param target_tl: transmission loss for which to calculate percentiles
        :type target_tl: int or float
        :param percentile: percentile above which we check for channels
        :type percentile: int
        :return: dictionaries with monthly values for fractions of days with
        ranges within 10% of the percentile range, the days within and outside
        of this 10%, the days with the largest ranges, and dictionaries of the
        ranges of the percentile and worst days.
        """
        monthly_ranges = {i: 0 for i in range(1,13)}
        monthly_fraction_inside = {i: 0 for i in range(1,13)}
        monthly_days_inside = {i: [] for i in range(1, 13)}
        monthly_days_outside = {i: [] for i in range(1,13)}
        worst_days = {i: 0 for i in range(1, 13)}
        max_ranges = {i: 0 for i in range(1, 13)}
        sort_days = {i: [] for i in range(1, 13)}
        for a in self.date_dict:
            month = self.date_dict[a].month
            sort_days[month].append(a)

        for b in sort_days:
            input_days = sort_days[b]
            r, bin_days, outlier_days = self.percentile(
                input_days, target_tl, percentile)

            total_days = bin_days + outlier_days
            extended_range = r*1.1
            monthly_ranges[b] = r

            for i in outlier_days:
                outlier_range = self.mean_range([i], target_tl)
                if outlier_range <= extended_range:
                    monthly_days_inside[b].append(i)
                elif outlier_range > max_ranges[b]:
                    max_ranges[b] = outlier_range
                    worst_days[b] = [i]
                    monthly_days_outside[b].append(i)
                else:
                    monthly_days_outside[b].append(i)

        for b in sort_days:
            monthly_fraction_inside[b] = (float(len(monthly_days_inside[b])) +
                                          len(bin_days)) / len(total_days)

        return monthly_fraction_inside, monthly_days_inside, \
            monthly_days_outside, monthly_ranges, worst_days, max_ranges