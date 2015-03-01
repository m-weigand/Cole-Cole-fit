#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2012-2015  Maximilian Weigand, 2012-2013 Gunnar Jansen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Command line wrapper for cc_fit

Todo
----
* restructure output methods to only take an output directory

END DOCUMENTATION
"""
import os
import lib_cc_fit.cc_fit as cls_cc_fit
import numpy as np
from optparse import OptionParser


def handle_cmd_options():
    description = """
Copyright (C) 2012-2015  Maximilian Weigand, 2012-2013 Gunnar Jansen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

    parser = OptionParser(description=description)

    parser.add_option("-f", "--frequency_file", dest="frequency_file",
                      type='string',
                      help="Frequency file (default: frequencies.dat)",
                      metavar="FILE", default='frequencies.dat')
    parser.add_option("--ignore", dest="ignore", type="string",
                      help="Frequencies indices to ignore. Separate by ';' " +
                      "and start with 1. For example:\"1;2;4;5\". " +
                      "Also allowed are ranges:\"2-10\", and open ranges:" +
                      "\"5-\" (default:-1 (all))",
                      default=-1)
    parser.add_option("-d", "--data_file", dest="data_file", type='string',
                      help="Data file (default: data.dat)",
                      metavar="FILE", default='data.dat')
    parser.add_option("-c", "--nr_cc_terms", dest="nr_cc_terms", type='int',
                      help="Number of Cole-Cole terms (default: 2)",
                      metavar="INT", default=2)
    parser.add_option("-p", "--plot_spectra", action="store_true",
                      dest="plot_spectra",
                      help="If set, plot spectra to spectrum_XX.png files " +
                      "(off by default)", default=False)
    parser.add_option("--mean", action="store_true", dest="save_mean",
                      help="Save mean CC parameters to file", default=False)
    parser.add_option("-o", "--output", type='string', metavar='FILE',
                      help="Output directory (default: results)",
                      default="results", dest="output")
    parser.add_option("--magslog10", action="store_true", dest="magslog10",
                      help="Magnitudes are in log10", default=False)
    parser.add_option("-m", "--starting_model", type="int",
                      dest="starting_model",
                      help="Starting model heuristic (see manual for valid " +
                      "values (default: 2)", default=2)
    parser.add_option("-s", "--init_file", type="string",
                      dest="starting_model_file",
                      help="Load intial CC parameters from file (default: " +
                      "None)", default=None)

    (options, args) = parser.parse_args()

    # extract the frequencies from the ignore string
    if options.ignore != -1:
        frequency_nr = np.loadtxt(options.frequency_file).size
        options.ignore = get_filter_ids(options.ignore, frequency_nr)
        if np.any(options.ignore < 0):
            raise Exception('Only indices larger than 0 are allowed')
        if np.any(options.ignore > (frequency_nr - 1)):
            raise Exception('The maximal allowed index is {0}'.format(
                frequency_nr))
    else:
        options.ignore = []

    return options


def get_filter_ids(filter_string, nr_frequencies=None):
    """
    If nr_frequencies is provided, then range can also have the form "%i-",
    e.g. "4-", and the end index will be set to the largest frequency.
    """
    sections = filter_string.split(';')

    filter_ids = []
    # now look for ranges and expand if necessary
    for section in sections:
        filter_range = section.split('-')
        if(len(filter_range) == 2):
            start = filter_range[0]
            end = filter_range[1]
            # check for an open range, e.g. 4-
            if(end == ''):
                if(nr_frequencies is not None):
                    end = nr_frequencies
                else:
                    continue
            filter_ids += range(int(start) - 1, int(end))
        else:
            filter_ids.append(int(section) - 1)
    return np.array(filter_ids)

if __name__ == '__main__':
    options = handle_cmd_options()

    cc_fit = cls_cc_fit.cc_fit()
    cc_fit.load_data(options.data_file, options.ignore, options.magslog10)
    cc_fit.load_frequencies(options.frequency_file, options.ignore)
    cc_fit.set_nr_cc_terms(options.nr_cc_terms)

    if options.starting_model_file is not None:
        init_pars = np.atleast_2d(np.loadtxt(options.starting_model_file))
        cc_fit.cc_pars_init = init_pars
    else:
        cc_fit.set_initial_values(options.starting_model)

    cc_fit.fit_all_spectra()

    # output
    pwdx = os.getcwd()
    if not os.path.isdir(options.output):
        os.makedirs(options.output)
    os.chdir(options.output)

    if options.plot_spectra:
        cc_fit.plot_all_spectra('./', prefix='')

    if options.save_mean:
        cc_fit.save_mean_of_cc_pars('cc_fits_mean.dat')
    else:
        cc_fit.save_non_essential_results('./')
        cc_fit.save_cc_pars('cc_fits.dat')
        cc_fit.save_cc_errors('cc_fits.dat.err')

    os.chdir(pwdx)
