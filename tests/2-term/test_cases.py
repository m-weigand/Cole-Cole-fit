#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Execute to run double Cole-Cole term characterization
"""
import sys
sys.path.append('..')
import test_helper as th
import os
import itertools
import numpy as np
import subprocess

testdir = 'data'


def _generate_cc_sets():
    """Generate multiple complex resistivity spectra by sampling a certain
    CC-parameter space
    """
    rho0_list = np.log(np.array((10, )))
    # values the first (low-frequency) term is constructed from
    m1_list = (0.05, 0.1)
    tau1_list = np.log(np.array((0.4, 1.0)))
    c1_list = (0.6, 0.8)

    # values for the second term
    m2_list = (0.1, )
    tau2_list = np.log(np.array((0.0004, 0.00001)))
    c2_list = (0.6, )

    for cc_set in itertools.product(rho0_list,
                                    m1_list,
                                    tau1_list,
                                    c1_list,
                                    m2_list,
                                    tau2_list,
                                    c2_list
                                    ):
        yield cc_set


def _get_frequencies():
    return np.logspace(-3, 4, 20)


def _fit_spectra():
    pwd = os.getcwd()
    for directory in th._get_cc_dirs(testdir):
        os.chdir(directory)
        cmd = 'cc_fit.py -p -c 2 -m 2'
        subprocess.call(cmd, shell=True)
        os.chdir(pwd)


if __name__ == '__main__':
    frequencies = _get_frequencies()
    for x in _generate_cc_sets():
        print(x)
    th._generate_spectra(frequencies, testdir, _generate_cc_sets)
    _fit_spectra()
    th._evaluate_fits(testdir)
