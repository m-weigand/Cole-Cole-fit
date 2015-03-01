#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Execute to run single Cole-Cole term characterization
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
    rho0_list = np.log(np.array((10, 100)))
    m_list = (0.1, 0.5)
    tau_list = np.log(np.array((0.04, 0.004)))
    c_list = (0.1, 0.8)

    for cc_set in itertools.product(rho0_list,
                                    m_list,
                                    tau_list,
                                    c_list):
        yield cc_set


def _get_frequencies():
    return np.logspace(-3, 4, 20)


def _fit_spectra():
    pwd = os.getcwd()
    for directory in th._get_cc_dirs(testdir):
        os.chdir(directory)
        cmd = 'cc_fit.py -p -c 1 -m 2'
        subprocess.call(cmd, shell=True)
        os.chdir(pwd)


if __name__ == '__main__':
    frequencies = _get_frequencies()
    th._generate_spectra(frequencies, testdir, _generate_cc_sets)
    _fit_spectra()
    th._evaluate_fits(testdir)
