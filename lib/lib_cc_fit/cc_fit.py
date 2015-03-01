#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Cole-Cole fit related functions
NOTE: We only fit linear magnitude data in extract.dat files!
"""
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from scipy.optimize import leastsq
import colecole


class cc_fit:
    """
    Procedure for a Cole-Cole fit:
    - load data
    - load frequencies
    - set nr of CC terms
    - set initial parameters:
      - defaults
      - user supplied
    - fit one or more spectra
      - optional: use results of previous spectrum as starting values
    - optional: plot spectra + fits
    - save fits and RMS values to file
    - optional save m / rho0 to file

    TODO:
        - Re-add checking of initial parameters
        - Add option to use previous fit results as starting values
        - Add options to change boundaries etc.
    """

    def __init__(self):
        self.data = None
        self.frequencies = None
        self.nr_cc_terms = 1
        self.cc_pars = None
        # starting values
        self.cc_pars_init = None
        self.magnitude_rms = None
        self.phase_rms = None
        self.errors = None
        self.stdout = []

    def load_data(self, filename, ignore, mags_are_in_log10=False):
        """
        Load magnitude and phase spectra from file.
        """
        tmp_data = self.__load_file(filename)

        # if we have only one spectrum to fit, correct the shape of the array
        datashape = tmp_data.shape
        if(len(datashape) == 1):
            tmp_data.shape = (1, datashape[0])

        # test if we got an id column, and delete if necessary
        if((tmp_data.shape[1] % 2) == 1):
            print('Found ID-column, deleting')
            tmp_data = np.delete(tmp_data, 0, 1)

        # if there is an ignore list, remove the correspoding entries from
        # magnitude and phase
        if(ignore != []):
            ignore_mag_ids = np.array(ignore, dtype=int)
            ignore_pha_ids = ignore_mag_ids + tmp_data.shape[1] / 2 - len(
                ignore)

            # delete mags
            tmp_data = np.delete(tmp_data, ignore_mag_ids, 1)
            tmp_data = np.delete(tmp_data, ignore_pha_ids, 1)

        if mags_are_in_log10 is True:
            tmp_data[:, 0:tmp_data.shape[1] / 2] = 10 ** (
                tmp_data[:, 0:tmp_data.shape[1] / 2])

        # we got linear data, we want log_e
        # todo: implement check

        tmp_data[:, 0:tmp_data.shape[1] / 2] = np.log(
            tmp_data[:, 0:tmp_data.shape[1] / 2])

        self.data = tmp_data

    def load_frequencies(self, filename, ignore):
        """
        Load frequencies from file.
        """
        self.frequencies = self.__load_file(filename)
        if(ignore is not None):
            ignore_ids = np.array(ignore, dtype=int)
            self.frequencies = np.delete(self.frequencies, ignore_ids)
        self.fin = np.hstack((self.frequencies, self.frequencies))

    def __load_file(self, filename):
        """
        Wrapper for np.loadtxt with exception handling
        """
        try:
            data = np.loadtxt(filename)
        except IOError:
            print('There was an error opening the file: {0}'.format(filename))
            data = np.array([])
        return data

    def set_nr_cc_terms(self, nr):
        """
        Set the number of Cole-Cole terms to fit to the data and allocate
        corresponding arrays.
        """
        self.nr_cc_terms = nr
        new_nr_of_pars = 1 + 3 * nr
        datasets, dummy = self.data.shape

        if self.cc_pars is not None:
            self.cc_pars = np.resize(self.cc_pars, (datasets, new_nr_of_pars))
            self.cc_pars_init = np.resize(self.cc_pars_init,
                                          (datasets, new_nr_of_pars))
        else:
            self.cc_pars = np.zeros((datasets, new_nr_of_pars))
            self.cc_pars_init = np.zeros((datasets, new_nr_of_pars))

        if(self.errors is not None):
            self.errors = np.resize(self.errors, (datasets, new_nr_of_pars))
        else:
            self.errors = np.zeros((datasets, new_nr_of_pars))

    def initial_parameter_heuristic_1(self, spectrum):
        """
        Heuristic 1

        Initialize a Cole-Cole parameter set with heuristically determined
        parameters.

        This heuristic uses magic numbers to set the initial values
        """
        nr_cc_pars = 1 + 3 * self.nr_cc_terms
        p0 = np.zeros((1, nr_cc_pars)).flatten()

        # empirically determine an m value
        rho0_init = spectrum[0]
        m_init = 0.1
        # m_init = -0.0031604504 * np.mean(spectrum[spectrum.shape[0]/2:])

        tmparray = np.array((rho0_init, m_init, np.log(0.01), 0.6)).flatten()
        p0[0:4] = tmparray
        for additional_cc_terms in range(0, self.nr_cc_terms - 1):
            p0[4 + additional_cc_terms * 3: 7 + additional_cc_terms * 3] = \
                self.get_values_for_secondary_terms()
        return p0

    def get_values_for_secondary_terms(self):
        """
        Magic numbers for secondary and upwards Cole-Cole terms.
        """
        return [0.2, -10.127, 0.6]

    def initial_parameter_heuristic_2(self, spectrum):
        """
        Heuristic 2

        Initialize a Cole-Cole parameter set with heuristically determined
        parameters.

        This heuristic uses the following procedure for the CC parameters:

            - rho0 is determined by the lowest frequency magnitude
            - m is determined using a simple line search which samples the
              space 0.01 - 1.0
            - tau values are determined by logspacing the number of terms over
              the frequency range.
            - c is set to 0.5
        """
        print('Heuristic 2')
        nr_cc_pars = 1 + 3 * self.nr_cc_terms
        p0 = np.zeros((1, nr_cc_pars)).flatten()
        # rho0
        p0[0] = spectrum[0]

        # c
        c_indices = [i * 3 + 3 for i in xrange(self.nr_cc_terms)]
        p0[c_indices] = 0.5
        # tau
        tau_min = 1.0 / (2 * np.pi * self.frequencies.max())
        tau_max = 1.0 / (2 * np.pi * self.frequencies.min())

        tau_indices = [i * 3 + 2 for i in xrange(self.nr_cc_terms)]
        # log-space sampling of tau range
        # use 2 more tau values for the boundaries:
        # we do not want to set an initial tau value to one the frequency
        # boundaries.
        tau_space = np.logspace(np.log10(tau_min),
                                np.log10(tau_max),
                                self.nr_cc_terms + 2)
        p0[tau_indices] = np.log(tau_space[1:-1])

        # m
        # m is selected by testing m values within a certain value range for
        # the smallest phase rms
        m_indices = [i * 3 + 1 for i in xrange(self.nr_cc_terms)]
        test_m_values = np.linspace(0.05, 0.9, 10)
        best_index = -1
        best_rms = np.inf

        for index, test_m in enumerate(test_m_values):
            p0[m_indices] = test_m
            test_forward = colecole.cole_log(self.fin, p0)
            test_phase_rms = self._phase_rms(spectrum, test_forward)
            if test_phase_rms < best_rms:
                best_index = index
                best_rms = test_phase_rms

        p0[m_indices] = test_m_values[best_index]

        return p0

    def set_initial_values(self, heuristic_number, options=None):
        """
        Set the initial parameters, using the selected heuristic
        """
        self.heuristic_number = heuristic_number

        # we set the starting parameters for each parameter set separately
        for index in range(0, self.data.shape[0]):
            data = self.data[index]
            if heuristic_number == 1:
                init_func = self.initial_parameter_heuristic_1
            elif heuristic_number == 2:
                init_func = self.initial_parameter_heuristic_2

            self.cc_pars_init[index] = init_func(data)

    def _phase_rms(self, spectrum, forward):
        pha_rms = np.sqrt(
            sum((spectrum[len(spectrum) / 2:len(spectrum)] -
                 forward[1, :]) ** 2))
        pha_rms /= (len(spectrum) / 2)
        return pha_rms

    def fit_spectrum(self, spectrum, p0):
        """
        Fit a spectrum, given the data, the frequencies, and the starting
        parameters.

        Return the fit results, the magnitude and phase RMS, and the forward
        response of the fit parameters.
        """
        nr_to_fit = len(spectrum) / 2

        y_meas = np.hstack(
            (spectrum[0:nr_to_fit],
             spectrum[len(spectrum) / 2:len(spectrum) / 2 + nr_to_fit]))
        x = np.hstack((self.fin[0:nr_to_fit],
                       self.fin[len(self.fin) / 2:len(self.fin) /
                                2 + nr_to_fit]))

        # the actual fitting routine
        plsq, cov, info, mesg, success = leastsq(
            self.__residuals, p0, args=(y_meas, x),
            full_output=1,
            maxfev=10000000)
        # plsq,cov,info,mesg,success = leastsq(
        #   self.__residuals, p0, Dfun=self.__Dres,
        #   args=(y_meas, x), maxfev= 10000000, full_output=1)

        # compute spectral response using fit results
        forward = colecole.cole_log(self.fin, plsq)

        # compute mag rms
        mag_rms = np.sqrt(
            sum((spectrum[0:len(spectrum) / 2] - forward[0, :]) ** 2))
        mag_rms /= (len(spectrum / 2))

        # compute pha rms
        pha_rms = self._phase_rms(spectrum, forward)

        # compute complex rms
        rms = np.sqrt(
            sum((spectrum[0:len(spectrum) / 2] - forward[0, :]) ** 2) +
            sum((spectrum[len(spectrum) / 2:len(spectrum)] -
                 forward[1, :]) ** 2))
        rms /= len(spectrum)

        # calculate final chi square
        chisq = sum(info["fvec"] * info["fvec"])

        # deegrees of freedom
        dof = len(self.fin) / 2 - len(p0)

        # residuals d - f
        # mx_fvec = self.__residuals(plsq, y_meas, x)

        # compute fit errors ('aymptotic standard error')
        np.seterr(all='raise')
        if cov is not None:
            errors = np.zeros((plsq.shape[0]))
            for i, pmin in enumerate(plsq):
                try:
                    # uncertainties are calculated as per gnuplot, "fixing" the
                    # result -                # for non unit values of the
                    # reduced chisq.
                    errors[i] = np.sqrt(cov[i, i]) * np.sqrt(chisq / dof)
                except FloatingPointError:
                    print('FloatingPointError:')
                    print(cov[i, i])
                    print(chisq / dof, chisq, dof)
                    errors[i] = np.nan

        else:
            errors = None

        # verbose output of fit results if option is selected
        plot = False
        if plot:
            if cov is not None:
                print("Fitted parameters at minimum, with 68% C.I.:")
                for i, pmin in enumerate(plsq):
                    print('{0:02} {1:.10} +/- {2:.10} ({3:.10}%)'.format(
                        i, pmin, errors[i], errors[i] / pmin * 100))
                print
                print "Correlation matrix"
                for i in range(len(plsq)):
                    for j in range(i + 1):
                        print('{0}'.format(
                            cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])))
                    print('')
                print('')

        # return values
        fit_parameters = plsq

        # TODO: Sort CC terms for descending Tau values
        #       sort fit_parameters, mag_rms, pha_rms, forward, errors

        return fit_parameters, mag_rms, pha_rms, forward, errors

    def fit_all_spectra(self):
        # resize arrays
        self.magnitude_rms = np.zeros(self.data.shape[0])
        self.phase_rms = np.zeros(self.data.shape[0])
        self.forward_response = np.zeros(self.data.shape)
        self.residuals = np.zeros(self.data.shape)

        for id in range(0, self.data.shape[0]):
            fit_results, mag_rms, pha_rms, forward_response, errors = \
                self.fit_spectrum(
                    self.data[id, :], self.cc_pars_init[id, :])
            self.cc_pars[id, :] = fit_results
            self.errors[id, :] = errors
            self.magnitude_rms[id] = mag_rms
            self.phase_rms[id] = pha_rms
            self.forward_response[id] = forward_response.flatten()
            self.residuals[id] = self.data[id] - self.forward_response[id]

    def __Dres(self, p, y, x):
        err = 1
        # rho0
        if(np.abs(np.exp(p[0]) - np.exp(y[0])) > 0.5):
            err *= 1e10

        # m
        if(p[1] < 0 or (len(p) > 4 and p[4] < 0)):
            err *= 1e10

        if(p[1] > 1 or (len(p) > 4 and p[4] > 1)):
            err *= 1e10

        # tau
        if(p[2] < -12):
            err += 1e10
        # c
        if(p[3] < 0 or p[3] > 1):
            err *= 1e10

        # 2 CC terms and tau
        if(len(p) > 4 and (p[6] < 0 or p[6] > 1)):
            err *= 1e10
        ccJac = -colecole.cc_jac(x, p) * err
        return ccJac

    def __residuals(self, p, y, x):
        """
        Compute residuals of the measured data y and the response of the
        Cole-Cole function cole_log for the x(frequency) values and the
        parameters p. Used for leastsq function.
        """

        try:
            erg = colecole.cole_log(x, p)
        except:
            raise('error')

        erg1 = np.hstack((erg[0, :], erg[1, :]))
        err = y - erg1

        # NaN check and set err to 1e10
        if(np.isnan(err).any()):
            err *= 10e10

        # y holds the original data
        # compute weighting factors according to the invers of the phase values
        # mag_ori = y[:y.shape[0] / 2]
        # pha_ori = y[y.shape[0] / 2:]

        # mean_all = np.mean(np.abs(y))

        #   mean_weight_all = np.abs(y) / mean_all

        # err /= mean_weight_all

        # we want both mag and pha to have equal mean values AFTER weighting
        # pha_mean_weight = np.abs(np.mean(mag_ori)) / np.abs(np.mean(pha_ori))
#       print('mean(mag)', np.mean(mag_ori))
#       print('mean(pha)', np.mean(pha_ori * pha_mean_weight))

#       pha_weighting = np.abs(pha_ori)**2
#       #       pha_weighting = (np.abs(pha_ori)**4)
#        pha_weighting = (np.abs(pha_ori))**2 * 10
#     #  print(pha_weighting)
#       err[err.shape[0]/2:] *= pha_weighting
#       err[err.shape[0]/2:] *= pha_mean_weight

#       # mag
#       err[0:err.shape[0]/2] *= 100

        # implement a simple form of limit checking. If the given Cole-Cole
        # parameters lie out of certain bounds, increase the residual to a very
        # large value. This will ensure that this parameter set will not be
        # used any more.

#       #rho0
#       if(p[0] > 10e6 or np.abs(np.exp(p[0]) - np.exp(y[0])) > 0.5):
#           err *= 1e10

        # m
        if(p[1] < 0 or (len(p) > 4 and p[4] < 0)):
            err *= 1e10

        if(p[1] > 1 or (len(p) > 4 and p[4] > 1)):
            err *= 1e10

        # tau
        if(p[2] < -12):
            err += 1e10
        # c
        if(p[3] < 0 or p[3] > 1):
            err *= 1e10

        # 2 CC terms and tau
        if(len(p) > 4 and (p[6] < 0 or p[6] > 1)):
            err *= 1e10

        return err

    def save_mean_of_cc_pars(self, filename):
        self.cc_pars_mean = np.mean(self.cc_pars, axis=0)[np.newaxis, :]
        # compute mean error

        self.cc_error_mean = np.zeros(self.errors.shape[1])

        for i in self.errors:
            self.cc_error_mean += i ** 2

        self.cc_error_mean = np.sqrt(self.cc_error_mean)

        self.cc_error_mean /= self.errors.shape[0]

        try:
            np.savetxt(filename, self.cc_pars_mean, fmt='%.5f')
        except:
            print('There was an error saving mean Cole-Cole parameters ' +
                  'to the file: {0}'.format(filename))

        print('MEAN ERROR', self.cc_error_mean)
        try:
            np.savetxt(filename + '.err', self.cc_error_mean[np.newaxis, :],
                       fmt='%.5f')
        except:
            print('There was an error saving mean Cole-Cole parameter ' +
                  'errors to the file: {0}'.format(filename + '.err'))

    def save_cc_pars(self, filename):
        try:
            np.savetxt(filename, self.cc_pars, fmt='%.5f')
        except:
            pass

    def save_cc_errors(self, filename):
        try:
            np.savetxt(filename, self.errors, fmt='%.5f')
        except:
            pass

    def save_non_essential_results(self, directory):
        """Save non-essential results, such as results and forward resonse
        """
        np.savetxt(directory + os.sep + 'forward_response.dat',
                   self.forward_response)

        np.savetxt(directory + os.sep + 'residuals.dat',
                   self.residuals)

    def plot_all_spectra(self, directory=None, prefix=None):
        if prefix is None:
            prefix = '_'
        for id in range(0, self.data.shape[0]):
            self.plot_spectrum(
                directory + os.sep + '{1}spectrum_{0:02}'.format(id + 1,
                                                                 prefix), id)

    def plot_spectrum(self, filename, id):
        """
        Plot the spectrum specified by id to a file.
        """
        print('Plotting spectrum {0} of {1}'.format(id + 1,
                                                    self.data.shape[0]))
        # use more frequencies
        f_e = np.logspace(np.log10(self.frequencies.min()),
                          np.log10(self.frequencies.max()),
                          100)
        fin_e = np.hstack((f_e, f_e))
        # generate forward response
        forward_orig_f = colecole.cole_log(self.fin, self.cc_pars[id])
        forward = colecole.cole_log(fin_e, self.cc_pars[id])
        forward_init = colecole.cole_log(fin_e, self.cc_pars_init[id])

        # generate forward response for each CC term
        forward_cc_terms = []
        for i in range(0, (len(self.cc_pars[id]) - 1) / 3):
            oneterm_cc = [0, (i * 3) + 1, (i * 3) + 2, (i * 3) + 3]
            forward_cc_terms.append(
                colecole.cole_log(fin_e, self.cc_pars[id][oneterm_cc]))

        fig, axes = plt.subplots(3, 1, figsize=(7, 6))

        # plot magnitude
        ax = axes[0]
        ax.semilogx(f_e,
                    (np.exp(forward_init[0, :])),
                    'b-',
                    linestyle='dashed',
                    label='initial parameters')
        ax.semilogx(f_e,
                    (np.exp(forward[0, :])),
                    'g-',
                    label='fit response')

        ax.semilogx(self.frequencies,
                    (np.exp(self.data[id, 0: len(self.data[id, :]) / 2])),
                    'r.', linewidth=2.0,
                    label='data')
        ax.set_title('magnitude RMS: {0:.3e}'.format(self.magnitude_rms[id]))
        ax.set_xlabel('frequency [Hz]')
        ax.set_ylabel(r'$|Z/\rho| [\Omega (m)]$')

        # plot phase
        ax = axes[1]
        ax.semilogx(f_e, -forward[1, :], 'g-', linewidth=2.0,
                    label='fit response')
        ax.semilogx(self.frequencies,
                    -self.data[id, len(self.data[id, :]) / 2:],
                    'r.', linewidth=2.0,
                    label='data')

        ax.semilogx(f_e,
                    -forward_init[1, :],
                    'b-', linestyle='dashed',
                    label='initial model')
        # plot single terms
        colors = ('k', 'gray')
        for nr, term in enumerate(forward_cc_terms):
            ax.semilogx(f_e,
                        -term[1, :],
                        '-', color=colors[nr % 2],
                        linestyle='dashed',
                        label='term {0}'.format(nr + 1)
                        )

        ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0, 0, 1, 1),
                  bbox_transform=fig.transFigure)
        ax.set_title('phase RMS: {0:.3e}'.format(self.phase_rms[id]))
        ax.set_xlabel('frequency [Hz]')
        ax.set_ylabel(r'$-\varphi~[mrad]$')

        # plot phase residuals
        ax = axes[2]
        pha_residuals = np.abs(forward_orig_f[1, :] -
                               self.data[id, len(self.data[id, :]) / 2:])
        ax.semilogx(self.frequencies, pha_residuals, '-')
        ax.set_xlabel('frequency [Hz]')
        ax.set_ylabel(r'abs ($\phi_{inv} - \phi_{ori}$)')
        ax.set_title(
            'RMS: {0}'.format(np.sqrt(np.sum(pha_residuals ** 2)) /
                              self.frequencies.shape[0]))

        for ax in axes:
            ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        fig.tight_layout()
        fig.subplots_adjust(top=0.8)
        fig.savefig('{0}.png'.format(filename))
        plt.close(fig)

    def check_and_correct_cc_parameter_set(p0):
        """
        Check the given Cole-Cole parameter set for certain, heuristically
        determined, parameter bounds. This function is usefull when the one fit
        result is used as the starting value of another.
        """
        nr_cc_pars = p0.size

        # if the m values are not between 0 and 1, set them to 0.01
        for index in range(1, nr_cc_pars, 3):
            if(p0[index] > 1 or p0[index] < 0):
                p0[index] = 0.01

        # if the log(tau) values are not between -1e12 and 6, reset them
        counter = 0
        for index in range(2, nr_cc_pars, 3):
            if(p0[index] > 6 or p0[index] < -12):
                p0[index] = np.log(0.01 / 10 ** (counter))
            counter += 1

        # reset c if outside of [0, 1]
        counter = 0
        for index in range(3, nr_cc_pars, 3):
            if p0[index] > 1 or p0[index] < 0:
                p0[index] = 0.4 + 0.2 * counter
            counter += 1

        return p0
