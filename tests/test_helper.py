import matplotlib as mpl
import os
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil
import lib_cc_fit.colecole as colecole
import glob


def _generate_spectra(frequencies, testdir, cc_generator):
    if os.path.isdir(testdir):
        shutil.rmtree(testdir)
    os.makedirs(testdir)
    fin = np.hstack((frequencies, frequencies))

    for nr, cc_set in enumerate(cc_generator()):
        outdir = testdir + os.sep + '{0:03}_'.format(nr) + \
            '{0}_{1}_{2}_{3}'.format(*cc_set) + os.sep
        os.makedirs(outdir)
        # generate data
        rlogmag_rpha = colecole.cole_log(fin, cc_set)
        data = np.vstack((np.exp(rlogmag_rpha[0, :]),
                          rlogmag_rpha[1, :])).flatten()
        data_2d = np.atleast_2d(data)

        # save data
        np.savetxt(outdir + 'frequencies.dat', frequencies)
        np.savetxt(outdir + 'data.dat', data_2d)

        cc_set_2d = np.atleast_2d(cc_set)
        np.savetxt(outdir + 'colecole_orig.dat', cc_set_2d)


def _get_cc_dirs(testdir):
    cc_dirs = sorted(glob.glob(testdir + '/*'))
    return cc_dirs


def plot_residuals(directory):
    residuals = []
    for directory in _get_cc_dirs(directory):
        residual_file = directory + os.sep + 'residuals.dat'
        residuals.append(np.loadtxt(residual_file))
    all_residuals = np.array(residuals)
    nr_f = all_residuals.shape[1] / 2
    fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))

    ax = axes[0]
    im = ax.imshow(all_residuals[:, 0:nr_f], interpolation='none')
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    fig.add_axes(ax_cb)
    fig.colorbar(im, cax=ax_cb)
    ax.set_xlabel('frequency')
    ax.set_title(r'Magnitude $\Delta(log(R))$')

    ax = axes[1]
    im = ax.imshow(all_residuals[:, nr_f:], interpolation='none')
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    fig.add_axes(ax_cb)
    fig.colorbar(im, cax=ax_cb)
    ax.set_xlabel('frequency')
    ax.set_title(r'Phase $\Delta \phi$')
    fig.tight_layout()
    fig.savefig('residuals.png', dpi=300)


def plot_result_diffs(testdir):
    cc_fit_residuals = []
    for directory in _get_cc_dirs(testdir):
        cc_orig = directory + os.sep + 'colecole_orig.dat'
        cc_fit = directory + os.sep + 'cc_fits.dat'
        cc_fit_residuals.append(np.loadtxt(cc_fit) - np.loadtxt(cc_orig))
    all_cc_fit_residuals = np.array(cc_fit_residuals)

    labels = [r'$\Delta log(\rho|R)$',
              r'$\Delta m$',
              r'$\Delta log(\tau)$',
              r'$\Delta c$'
              ]

    fig, axes = plt.subplots(4, 1, figsize=(6, 5))

    for ax, parameter, label in zip(axes, all_cc_fit_residuals.T, labels):
        ax.plot(parameter, '.-')
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        ax.set_ylabel(label)
        # ax.set_xlabel('frequency')
        # ax.set_title(r'Magnitude $\Delta(log(R))$')

    fig.tight_layout()
    fig.savefig('cc_fit_residuals.png', dpi=300)


def _evaluate_fits(testdir):
    """ Create evaluation plots of fit results
    """
    plot_residuals(testdir)
    plot_result_diffs(testdir)
