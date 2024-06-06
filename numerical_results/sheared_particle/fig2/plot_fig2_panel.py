#!/usr/bin/env python3

import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from argparse import ArgumentParser
from matplotlib import colors, colormaps

# set matplotlib global font settings
# (change depending on font and LaTeX availability)
mpl.rc('font', **{'family': 'serif', 'serif': ['Times']})
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')

# set matplotlib global axis linewidths thicker than default
axis_linewidth = 1.5
mpl.rcParams['axes.linewidth'] = axis_linewidth
mpl.rcParams['xtick.major.width'] = axis_linewidth
mpl.rcParams['xtick.minor.width'] = axis_linewidth
mpl.rcParams['ytick.major.width'] = axis_linewidth
mpl.rcParams['ytick.minor.width'] = axis_linewidth

# key figure dimensions in inches
# ("plot" refers to the area spanned by the 
#  figure axis and axis labels, and does not
#  include the space occupied by the title and legend)
plot_width = 1.5
plot_height = 1.5
colorbar_width = 0.75
title_height = 0.25

# these adjust fractions assume figured is comprised
# by just the plot (defined above), and are corrected
# below to account for title and legend dimensions
plot_top_adjust = 0.95
plot_left_adjust = 0.30
plot_right_adjust = 0.95
plot_bottom_adjust = 0.30

# fontsizes (pt)
ticklabel_fontsize = 10
axislabel_fontsize = 11


def read_theta_from_simulation_output(filename):
    '''Reads policy parameter vector in the format
       output by its optimizer'''

    theta_list = []

    with open(filename, 'r') as file:
        for line in file:
            this_line = line.strip().split()
            theta_values = [float(v) for v in this_line
                            if v != '\\']
            theta_list.extend(theta_values)
            if this_line[-1] == '\\':
                break

    return np.asarray(theta_list, dtype=np.float64)


def squareexp_xy(nx, ny, sigma, x, y):
    '''Evaluates the squared exponential control policy basis on the
       2d meshgrid spanned by x and y, assigning the spacing between
       the basis functions from the scale parameter sigma'''

    if (x.ndim != 1) or (y.ndim != 1):
        raise ValueError(u'Input ``x`` and ``y`` must be 1d arrays')

    if (sigma <= 0.0):
        raise ValueError(u'Input ``sigma`` must be positive scalar')

    X, Y = np.meshgrid(x, y)

    σx = sigma / nx
    μx = σx * (np.arange(nx, dtype=np.float64) - (nx - 1) * 0.5)

    σy = sigma / ny
    μy = σy * (np.arange(ny, dtype=np.float64) - (ny - 1) * 0.5)

    U, V = np.meshgrid(μx, μy)

    XX = (U[:, :, None, None] - X[None, None, :, :]) / σx
    YY = (V[:, :, None, None] - Y[None, None, :, :]) / σy

    fxy = np.exp(-0.5 * (XX ** 2 + YY ** 2))
    dfxydx = fxy * XX / σx
    dfxydy = fxy * YY / σy

    return (dfxydx, dfxydy, fxy,)


def eval_doob(x, y, theta, sigma, grad=False):
    '''Evaluates the squared exponential control policy ansatz on the
       2d meshgrid spanned by x and y at control parameter vector theta
       and basis scale parameter sigma'''

    if (theta.ndim != 2):
        raise ValueError(u'Input ``theta`` must be a 2d array')

    (nx, ny,) = theta.shape
    (dfxydx, dfxydy, fxy,) = squareexp_xy(nx, ny, sigma, x, y)

    def contract_basis(xy_basis):
        return np.sum(theta[:, :, None, None] * xy_basis, axis=(0, 1,))

    if grad:
        xy_basis = (dfxydx, dfxydy, -fxy,)
    else:
        xy_basis = (-fxy,)

    return tuple(map(contract_basis, xy_basis))


def eval_pes(x, y):
    '''Plots the potential energy surface that confines the
       2d sheared particle'''

    if (x.ndim != 1) or (y.ndim != 1):
        raise ValueError(u'Input ``x`` and ``y`` must be 1d arrays')

    X, Y = np.meshgrid(x, y)
    pes = 2.0 * (6.0 + 4.0 * X ** 4 - 6.0 * Y ** 2 +
                 3.0 * Y ** 4 + 10.0 * X ** 2 * (Y ** 2 - 1.0))

    return pes


def eval_neqpes(x, y, theta, sigma):
    '''Plots the effective nonequilibrium potential associated with
       the 2d sheared particle in a given NESS'''

    return (eval_pes(x, y) + eval_doob(x, y, theta, sigma)[-1],)


# colormap for effective equilibrium potential
neqpes_cmap = colormaps['turbo_r']
neqpes_norm = colors.Normalize(vmin=0.0, vmax=10.0)
neqpes_clim = (0.0, 10.0,)

# colormap for control potential
doob_cmap = colormaps['magma']
doob_norm = colors.Normalize(vmin=0.0, vmax=4.0)
doob_clim = (0.0, 4.0,)

# simple argument parser
parser = ArgumentParser()
parser.add_argument('pickle', type=str,)
parser.add_argument('-plot_neqpes', action='store_true')

 # run this script in its current directory as, e.g.,:
 # ```python3 plot_fig2_panel.py './fig2a_data.pickle' -plot_neqpes```
if __name__ == '__main__':

    args = parser.parse_args()

    # load pickle file contents
    data_dict = pickle.load(open(args.pickle, 'rb'))

    # control policy ansatz hyper-parameters
    nx, ny, = data_dict['basis_shape']  # number of basis function along x and y axes
    sigma = data_dict['basis_scale']   # scale parameter for all basis functions

    # control policy parameter vector (amplitude coefficients for basis functions)
    theta, = data_dict['theta']
    theta.shape = (nx, ny,)

    # create figure with set panel dimensions
    fig_height = plot_height + title_height
    fig_width = plot_width + colorbar_width
    fig, ax, = plt.subplots(nrows=1, ncols=1, figsize=(fig_width, fig_height,), layout='none')

    # correct adjust fractions to take into account title and colorbar dimensions
    top_adjust = 1.0 - ((1.0 - plot_top_adjust) * plot_height + title_height) / fig_height
    left_adjust = plot_left_adjust * plot_width / fig_width
    right_adjust = 1.0 - ((1.0 - plot_right_adjust) * plot_width + colorbar_width) / fig_width
    bottom_adjust = plot_bottom_adjust * plot_height / fig_height

    fig.subplots_adjust(top=top_adjust,
                        left=left_adjust,
                        right=right_adjust,
                        bottom=bottom_adjust,)

    # set axis (tick)labels and figure title
    ax.set_xticks([-1.2, 0.0, 1.2,])
    ax.set_yticks([-1.2, 0.0, 1.2,])

    ax.xaxis.set_tick_params(labelsize=ticklabel_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticklabel_fontsize)

    ax.set_xlabel(r'$x_1$', fontsize=axislabel_fontsize, labelpad=0.0)
    ax.set_ylabel(r'$x_2$', fontsize=axislabel_fontsize, labelpad=-3.5)
    ax.set_title(data_dict['panel_title'], loc='center', pad=4.0, fontsize=axislabel_fontsize)

    # choose colormap depending on whether plotting the nonequilibrium pes
    cmap = neqpes_cmap if args.plot_neqpes else doob_cmap
    norm = neqpes_norm if args.plot_neqpes else doob_norm
    clim = neqpes_clim if args.plot_neqpes else doob_clim

    # render the control potential or the nonequilibrium pes on a square grid
    grid_lims = (-1.5, 1.5,)
    grid_size = 401
    x = y = (max(grid_lims) - min(grid_lims)) * \
        (np.arange(grid_size) - (grid_size - 1) * 0.5) / grid_size

    if args.plot_neqpes:
        heatmap = eval_neqpes(x, y, theta, sigma)[-1]
    else:
        heatmap = eval_doob(x, y, theta, sigma)[-1]

    im = ax.imshow(heatmap - heatmap.min(), extent=2 * grid_lims,
                   interpolation='none', origin='lower', aspect='equal',
                   norm=norm, cmap=cmap,)

    # fit a colorbar into the remaining space of the panel
    left, bottom, width, height = plt.gca().get_position().bounds
    cax = fig.add_axes((0.2 * colorbar_width / fig_width + left + width, bottom,
                        0.2 * colorbar_width / fig_width, height,))
    cb = fig.colorbar(im, cax=cax, format=lambda x, _: f"{x:2.0f}")

    # choose colorbar label depending on whether plotting the nonequilibrium pes
    if args.plot_pdf:
        cb_label = r'$\beta \bigl(U + V_{\boldsymbol{\theta}}\bigr)(X)$'
    else:
        cb_label = r'$\beta V_{\boldsymbol{\theta}}(X)$'

    cb.set_label(cb_label, size=ticklabel_fontsize, labelpad=1.0)
    cb.ax.set_ylim(clim)
    cb.set_ticks(np.linspace(*cb.ax.get_ylim(), 3).tolist())
    cb.ax.tick_params(labelsize=ticklabel_fontsize)

    plt.show()
