#!/usr/bin/env python3

import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from argparse import ArgumentParser

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
plot_width = 1.6
plot_height = 1.4
legend_width = 0.75
title_height = 0.25

# these adjust fractions assume figured is comprised
# by just the plot (defined above), and are corrected
# below to account for title and legend dimensions
plot_top_adjust = 0.95
plot_left_adjust = 0.25
plot_right_adjust = 0.95
plot_bottom_adjust = 0.25

# fontsizes (pt)
legend_fontsize = 9
ticklabel_fontsize = 10
axislabel_fontsize = 11

# linewidth, markersize, and alpha
plot_linewidth = 1.5
plot_markersize = 4.0
fill_alpha = 0.25

feq_ref_color = '0.7'
fneq_ref_color = 'black'
fneq_Q_color = 'red'
fneq_T_color = 'darkturquoise'
fneq_V_color = 'blue'

feq_ref_label = r'$F^\mathrm{eq}$'
fneq_ref_label = r'$F^\mathrm{neq}$'
fneq_Q_label = r'$\hat{F}^\mathrm{neq}_{\mathcal{Q}}$'
fneq_T_label = r'$\hat{F}^\mathrm{neq}_{\mathcal{T}}$'
fneq_V_label = r'$\hat{F}^\mathrm{neq}_{\mathrm{VTR}}$'


plot_kwargs = dict(linestyle='-',
                   linewidth=plot_linewidth,
                   marker='.',
                   markerfacecolor='white',
                   markeredgewidth=0.5 * plot_linewidth,
                   markevery=2,)


def ax_plot_line(ax, x, y, yshift=0.0, **plot_kwargs):
    '''Thin wrapper for matplotlib.axes.Axes.plot
       allowing for shifts in y-axis baseline'''

    line, = ax.plot(x, y - yshift, **plot_kwargs)

    return line


def ax_plot_linefill(ax, x, y, yerror, yshift=0.0, alpha=1.0, **plot_kwargs):
    '''Plots line with surrounding y-axis filled areas representing errorbars
       allowing for shifts in y-axis baseline and arbitrary errorbar opacity'''

    line = ax_plot_line(ax, x, y, yshift=yshift, **plot_kwargs)

    fill = ax.fill_between(x, y - yshift - yerror, y - yshift + yerror,
                           facecolor=line.get_color(), edgecolor=line.get_color(),
                           linewidth=0.25 * line.get_linewidth(), alpha=alpha,
                           zorder=line.get_zorder()-0.5,)

    return (line, fill,)


hist_data_dirname = '/Users/jrosa/Dropbox/Remote_Code/filament_build/with_mpi/varying_dims/new_loss/deterministic_active_force/dvr_solver_sparse_output/ShearedDoubleWell2D/'
data_1d_dirname = '/Volumes/Scratch2/deterministic_active_force/cleaner_code/one-dimensional_basis/gaussian/zero_initialization/relax/'
data_2d_dirname = '/Volumes/Scratch2/deterministic_active_force/cleaner_code/two-dimensional_basis/gaussian/zero_initialization/relax/'

# simple argument parser
parser = ArgumentParser()
parser.add_argument('pickle', type=str,)

 # run this script in its current directory as, e.g.,:
 # ```python3 plot_fig1_panel.py './fig1a_data.pickle'```
if __name__ == '__main__':

    args = parser.parse_args()

    # load pickle file contents
    data_dict = pickle.load(open(args.pickle, 'rb'))

    # numerically exact FES estimates
    feq_ref = data_dict['feq_ref']    # equilibrium (from quadrature of exact density)
    fneq_ref = data_dict['fneq_ref']  # nonequilibrium (from trajectory histogramming)

    # Kawasaki--Crooks nonequilibrium FES estimates
    fneq_Q = data_dict['fneq_Q']  # from heat evaluation along nonequilibrium relaxation trajectories
    fneq_T = data_dict['fneq_T']  # from traffic evaluation along equilibrium relaxation trajectories
    fneq_V = data_dict['fneq_V']  # from variational time reversal of relaxation trajectories

    # create figure with set panel dimensions
    fig_height = plot_height + title_height
    fig_width = plot_width + legend_width
    fig, ax, = plt.subplots(nrows=1, ncols=1, figsize=(fig_width, fig_height,), layout='none')

    # correct adjust fractions to take into account title and legend dimensions
    top_adjust = 1.0 - ((1.0 - plot_top_adjust) * plot_height + title_height) / fig_height
    left_adjust = plot_left_adjust * plot_width / fig_width
    right_adjust = 1.0 - ((1.0 - plot_right_adjust) * plot_width + legend_width) / fig_width
    bottom_adjust = plot_bottom_adjust * plot_height / fig_height

    fig.subplots_adjust(top=top_adjust,
                        left=left_adjust,
                        right=right_adjust,
                        bottom=bottom_adjust,)

    # set axis limits, axis (tick)labels, and figure title
    ax.set_xlim([-1.8, 1.8,])
    ax.set_xticks([-1.5, 0.0, 1.5,])
    ax.set_ylim([-1.0, 7.0,])
    ax.set_yticks([0.0, 3.0, 6.0,])
    ax.axhline(y=0.0, c='black', lw=0.5, zorder=-10)

    ax.xaxis.set_tick_params(labelsize=ticklabel_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticklabel_fontsize)

    ax.set_xlabel(r'$r$', fontsize=axislabel_fontsize, labelpad=0.5)
    ax.set_ylabel(r'$F(r)$', fontsize=axislabel_fontsize, labelpad=2.0)
    ax.set_title(data_dict['panel_title'], loc='center', pad=3.0, fontsize=axislabel_fontsize)

    # plot numerically exact FES estimates
    feq_ref_line = ax_plot_line(
        ax, feq_ref['x'], feq_ref['y'], yshift=feq_ref['y'].min(),
        color=feq_ref_color, label=feq_ref_label, **plot_kwargs
    )
    fneq_ref_line = ax_plot_line(
        ax, fneq_ref['x'], fneq_ref['y'], yshift=fneq_ref['y'].min(), ms=plot_markersize,
        color=fneq_ref_color, label=fneq_ref_label, **plot_kwargs
    )

    # plot Kawasaki--Crooks nonequilibrium FES estimates
    ax_plot_linefill(
        ax, fneq_Q['x'], fneq_Q['y'], fneq_Q['yerr'], yshift=0.0, alpha=fill_alpha,
        color=fneq_Q_color, label=fneq_Q_label, **plot_kwargs
    )
    ax_plot_linefill(
        ax, fneq_T['x'], fneq_T['y'], fneq_T['yerr'], yshift=0.0, alpha=fill_alpha,
        color=fneq_T_color, label=fneq_T_label, **plot_kwargs
    )
    ax_plot_linefill(
        ax, fneq_V['x'], fneq_V['y'], fneq_V['yerr'], yshift=0.0, alpha=fill_alpha,
        color=fneq_V_color, label=fneq_V_label, **plot_kwargs
    )

    # draw legend to the right of plot area
    left, bottom, width, height = plt.gca().get_position().bounds
    bbox_to_anchor = [left + width, 0.0, legend_width / fig_width, 1.0]

    fig.legend(frameon=True, reverse=False,
               loc='center', bbox_to_anchor=bbox_to_anchor, borderpad=0.2,
               facecolor='white', edgecolor='white', framealpha=1.0,
               markerscale=1.0, labelspacing=0.2, fontsize=legend_fontsize,
               handlelength=0.4, handleheight=1.0, handletextpad=0.3)

    plt.show()
