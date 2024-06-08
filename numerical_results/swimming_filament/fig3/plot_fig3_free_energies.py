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

# plot dimensions in inches
# ("plot" refers to the area spanned by the 
#  figure axis and axis labels, and does not
#  include the space occupied by the top title,
#  bottom legend, and leftmost axis label)
plot_height = 0.8
title_height = 0.2
legend_height = 0.3

plot_width = 1.8
label_width = 0.2

# these adjust fractions assume the figure is comprised
# by just the plot (defined above), and are corrected
# below to account for the dimensions of features
plot_top_adjust = 0.97
plot_left_adjust = 0.35
plot_right_adjust = 0.97
plot_bottom_adjust = 0.40

# fontsizes (pt)
legend_fontsize = 9
ticklabel_fontsize = 10
axislabel_fontsize = 11

# linewidth, markersize, and alpha
plot_linewidth = 2.0
plot_markersize = 4.0
fill_alpha = 0.15

feq_ref_color = '0.7'
fneq_ref_color = 'black'
fneq_Q_color = 'red'
fneq_T_color = 'darkturquoise'
fneq_V_color = 'blue'

feq_ref_label = r'$F^\mathrm{eq}$'
fneq_ref_label = r'$F^\mathrm{neq}$'
fneq_Q_label = r'$\hat{F}^\mathrm{neq}_{Q}$'
fneq_T_label = r'$\hat{F}^\mathrm{neq}_{T}$'
fneq_V_label = r'$\hat{F}^\mathrm{neq}_{\mathrm{VTR}}$'

# default kwargs for matplotlib.axes.Axes.plot()
plot_kwargs = dict(linestyle='-',
                   linewidth=plot_linewidth,
                   marker='.',
                   markerfacecolor='white',
                   markeredgewidth=0.5 * plot_linewidth,
                   markevery=1,)


def ax_plot_shifted_line(ax, x, y, **plot_kwargs):
    '''Thin wrapper for matplotlib.axes.Axes.plot
       that y-shifts xydata minimum to zero'''

    y_min = y[~np.isnan(y)].min()

    line, = ax.plot(x, y - y_min, **plot_kwargs)
    return line


def ax_plot_shifted_linefill(ax, x, y, yerror, alpha=0.0, fill_scale=1.0, **plot_kwargs):
    '''Plots y-shifted xy-data with surrounding areas filled to represent errorbars'''

    line = ax_plot_shifted_line(ax, x, y, **plot_kwargs)

    fill = ax.fill_between(line.get_xdata(),
                           line.get_ydata() - yerror.max() * abs(fill_scale),
                           line.get_ydata() + yerror.max() * abs(fill_scale),
                           facecolor=line.get_color(), edgecolor=line.get_color(),
                           linewidth=0.25 * line.get_linewidth(), alpha=alpha,
                           zorder=line.get_zorder()-0.5,)
    return (line, fill,)


parser = ArgumentParser(description='''Plot nonequilibrium free-energy surfaces
                                       for the N-mer swimming filament''')
parser.add_argument('pickle', type=str,
                    metavar='path-to-data-dict-pickle',
                    help='location of pickle containing FES estimates')
parser.add_argument('peclet', type=str,
                    metavar='swimming-speed',
                    help='swimming speed for which plot FES estimates (either 10, 20, 40, or 60)')

if __name__ == '__main__':

    args = parser.parse_args()

    # load pickle file contents
    data_dict = pickle.load(open(args.pickle, 'rb'))
    assert (args.peclet in data_dict['peclet_keys'])

    # numerically exact FES estimates
    feq_ref = data_dict['feq_ref'][args.peclet]    # equilibrium (from WHAM on umbrella-biased trajectories)
    fneq_ref = data_dict['fneq_ref'][args.peclet]  # nonequilibrium (from trajectory histogramming)

    # Kawasaki--Crooks nonequilibrium FES estimates
    fneq_Q = data_dict['fneq_Q'][args.peclet]  # from heat evaluation along nonequilibrium relaxation trajectories
    fneq_T = data_dict['fneq_T'][args.peclet]  # from traffic evaluation along equilibrium relaxation trajectories
    fneq_V = data_dict['fneq_V'][args.peclet]  # from variational time reversal of relaxation trajectories

    fig_width = plot_width + label_width
    fig_height = plot_height + title_height + legend_height

    fig, ax, = plt.subplots(ncols=1, nrows=1, layout='none',
                            figsize=(fig_width, fig_height,),)

    top_adjust = 1.0 - ((1.0 - plot_top_adjust) * plot_height + title_height) / fig_height
    left_adjust = plot_left_adjust * plot_width / fig_width
    right_adjust = 1.0 - ((1.0 - plot_right_adjust) * plot_width + label_width) / fig_width
    bottom_adjust = (plot_bottom_adjust * plot_height + legend_height) / fig_height

    fig.subplots_adjust(top=top_adjust,
                        left=left_adjust,
                        right=right_adjust,
                        bottom=bottom_adjust,)

    # set axis limits and figure title
    ax.set_xlim(data_dict['panel_xlim'])
    ax.set_ylim(data_dict['panel_ylim'])

    ax.set_xlabel(r'$r$', fontsize=axislabel_fontsize, labelpad=0.0)
    ax.set_ylabel(r'$F(r)$', fontsize=axislabel_fontsize, labelpad=4.0)
    ax.set_title(data_dict['panel_title'], loc='center', pad=3.0, fontsize=axislabel_fontsize)

    # set axis ticklabels
    xticks = np.linspace(*tuple(ax.get_xlim()), num=3, endpoint=True).tolist()
    ax.set_xticks(xticks)
    yticks = [0.0, 3.0, 6.0]
    ax.set_yticks(yticks)

    ax.xaxis.set_tick_params(labelsize=ticklabel_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticklabel_fontsize)

    # set figure label (appears to the right of the plot area)
    ax.text(1.1, 0.5, r'$\mathrm{{Pe}} = {:s}$'.format(args.peclet),
            fontsize=axislabel_fontsize, ha='left', va='center',
            transform=ax.transAxes, rotation=90)

    ax.axhline(y=0.0, c='k', lw=0.5, zorder=-10)

    # plot equilibrium and nonequilibrium ground-truth FES estimates
    ax_plot_shifted_line(
        ax, feq_ref['x'], feq_ref['y'], ms=0.0,
        color=feq_ref_color, label=feq_ref_label, **plot_kwargs
    )
    ax_plot_shifted_line(
        ax, fneq_ref['x'], fneq_ref['y'], ms=plot_markersize,
        color=fneq_ref_color, label=fneq_ref_label, **plot_kwargs
    )

    # plot Kawasaki--Crooks nonequilibrium FES estimates
    ax_plot_shifted_linefill(
        ax, fneq_Q['x'], fneq_Q['y'], fneq_Q['yerr'], alpha=fill_alpha, fill_scale=0.5,
        color=fneq_Q_color, label=fneq_Q_label, ms=0.0, **plot_kwargs
    )
    # ax_plot_shifted_linefill(
    #     ax, fneq_T['x'], fneq_T['y'], fneq_T['yerr'], alpha=fill_alpha, fill_scale=0.5,
    #     color=fneq_T_color, label=fneq_T_label, ms=0.0, **plot_kwargs
    # )
    ax_plot_shifted_linefill(
        ax, fneq_V['x'], fneq_V['y'], fneq_V['yerr'], alpha=fill_alpha, fill_scale=0.5,
        color=fneq_V_color, label=fneq_V_label, ms=0.0, **plot_kwargs
    )

    # draw legend at the bottom of plot area
    legend_bbox = [0.2 * left_adjust, 0.0, 1.0, 0.8 * bottom_adjust]
    fig.legend(frameon=True, reverse=False, ncols=4, columnspacing=0.5,
               loc='lower left', bbox_to_anchor=legend_bbox, borderpad=0.0,
               facecolor='white', edgecolor='white', framealpha=1.0,
               markerscale=1.0, labelspacing=0.2, fontsize=legend_fontsize,
               handlelength=0.6, handleheight=1.0, handletextpad=0.4)
    plt.show()
