#!/usr/bin/env python3

import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from argparse import ArgumentParser
from matplotlib import colors, colormaps

# NOTE: makes use of `brokenaxes` library; see
# `https://github.com/bendichter/brokenaxes`
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec

# set matplotlib global font settings
# (change depending on font and LaTeX availability)
mpl.rc('font', **{'family': 'serif',  'serif': ['Times']})
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
#  include the space occupied by the title and legend)
plot_width = 1.2
plot_height = 1.2
legend_width = 1.0
title_height = 0.2

# these adjust fractions assume the figure is comprised
# by just the plot (defined above), and are corrected
# below to account for title and legend dimensions
plot_top_adjust = 0.9
plot_left_adjust = 0.3
plot_right_adjust = 0.9
plot_bottom_adjust = 0.3

# fontsizes (pt)
legend_fontsize = 10
ticklabel_fontsize = 10
axislabel_fontsize = 11

mpl.rcParams["xtick.labelsize"] = ticklabel_fontsize
mpl.rcParams["ytick.labelsize"] = ticklabel_fontsize

parser = ArgumentParser(description='''Plot optimized basis parameters for 
                                       swimming N-mer at various Péclet''')

parser.add_argument('pickle', type=str)
parser.add_argument('-theta_label_prefix', type=str, default=r'$\mathrm{Pe}')

parser.add_argument('-k_bond_ylim', type=float, nargs=2, default=[-4, 4])      # default=[None, None,])
parser.add_argument('-l_bond_ylim', type=float, nargs=2, default=[-3, 0.5])    # default=[None, None,])
parser.add_argument('-k_bend_ylim', type=float, nargs=2, default=[-1, 1])      # default=[None, None,])
parser.add_argument('-l_bend_ylim', type=float, nargs=2, default=[-0.6, 0.1])  # default=[None, None,])

parser.add_argument('-k_bond_ref', type=float, default=1000.0)
parser.add_argument('-l_bond_ref', type=float, default=1.0)
parser.add_argument('-k_bend_ref', type=float, default=10.0)
parser.add_argument('-l_bend_ref', type=float, default=0.0)

if __name__ == '__main__':

    args = parser.parse_args()

    # load pickle file contents
    data_dict = pickle.load(open(args.pickle, 'rb'))

    theta_labels = sorted([key for key in data_dict.keys()
                           if key.startswith(args.theta_label_prefix)])
    theta_data = {label: data_dict[label] for label in theta_labels}

    N = int(data_dict['N'])  # args.monomer_number
    bond_indices = np.arange(N - 1) + 1
    bend_indices = np.arange(N - 2) + 1

    fig_height = plot_height
    fig_width = 4 * plot_width + legend_width

    top_adjust = 1.0 - ((1.0 - plot_top_adjust) * plot_height) / fig_height
    left_adjust = plot_left_adjust * plot_width / fig_width
    right_adjust = 1.0 - ((1.0 - plot_right_adjust) * plot_width + legend_width) / fig_width
    bottom_adjust = plot_bottom_adjust
    wspace_adjust = (plot_left_adjust + 1.0 - plot_right_adjust) / (plot_right_adjust - plot_left_adjust)
    hspace_adjust = 0.0

    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    gs = GridSpec(nrows=1, ncols=4, figure=fig,
                  top=top_adjust,
                  left=left_adjust,
                  right=right_adjust,
                  bottom=bottom_adjust,
                  wspace=wspace_adjust,
                  hspace=hspace_adjust,
                  width_ratios=None, height_ratios=None)

    ax_k_bond_xlims = ((1, 5), (N-5, N-1))
    ax_l_bond_xlims = ((1, 5), (N-5, N-1))
    ax_k_bend_xlims = ((1, 4), (N-5, N-2))
    ax_l_bend_xlims = ((1, 4), (N-5, N-2))

    brokenaxes_kwargs = dict(tilt=45, d=0.008, wspace=0.6, despine=False)
    plot_kwargs = dict(linestyle='-', marker='o', markerfacecolor='white',
                       linewidth=1.5, markersize=3.0,)

    ax_k_bond = brokenaxes(xlims=ax_k_bond_xlims, subplot_spec=gs[0], **brokenaxes_kwargs)
    ax_l_bond = brokenaxes(xlims=ax_l_bond_xlims, subplot_spec=gs[1], **brokenaxes_kwargs)
    ax_k_bend = brokenaxes(xlims=ax_k_bend_xlims, subplot_spec=gs[2], **brokenaxes_kwargs)
    ax_l_bend = brokenaxes(xlims=ax_l_bend_xlims, subplot_spec=gs[3], **brokenaxes_kwargs)

    label_colornorm = colors.Normalize(vmin=0.0, vmax=1.0)
    label_colormap = lambda c: colormaps['winter'](label_colornorm(c))

    n_label = len(theta_labels)
    label_color = {label: label_colormap(i / n_label)
                   for (i, label,) in enumerate(theta_labels, start=1)}

    for label in theta_labels:
        data = theta_data[label]  # read_theta_from_simulation_output(theta_files[label])

        k_bond_data = data[0 * N + 1 : 1 * N]
        l_bond_data = data[1 * N + 1 : 2 * N]
        k_bend_data = data[2 * N + 2 : 3 * N]
        l_bend_data = data[3 * N + 2 : 4 * N]

        # k_bond_data = (k_bond_data - args.k_bond_ref) / args.k_bond_ref
        # l_bond_data = (l_bond_data - args.l_bond_ref)  # / args.l_bond_ref
        # k_bend_data = (k_bend_data - args.k_bend_ref) / args.k_bend_ref
        # l_bend_data = (l_bend_data - args.l_bend_ref)  # / args.l_bend_ref

        ax_k_bond.plot(bond_indices, k_bond_data, color=label_color[label], **plot_kwargs,)
        ax_l_bond.plot(bond_indices, l_bond_data, color=label_color[label], **plot_kwargs,)
        ax_k_bend.plot(bend_indices, k_bend_data, color=label_color[label], **plot_kwargs,)
        ax_l_bend.plot(bend_indices, l_bend_data, color=label_color[label], **plot_kwargs,
                       label=r'{:s}'.format(label.strip()))

    left, bottom, width, height = plt.gca().get_position().bounds
    legend_bbox = [left + width + 0.02, 0.9 * bottom, legend_width / fig_width, bottom]

    fig.legend(*ax_l_bend.axs[0].get_legend_handles_labels(),
               frameon=True, reverse=False,
               loc='lower left', bbox_to_anchor=legend_bbox,
               facecolor='white', edgecolor='white', framealpha=1.0,
               markerscale=1.0, labelspacing=0.2, borderpad=0.0, fontsize=legend_fontsize,
               handlelength=0.6, handleheight=0.5, handletextpad=0.4)

    ax_l_bond.axhline(y=0.0, c='k', lw=0.5, zorder=-10)
    ax_k_bond.axhline(y=0.0, c='k', lw=0.5, zorder=-10)
    ax_l_bend.axhline(y=0.0, c='k', lw=0.5, zorder=-10)
    ax_k_bend.axhline(y=0.0, c='k', lw=0.5, zorder=-10)

    ax_k_bond.set_xlabel(r'bond $i$', fontsize=ticklabel_fontsize, labelpad=16)
    ax_l_bond.set_xlabel(r'bond $i$', fontsize=ticklabel_fontsize, labelpad=16)
    ax_k_bond.set_ylabel(r'$k_i^{\mathrm{bond}}$', fontsize=axislabel_fontsize+1, labelpad=4)
    ax_l_bond.set_ylabel(r'$\ell_i^{\mathrm{bond}}$', fontsize=axislabel_fontsize+1, labelpad=4)

    ax_k_bend.set_xlabel(r'angle $i$', fontsize=ticklabel_fontsize, labelpad=16)
    ax_l_bend.set_xlabel(r'angle $i$', fontsize=ticklabel_fontsize, labelpad=16)
    ax_k_bend.set_ylabel(r'$k_i^{\mathrm{bend}}$', fontsize=axislabel_fontsize+1, labelpad=4)
    ax_l_bend.set_ylabel(r'$\ell_i^{\mathrm{bend}}$', fontsize=axislabel_fontsize+1, labelpad=4)

    ax_k_bond.big_ax.set_xlim(ax_k_bond_xlims[0][0] - 1, ax_k_bond_xlims[-1][-1] + 1)
    ax_l_bond.big_ax.set_xlim(ax_l_bond_xlims[0][0] - 1, ax_l_bond_xlims[-1][-1] + 1)
    ax_k_bend.big_ax.set_xlim(ax_k_bend_xlims[0][0] - 1, ax_k_bend_xlims[-1][-1] + 1)
    ax_l_bend.big_ax.set_xlim(ax_l_bend_xlims[0][0] - 1, ax_l_bend_xlims[-1][-1] + 1)

    ax_k_bond.axs[0].set_ylim(args.k_bond_ylim)
    ax_k_bond.axs[1].set_ylim(args.k_bond_ylim)
    ax_k_bond.axs[0].set_xticks(ax_k_bond_xlims[0])
    ax_k_bond.axs[1].set_xticks(ax_k_bond_xlims[1])

    ax_l_bond.axs[0].set_ylim(args.l_bond_ylim)
    ax_l_bond.axs[1].set_ylim(args.l_bond_ylim)
    ax_l_bond.axs[0].set_xticks(ax_l_bond_xlims[0])
    ax_l_bond.axs[1].set_xticks(ax_l_bond_xlims[1])

    ax_k_bend.axs[0].set_ylim(args.k_bend_ylim)
    ax_k_bend.axs[1].set_ylim(args.k_bend_ylim)
    ax_k_bend.axs[0].set_xticks(ax_k_bend_xlims[0])
    ax_k_bend.axs[1].set_xticks(ax_k_bend_xlims[1])

    ax_l_bend.axs[0].set_ylim(args.l_bend_ylim)
    ax_l_bend.axs[1].set_ylim(args.l_bend_ylim)
    ax_l_bend.axs[0].set_xticks(ax_l_bend_xlims[0])
    ax_l_bend.axs[1].set_xticks(ax_l_bend_xlims[1])

    ax_l_bond.axs[0].set_yticks([-3.0, 0.0])
    ax_k_bond.axs[0].set_yticks([-4.0, 4.0])
    ax_l_bend.axs[0].set_yticks([-0.6, 0.0])
    ax_k_bend.axs[0].set_yticks([-1.0, 1.0])

    plt.show()
