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

# plot dimensions in inches
# ("plot" refers to the area spanned by the 
#  figure axis and axis labels)
plot_width = 1.5
plot_height = 1.33

# adjust fractions that locate figure axis within the plot
top_adjust = 0.97
left_adjust = 0.30
right_adjust = 0.97
bottom_adjust = 0.30

# fontsizes (pt)
legend_fontsize = 9
ticklabel_fontsize = 10
axislabel_fontsize = 11

# linewidth and markersize
plot_linewidth = 2.0
plot_markersize = 4.0

# peclet colorscale to generate colored version of this figure
peclet_min = 0.0
peclet_max = 60.0
peclet_cmap = colormaps['winter']
peclet_colorscale = colors.Normalize(vmin=peclet_min, vmax=peclet_max,)
peclet_cmap_lambda = lambda peclet: peclet_cmap(peclet_colorscale(np.float_(peclet)))

# default kwargs for matplotlib.axes.Axes.plot()
plot_kwargs = dict(linestyle='-', markerfacecolor='white',
                   linewidth=plot_linewidth, markeredgewidth=0.5 * plot_linewidth,)


def ax_plot_line(ax, x, y, yshift=0.0, **plot_kwargs):
    '''Thin wrapper for matplotlib.axes.Axes.plot
       allowing for shifts in y-axis baseline'''

    line, = ax.plot(x, y - yshift, **plot_kwargs)

    return line


def interpolate_xydata(x, y, n=3):
    '''Linearly interpolates an xy dataset adding n new points
       between each xy datapoint'''

    assert(x.ndim == y.ndim == 1 and x.size == y.size and n > 0)

    r = np.column_stack([x, y,])
    r_interp = np.zeros((n * (r.shape[0] - 1) + 1, r.shape[1],), dtype=r.dtype)
    r_interp[::n, :] = r

    for seg in range(r.shape[0] - 1):
        seg_sta = n * seg
        seg_end = n + seg_sta

        r_interp[seg_sta:seg_end] = np.column_stack([
            np.linspace(r[seg, 0], r[seg+1, 0], n, endpoint=False),
            np.linspace(r[seg, 1], r[seg+1, 1], n, endpoint=False),
        ])
    return r_interp.T


# simple argument parser
parser = ArgumentParser()
parser.add_argument('pickles', type=str, nargs='+')

markers = ['^', 'v', 'd', 'o',]

 # run this script in its current directory as, e.g.,:
 # ```python3 plot_fig3_averages.py './fig3a_10mer_data.pickle' './fig3a_20mer_data.pickle'```
if __name__ == '__main__':

    args = parser.parse_args()

    # each data dictionary should contain:
    # - monomer number `N` (str)
    # - swimming speed `Pe` (list[float, ...])
    # - end-to-end distance mean `<R>_mean` (list[float, ...])
    # - end-to-end distance standard error `<R>_std` (list[float, ...])
    # - end-to-end distance histogram (`r`, `P(r)`) (list[numpy.ndarray[numpy.float64]])
    data_dicts = [pickle.load(open(p, 'rb')) for p in args.pickles]

    # create figure with set panel dimensions
    fig_height = plot_height
    fig_width = plot_width
    fig, ax, = plt.subplots(nrows=1, ncols=1, layout='none',
    	                    figsize=(fig_width, fig_height,),)

    fig.subplots_adjust(top=top_adjust,
                        left=left_adjust,
                        right=right_adjust,
                        bottom=bottom_adjust,)

    # set axis limits and axis (tick)labels
    ax.set_xlim([-2.5, 62.5,])
    ax.set_xticks([0.0, 20.0, 40.0, 60.0,])
    ax.set_ylim([0.67, 1.03,])
    ax.set_yticks([0.7, 0.8, 0.9, 1.0,])

    ax.xaxis.set_tick_params(labelsize=ticklabel_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticklabel_fontsize)

    ax.set_xlabel(r'$\mathrm{Pe}$',
    	          fontsize=axislabel_fontsize, labelpad=1.0)
    ax.set_ylabel(r'$\langle R \rangle^\mathrm{neq}'
    	          r'/\langle R \rangle^\mathrm{eq}$',
                  fontsize=axislabel_fontsize, labelpad=1.0)

    # grab Pe values for each N-mer
    peclet_keys = {
        d['N']: d['peclet_keys']
        for d in data_dicts
    }
    peclet_values = {
        d['N']: np.asarray(d['peclet_values'], dtype=np.float64)
        for d in data_dicts
    }

    # grab mean end-to-end distances estimated from averaging over
    # long-time, steady-state trajectories
    mean_distances_traj = {
        d['N']: np.asarray(
            [d[k]['<R>'] for k in peclet_keys[d['N']]], dtype=np.float64
        ) for d in data_dicts
    }

    # estimate mean end-to-end distance by trapezoidal integration of 
    # the respective histogram at each Pe, for each N-mer
    def trapz_mean(x, f):
        return np.trapz(x * f, x=x)

    mean_distances_hist = {
        d['N']: np.asarray(
            [trapz_mean(d[k]['r'], d[k]['P(r)']) for k in peclet_keys[d['N']]], dtype=np.float64
        ) for d in data_dicts
    }

    # plot mean end-to-end distance values at each Pe, for each N-mer
    for marker, data_dict in zip(markers, data_dicts):
        N = data_dict['N']
        ax.plot(peclet_values[N], mean_distances_traj[N] / mean_distances_hist[N][0],
                label=f'$N = {N}$', color='k', marker=marker, **plot_kwargs)
  
    legend_bbox = [-0.03, -0.07, 0.0, 0.0]
    ax.legend(frameon=False, reverse=False, fontsize=legend_fontsize,
              loc='lower left', bbox_to_anchor=legend_bbox,
              markerscale=1.0, labelspacing=0.2, borderpad=0.2,
              handlelength=0.5, handleheight=1.0, handletextpad=0.3,)

    plt.show()
