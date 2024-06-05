Enclosed are the numerical results for the sheared particle in a confining potential, together with plotting scripts to generate Figs. 1-2 of [the article on arXiv](http://arxiv.org/abs/2406.01582).

The folder `fig1` (resp. `fig2`) contains the python script `plot_fig1_panel.py` (resp. `plot_fig2_panel.py`) used to generate each panel in Fig. 1 (resp. Fig. 2).  Raw numerical data plotted with these scripts are included separately for each panel as python `pickle`s (e.g., `fig1a_data.pickle`).  Each data pickle contains a `dict` comprised of:

- the panel title (e.g., "Pe = 16");
- a `dict` with the plotted data in a `numpy.NDArray`, whether it be an *xy* dataset for the curves in Fig. 1 or basis function amplitudes for the heatmaps in Fig. 2; and
- a `dict` of plotting arguments passed to `matplotlib.axes.Axes.plot` to plot each curve in the panel.

Panels generated with these scripts were separately stored and manually arranged in a vector graphics editor to form Figs. 1-2.
