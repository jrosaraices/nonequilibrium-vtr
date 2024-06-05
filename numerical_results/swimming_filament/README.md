Enclosed are the numerical results for the swimming filament, presented in Figs. 3-4 of [the article on arXiv](http://arxiv.org/abs/2406.01582).

The folder `fig3` contains the python scripts `plot_fig3_averages.py`, which generates Fig. 3B from the data stored in `fig3b_data.pickle`, and `plot_fig3_free_energies.py`, which can used to generate both Fig. 3C from the data stored in `fig3c_data.pickle`, and Fig. 3D from the data stored in `fig3d_data.pickle`.  All pickled data is stored in `numpy.NDArray` objects in a simple dictionary structure.

Similarly, the folder `fig4` contains `plot_fig4_averages.py`, which plots optimized policy parameters stored in pickles `fig4b_data.pickle` and `fig4c_data.pickle` for the swimming 10-mer and 20-mer, respectively.

Panels generated with these scripts were separately saved and manually arranged in a vector graphics editor to generate Figs. 3-4.
