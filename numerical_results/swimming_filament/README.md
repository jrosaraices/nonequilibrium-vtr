Enclosed are the numerical results for the swimming filament, together with plotting scripts to generate Figs. 3-4 of [the article on arXiv](http://arxiv.org/abs/2406.01582).

The folder `fig3` contains the following python scripts:
- `plot_fig3_averages.py` generates Fig. 3B from the data stored in `fig3b_averages_10mer.pickle` and `fig3b_averages_20mer.pickle`
- `plot_fig3_free_energies.py` generates Fig. 3C from the data stored in `fig3_free_energies_10mer.pickle` and `fig3_free_energies_20mer.pickle`.  The pickled data is stored in `numpy.NDArray` objects in a simple dictionary structure.

Similarly, the folder `fig4` contains `plot_fig4_averages.py`, which plots optimized policy parameters pickled in `fig4_policy_params_10mer.pickle` and `fig4c_policy_params_20mer.pickle` for the swimming 10-mer and 20-mer, respectively.

Panels generated with these scripts were separately stored and manually arranged in a vector graphics editor to form Figs. 3-4.
