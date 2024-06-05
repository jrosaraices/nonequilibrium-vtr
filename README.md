# Data supplement: _Variational time reversal for free energy estimation in nonequilibrium steady states_

> __NOTE__: This repository is still under construction, expected to be complete by _Wednesday, Jun 5 at 9PM PST_.

This repository contains the numerical results, with generating code, presented in the journal article titled ___Variational time reversal for free energy estimation in nonequilibrium steady states___, authored by [Jorge L. Rosa-Raíces](mailto:jrosaraices@berkeley.edu) and [David T. Limmer](mailto:dlimmer@berkeley.edu) at UC Berkeley and submitted for publication in Physical Review E.  The article is currently available [on arXiv](https://arxiv.org/abs/2406.01582).  This repository is organized as follows:

- Folder `numerical_results` encloses post-processed simulation data plotted throughout the article, along with the _Python 3_ scripts used to generate figures for the article.
- Folder `simulation_code` encloses a pure _Fortran_ application written to comply with the [2018 language standard](https://wg5-fortran.org/f2018.html), and with no external library dependencies for portability.  The Makefile provided compiles the simulation code on macOS Sonoma with Command Line Tools ≤ 14.3.1 using either [Intel's legacy Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) `ifort`, or [NAG's Fortran Compiler](https://nag.com/fortran-compiler/) `nagfor`.  The [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) has not been tested to compile the provided code.

The raw trajectory data produced by the simulation code and post-processed to generate the numerical results is available upon reasonable request from the authors. <!--Input files are provided to generate the same data, together with processing scripts, are included together with the processed data files.-->
