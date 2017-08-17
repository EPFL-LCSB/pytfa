Quick start
===========

Three tutorial files detail thoroughly normal usages of the pytfa package. They
can be found at::

    pytfa
    └── tutorials
        ├── figure_paper.py
        ├── glycolysis_example.py
        └── tutorial_sampling.py

`figure_paper.py` details how to get the figure from our paper [CITE], a simple
use case for TFA

`glycolysis_example.py` shows a more realistic case with a reduced model of
yeast, highlights differences betweeen FBA and TFA. It also cycles through several
solvers (if more are installed), to show how simple it is to change your solver
(thanks to [optlang](https://github.com/biosustain/optlang) ). Finally it generates
plots to visualize results.

`tutorial_sampling.py` shows how to sample a variable, for example thermodynamic
displacement, and generate plots to visualize the results.

The next sections give more details on how the thermodynamic model is
structured, and how data is managed.

Cheers,

The py.TFA team