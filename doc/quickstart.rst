Quick start
===========

Three tutorial files detail thoroughly normal usages of the pytfa package. They
can be found at::

    pytfa
    └── tutorials
        ├── figure_paper.py
        ├── tutorial_basics.py
        └── tutorial_sampling.py

`figure_paper.py` details how to get the figure from our paper [1]_, a simple
use case for TFA on a reduced *Escherichia coli*. We show that adding thermodyamics
constraints and simple concentration data allow to substantially reduce the flux space.

`tutorial_basics.py` shows a more realistic case with two models (reduced or full genome-scale) of
*Escherichia coli*. It also cycles through several
solvers (if more are installed), to show how simple it is to change your solver
(thanks to `optlang`_).

`tutorial_sampling.py` shows how to sample a variable, for example thermodynamic
displacement, and generate plots to visualize the results.

**If you plan to run the tutorials with full genome-scale models, we recommend you to get a commercial
solver, as it has been seen that GLPK's lack of parallelism significantly increases solving time**

The next sections give more details on how the thermodynamic model is
structured, and how data is managed.

Cheers,

The py.TFA team

.. _optlang: https://github.com/biosustain/optlang

.. [1]
   Salvy, P., Fengos, G., Ataman, M., Pathier, T., Soh, K. C., & Hatzimanikatis, V. (2018). 
   pyTFA and matTFA: a Python package and a Matlab toolbox for Thermodynamics-based Flux Analysis. Bioinformatics, 35(1), 167-169.
