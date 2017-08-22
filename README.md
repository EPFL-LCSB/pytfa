# pyTFA

[![PyPI](https://img.shields.io/pypi/v/pytfa.svg)](https://pypi.org/project/pytfa/)
[![Documentation Status](https://readthedocs.org/projects/pytfa/badge/?version=latest)](http://pytfa.readthedocs.io/en/latest/?badge=latest)
[![license](http://img.shields.io/badge/license-APACHE2-blue.svg)](https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt)


[![Build Status](https://travis-ci.org/EPFL-LCSB/pytfa.svg?branch=master)](https://travis-ci.org/EPFL-LCSB/pytfa)
[![Codecov](https://img.shields.io/codecov/c/github/codecov/pytfa.svg)](https://codecov.io/gh/EPFL-LCSB/pytfa)
[![Codacy branch grade](https://img.shields.io/codacy/grade/d8fd67ee134d46a69115c9b39c19be26/master.svg)]()

Thermodynamics-based Flux Analysis, in Python.

Implements:
Henry, Christopher S., Linda J. Broadbelt, and Vassily Hatzimanikatis.
"Thermodynamics-based metabolic flux analysis."
Biophysical journal 92.5 (2007): 1792-1805.
[DOI: https://doi.org/10.1529/biophysj.106.093138](https://doi.org/10.1529/biophysj.106.093138)

## Requirements

**NOTE :** This module requires Python 3.5. Other Python versions might work, but are not officially supported.

This module requires [`cobrapy`](https://github.com/opencobra/cobrapy/), as well as [`optlang`](https://github.com/biosustain/optlang) to work properly. The installer should take care of that for you. You might also want to install
a dedicated solver. GLPK, CPLEX and Gurobi are supported.

## Container-based install

You might want to use this program inside of a container. The [`docker/`](https://github.com/EPFL-LCSB/pytfa/tree/master/docker)
subfolder has all the necessary information and source files to set it up.

# Setup

**NOTE :** This step is not required if you're using the container, which bundles all this.

You can install this module with `pip`:

**NOTE :** Because this module requires Python 3, you might have to use `pip3` instead of `pip`

```bash
pip3 install pytfa
```
or from source
```bash
git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
pip3 install -e /path/to/pytfa
```

# Quick start
Three tutorial files detail thoroughly normal usages of the pytfa package. They can be found at:
```
pytfa
└── tutorials
    ├── figure_paper.py
    ├── relaxation_example.py
    └── tutorial_sampling.py
```

More information can be found [here](http://pytfa.readthedocs.io/en/latest/quickstart.html).

# Documentation

Documentation is hosted at [Read the Docs](http://pytfa.readthedocs.io/en/latest/index.html)

Alternatively you can also generate the docs locally.

Make sure [sphinx](https://www.sphinx-doc.org/en/stable/) is installed, and
install as well [the theme](https://github.com/rtfd/sphinx_rtd_theme) (this is
already bundled with the container):

```
pip install sphinx sphinx-rtd-theme
```

Then begin the generation :

```
cd work/pytfa/doc && make html
```

The resulting HTML files are located in `work/pytfa/doc/_build`.

## Testing the code

We recommand using the [Docker container](https://github.com/EPFL-LCSB/pytfa/tree/master/docker) for testing the code, as it comes with everything bundled.

Install [pytest](https://docs.pytest.org/en/latest/) if you don't already have
it (`pip install pytest`, already included in the container), then start the
tests with the `pytest` command.

# Usage

First, create your COBRApy model. Make sure to define the additional values
required by pyTFA, as said in the "Models" page of the documentation.

If you already have a Matlab model with thermodynamic data, you might want to
use `pytfa.io.import_matlab_model`. Otherwise, have a look at the [COBRApy
documentation](https://cobrapy.readthedocs.io/en/latest/io.html#MATLAB), then
add the required properties.

If you're using a specific solver, don't forget to tell COBRApy about it by
setting the `solver` property of your model to the name of your solver. See the
[COBRApy documentation](https://cobrapy.readthedocs.io/en/latest/solvers.html)
for more information about this.

## Thermodynamic database

You also need a thermodynamic database. Use `thermoDBconverter.py` if you have
a thermodynamic database from Matlab you wish to import to Python.

Thermodynamic databases are stored in `.thermodb` files and can be easily loaded
with `pytfa.io.load_thermoDB`.

## Find a solution

First, prepare your model for TFBA analysis with `pytfa.prepModelForTFBA`. Then,
call `pytfa.convToTFA` to add the thermodynamic-based constraints.

Finally, use the `optimize` method of your model to find an optimal solution.

## Example script

Here is an example script :

```python
import pytfa
from pytfa.io import import_matlab_model, load_thermoDB


cobra_model = import_matlab_model('../models/small_yeast.mat')

thermo_data = load_thermoDB('../data/thermo_data.thermodb')

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.solver = 'optlang-cplex'

## TFA conversion
mytfa.prepare()
mytfa.convert()

## Info on the model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()
```
