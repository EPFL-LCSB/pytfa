# pyTFA

Thermodynamic constraints for flux-based reactions, in Python.

## Requirements

This module requires `cobrapy` to work properly. You might also want to install
a dedicated solver. GLPK, CPLEX and Gurobi are supported.

## Container-based install

You might want to use this program inside of a container. the `docker/`
subfolder has all the necessary information and source files to set it up.

# Setup

**NOTE :** This step is not required if you're using the container, which bundles all this

You can install this module with `pip`:

```
pip install -e /path/to/pytfa
```

## Generating the docs

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
