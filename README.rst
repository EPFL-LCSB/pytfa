pyTFA
=====
|PyPI| |Documentation Status| |Build Status| |Codecov| |Codacy branch grade| |license| 

Thermodynamics-based Flux Analysis, in Python.
Paper : Pierre Salvy, Georgios Fengos, Meric Ataman, Thomas Pathier, Keng C Soh, Vassily Hatzimanikatis. "pyTFA and matTFA: a Python package and a Matlab toolbox for Thermodynamics-based Flux Analysis"
Bioinformatics (2018), bty499, `DOI:
https://doi.org/10.1093/bioinformatics/bty499 <https://doi.org/10.1093/bioinformatics/bty499>`_

Implements: Christopher S. Henry, Linda J. Broadbelt, and Vassily
Hatzimanikatis. "Thermodynamics-based metabolic flux analysis."
Biophysical journal 92.5 (2007): 1792-1805. `DOI:
https://doi.org/10.1529/biophysj.106.093138 <https://doi.org/10.1529/biophysj.106.093138>`__

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
    cd /path/to/pytfa
    git lfs install
    git lfs pull

**This module was developed in Python 3.5, and it is recommended to run Python 3.5 
to run commercial solvers such as Gurobi and CPLEX.**
Other Python versions (2.7, 3.4) should also work (see the `CI builds <https://travis-ci.org/EPFL-LCSB/pytfa>`_)


This module requires
`COBRApy <https://github.com/opencobra/cobrapy/>`_, as well as
`optlang <https://github.com/biosustain/optlang>`_ to work
properly. The installer should take care of that for you. You might also
want to install a dedicated solver. GLPK, CPLEX and Gurobi are
supported.

Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/pytfa/tree/master/docker

Setup
=====

*This step is not required if you're using the container, which bundles all this.*

You can install this module with ``pip``:

*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    pip3 install pytfa

or from source

.. code:: bash

    git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
    pip3 install -e /path/to/pytfa

Quick start
===========

Three tutorial files detail thoroughly normal usages of the pytfa
package. They can be found at:

::

    pytfa
    └── tutorials
        ├── figure_paper.py
        ├── tutorial_basics.py
        └── tutorial_sampling.py

More information can be found
`here <http://pytfa.readthedocs.io/en/latest/quickstart.html>`__.

Documentation
=============

Documentation is hosted at `Read the
Docs <http://pytfa.readthedocs.io/en/latest/index.html>`__

Alternatively you can also generate the docs locally.

Make sure `sphinx <https://www.sphinx-doc.org/en/stable/>`__ is
installed, and install as well `the
theme <https://github.com/rtfd/sphinx_rtd_theme>`__ (this is already
bundled with the container):

::

    pip install sphinx sphinx-rtd-theme

You can then generate the documentation with this command:

::

    cd work/pytfa/doc && make html

The resulting HTML files will be located in ``work/pytfa/doc/_build``.

Testing the code
----------------

We recommend using the `Docker
container <https://github.com/EPFL-LCSB/pytfa/tree/master/docker>`__ for
testing the code, as it comes with everything bundled.

Install `pytest <https://docs.pytest.org/en/latest/>`__ if you don't
already have it (``pip install pytest``, already included in the
container), then start the tests with the ``pytest`` command.

Usage
=====

First, create your COBRApy model. Make sure to define the additional
values required by pyTFA, as said in the "Models" page of the
documentation.

If you already have a Matlab model with thermodynamic data, you might
want to use ``pytfa.io.import_matlab_model``. Otherwise, have a look at
the `COBRApy
documentation <https://cobrapy.readthedocs.io/en/latest/io.html#MATLAB>`__,
then add the required properties.

If you're using a specific solver, don't forget to tell COBRApy about it
by setting the ``solver`` property of your model to the name of your
solver. See the `COBRApy
documentation <https://cobrapy.readthedocs.io/en/latest/solvers.html>`__
for more information about this.

Thermodynamic database
----------------------

You also need a thermodynamic database. Use ``thermoDBconverter.py`` if
you have a thermodynamic database from Matlab you wish to import to
Python.

Thermodynamic databases are stored in ``.thermodb`` files and can be
easily loaded with ``pytfa.io.load_thermoDB``.

Example script
--------------

Here is an example script :

.. code:: python

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

.. |PyPI| image:: https://img.shields.io/pypi/v/pytfa.svg
   :target: https://pypi.org/project/pytfa/
.. |Documentation Status| image:: https://readthedocs.org/projects/pytfa/badge/?version=latest
   :target: http://pytfa.readthedocs.io/en/latest/?badge=latest
.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt
.. |Build Status| image:: https://travis-ci.org/EPFL-LCSB/pytfa.svg?branch=master
   :target: https://travis-ci.org/EPFL-LCSB/pytfa
.. |Codecov| image:: https://img.shields.io/codecov/c/github/EPFL-LCSB/pytfa.svg
   :target: https://codecov.io/gh/EPFL-LCSB/pytfa
.. |Codacy branch grade| image:: https://img.shields.io/codacy/grade/d8fd67ee134d46a69115c9b39c19be26/master.svg
   :target: https://www.codacy.com/app/realLCSB/pytfa
.. |Code climate| image:: https://img.shields.io/codeclimate/github/EPFL-LCSB/pytfa.svg
   :target: https://codeclimate.com/github/EPFL-LCSB/pytfa
   
   
License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt>`_ file for more details.
