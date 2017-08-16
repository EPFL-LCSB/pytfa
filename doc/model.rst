pyTFA models
============

pyTFA models are based on :class:`COBRApy models <cobra.core.model.Model>`, with
additional values.

Compartment data
----------------
This is the ``compartments`` attribute of the model. It is a :any:`dict` where
each key is the symbol of a compartment, and the value is another :any:`dict`
with the following keys :

+-------------+---------------------------------------------------------------+
| c_min       | :any:`float` The minimum concentration for each               |
|             | metabolite in the compartment, in mol.L-1                     |
+-------------+---------------------------------------------------------------+
| c_max       | :any:`float` The maximum concentration for each               |
|             | metabolite in the compartment, in mol.L-1                     |
+-------------+---------------------------------------------------------------+
| ionicStr    | :any:`float` The ionic strength of the compartment (mV)       |
+-------------+---------------------------------------------------------------+
| membranePot | :py:class:`dict` A dictionnary representing the membrane      |
|             | potential between this compartment (which is the source       |
|             | compartment) and the others.                                  |
|             |                                                               |
|             | Each key is the symbol of another compartment (which is the   |
|             | destination compartment), and the value is the potential      |
|             | (in mV) from the source to the destination.                   |
+-------------+---------------------------------------------------------------+
| name        | :any:`string` The name of the compartment                     |
+-------------+---------------------------------------------------------------+
| pH          | :any:`float` The pH in the compartment                        |
+-------------+---------------------------------------------------------------+
| symbol      | :any:`string` The symbol of the compartment (which is         |
|             | the key of this dictionnary)                                  |
+-------------+---------------------------------------------------------------+

Here is an example::

  cobra_model.compartments['c'] = {
    'c_max': 0.01,
    'c_min': 9.9999999999999995e-08,
    'ionicStr': 0.25,
    'membranePot': {
      'c': 0,
      'e': 60,
      'g': 0,
      'm': -180,
      'n': 0,
      'p': 0,
      'r': 0,
      't': 0,
      'v': 0,
      'x': 0
    },
    'name': 'Cytosol',
    'pH': 7.0,
    'symbol': 'c'
  }


Metabolites
-----------
Each metabolite must be annotated with its ``SeedID``, which will be used to get
the thermodynamic values from the :doc:`thermoDB`. In order to do this, use the
``annotation`` attribute of each
:class:`metabolite <cobra.core.metabolite.Metabolite>`. Here is an example::

  cobra_model.metabolites[0].annotation = {
    'SeedID': 'cpd00018'
  }

pyTFA will also define a ``thermo`` a thermo attribute for each metabolite,
which is a :class:`pytfa.thermo.MetaboliteThermo`.

Reactions
---------
pyTFA will define a ``thermo`` attribute for each reaction. It is a
:py:class:`dict` with the following attributes:

+------------+-----------------------------------------------------------------+
| computed   | :any:`bool` Whether the thermodynamic values were computed or   |
|            | not.                                                            |
+------------+-----------------------------------------------------------------+
| deltaGR    | :any:`float` Sum of the non-concentration terms for the         |
|            | reaction. Used as the right-hand side value of a constraint.    |
|            |                                                                 |
|            | If the thermodynamic values were not computed, this is          |
|            | ``10**7``.                                                      |
+------------+-----------------------------------------------------------------+
| deltaGRerr | :any:`float` Error on ``deltaGR``                               |
|            |                                                                 |
|            | If the thermodynamic values were not computed, this is          |
|            | ``10**7``.                                                      |
+------------+-----------------------------------------------------------------+
| deltaGrxn  | :any:`float` Sum of the deltaGF of all the metabolites in the   |
|            | reaction. **Not defined if computed is False !**                |
+------------+-----------------------------------------------------------------+
| isTrans    | :any:`bool` Whether the reaction is a transport reaction or not |
+------------+-----------------------------------------------------------------+

Here are some examples::

  cobra_model.reactions[0].thermo = {
    'computed': False,
    'deltaGR': 10000000,
    'deltaGRerr': 10000000,
    'isTrans': False
  }

  cobra_model.reactions[99].thermo = {
    'computed': True,
    'deltaGR': 1.161097833014658,
    'deltaGRerr': 2,
    'deltaGrxn': 0,
    'isTrans': True,
  }
