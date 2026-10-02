Thermodynamic Databases
=======================


Converting a Matlab thermodynamic database
------------------------------------------

If you have a Matlab thermodynamic database, you can easily convert it to a
Python database thanks to the script ``thermoDBconverter.py``::

  python thermoDBconverter.py database.mat converted_database.thermodb

Loading a thermodynamic database
--------------------------------

Thermodynamic databases are compressed through :any:`zlib` and binary-encoded
with :any:`pickle`. In order to load them, you need first to uncompress them
with :any:`zlib.decompress` then load the result into memory with
:any:`pickle.loads`::

  import pickle
  import zlib
  with open('thermoDatabases/DB_AlbertyUpdate.thermodb', 'rb') as file:
      ReactionDB = pickle.loads(zlib.decompress(file.read()))

.. warning::
  Since the file is compressed, you **MUST** load it as a binary file by
  calling :any:`open` with the ``b`` flag, otherwise Python will try to
  decode it as unicode and raise an exception !

Building a database from the ModelSEEDDatabase
----------------------------------------------

:mod:`pytfa.thermo.modelseed` builds a thermodynamic database from a local
clone of the `ModelSEEDDatabase <https://github.com/ModelSEED/ModelSEEDDatabase>`_,
located with the ``db_path`` argument, the ``MODELSEED_DB_PATH`` environment
variable, or the configuration of modelseedpy (``pip install .[modelseed]``)::

  from pytfa.io import apply_compartment_data
  from pytfa.thermo.modelseed import annotate_seed_ids, \
      build_thermo_from_modelseed, build_compartment_data

  annotate_seed_ids(cobra_model)
  thermo_data = build_thermo_from_modelseed(cobra_model)
  tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
  apply_compartment_data(tmodel, build_compartment_data(tmodel))

* ``annotate_seed_ids`` sets the ``seed_id`` annotation from ModelSEED
  metabolite ids (``cpd00002_c0``) or from ``seed.compound`` annotations.
* The energies are the ``Group contribution`` estimates of the compounds, in
  kcal/mol. Each is paired with the protonation state of the metabolite in the
  model, shifted along the compound's pKas so that its transformed energy is
  unchanged. :meth:`~pytfa.ThermoModel.prepare` then leaves the reactions as the
  model balances them.
* The pKas are Marvin's acidic pKas.
* There are no cues: reaction errors combine the compound errors in quadrature.
* ``build_compartment_data`` defaults to modelseedpy's values: pH 7 and -160 mV
  in the cytosol, pH 6.5 and 0 mV in the extracellular space, and
  concentrations between 1e-6 and 0.02 M.

See ``tutorials/tutorial_modelseed.py``.

Structure of a thermodynamic database
-------------------------------------

A thermodynamic database is a :any:`dict` with the following fields:

  * ``name`` : :any:`string` The name of the database
  * ``units`` : :any:`string` The unit of the energies in the database. Can be
    ``kcal/mol`` or ``kJ/mol``.
  * ``metabolites`` : :py:class:`dict` A dictionnary containing the metabolites'
    thermodynamic data. See :ref:`metabolites` for more information.
  * ``cues`` :  :py:class:`dict` A dictionnary containing the cues'
    thermodynamic data. See :ref:`cues` for more information.


.. _metabolites:

Metabolites
^^^^^^^^^^^
This is a dictionnary storing various thermodynamic data about metabolites. It
is stored as a :py:class:`dict` where each key is a ``SeedID``. The values are
others :py:class:`dict` with the following keys.

+-------------+---------------------------------------------------------------+
| id          | :any:`string` SeedID of the metabolite                        |
+-------------+---------------------------------------------------------------+
| charge_std  | :any:`float` Charge of the metabolite (mV) in standard        |
|             | conditions                                                    |
+-------------+---------------------------------------------------------------+
| deltaGf_std | :any:`float` Transformed Gibbs energy of formation of the     |
|             | metabolite, in standard conditions.                           |
+-------------+---------------------------------------------------------------+
| deltaGf_err | :any:`float` Error on the transformed Gibbs energy of         |
|             | formation of the metabolite, in standard conditions           |
+-------------+---------------------------------------------------------------+
| mass_std    | :any:`float` Mass of the metabolite (g.mol-1)                 |
+-------------+---------------------------------------------------------------+
| nH_std      | :any:`int` Number of protons of the metabolite, in standard   |
|             | conditions                                                    |
+-------------+---------------------------------------------------------------+
| error       | :any:`string` Error on the metabolite's thermodynamic data.   |
|             | Thermodynamic values will be computed only if this equals to  |
|             | 'Nil'.                                                        |
+-------------+---------------------------------------------------------------+
| formula     | :any:`string` Formula of the metabolite.                      |
+-------------+---------------------------------------------------------------+
| nH_std      | :any:`int` Number of protons in the metabolite's formula      |
+-------------+---------------------------------------------------------------+
| name        | :any:`string` Name of the metabolite                          |
+-------------+---------------------------------------------------------------+
| other_names | :any:`list` (:any:`string`) Other names of the metabolite     |
+-------------+---------------------------------------------------------------+
| pKa         | :any:`list` (:any:`float`) pKas of the metabolite             |
+-------------+---------------------------------------------------------------+
| struct_cues | :py:class:`dict` (:any:`int`) cues of the metabolite          |
|             |                                                               |
|             | The keys of the array are the names of the cues, and the      |
|             | values the number of cues of this type that are part of the   |
|             | structure.                                                    |
+-------------+---------------------------------------------------------------+

Here is an example::

  ReactionDB['metabolites']['cpd00001'] = {
    'charge_std': 0,
    'deltaGf_err': 0.5,
    'deltaGf_std': -56.686999999999998,
    'error': 'Nil',
    'formula': 'H2O',
    'id': 'cpd00001',
    'mass_std': 18.0,
    'nH_std': 2,
    'name': 'H2O',
    'other_names': ['H2O', 'Water', 'HO-', 'OH-', 'h2o'],
    'pKa': [15.7],
    'struct_cues': {'H2O': 1}
  }


.. _cues:

Cues
^^^^

This is a dictionnary storing various thermodynamic data about cues. It
is stored as a :py:class:`dict` where each key is the cue ID, as referrenced
in the ``struct_cues`` attribute of :ref:`metabolites`. The values are
others :py:class:`dict` with the following keys.

+---------+-------------------------------------------------------------------+
| id      | :any:`string` ID of the cue                                       |
+---------+-------------------------------------------------------------------+
| charge  | :any:`float` The charge (mV) of the cue in standard conditions    |
+---------+-------------------------------------------------------------------+
| datfile | :any:`string` The dat file from which the data was imported.      |
|         | `Optional`                                                        |
+---------+-------------------------------------------------------------------+
| energy  | :any:`float` Transformed Gibbs energy of formation of the cue, in |
|         | standard conditions.                                              |
+---------+-------------------------------------------------------------------+
| error   | :any:`float` The error on the transformed Gibbs energy of         |
|         | formation of the cue, in standard conditions.                     |
+---------+-------------------------------------------------------------------+
| formula | :any:`string` Formula of the cue                                  |
+---------+-------------------------------------------------------------------+
| names   | :any:`list` (:any:`string`) Other names of the cue                |
+---------+-------------------------------------------------------------------+
| small   | :any:`bool` Whether this is a small cue or not                    |
+---------+-------------------------------------------------------------------+

Here is an example::

  ReactionDB['cues']['H2O'] = {
    'charge': 0,
    'datfile': 'H2O.gds',
    'energy': -56.686999999999998,
    'error': 0.5,
    'formula': 'H2O',
    'id': 'H2O',
    'names': ['H2O', 'OH-', 'HO-'],
    'small': True
  }
