Models
======

This folder contains the models used in the tutorials, as well as some
pre-curated thermodynamic information. Models that have such pre-curated
thermodynamic information will have a subfolder to their name, in which the
information is stored.

Pre-curated models
------------------

These models have been pre-curated for thermodynamics analysis, which mean they already have the necessary annotations to be mapped to our thermodynamic data. Here is a short example for ``small_ecoli.mat``:


.. code-block:: python

    from cobra.io import load_matlab_model

    from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data


    # Load the model
    cobra_model = load_matlab_model('../models/small_ecoli.mat')
    
    # Load reaction DB
    thermo_data = load_thermoDB('../data/thermo_data.thermodb')
    lexicon = read_lexicon('../models/small_ecoli/lexicon.csv')
    compartment_data = read_compartment_data('../models/small_ecoli/compartment_data.json')

    # Initialize the model
    tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
    tmodel.name = 'tutorial'
    
    # Annotate the model
    annotate_from_lexicon(tmodel, lexicon)
    apply_compartment_data(tmodel, compartment_data)

    ## TFA conversion
    tmodel.prepare()
    tmodel.convert()

Publications
------------

Some models have been published:

+-----------------+-------------------------------------------------------------------------+
| small_ecoli.mat | Ataman, M., Gardiol, D. F. H., Fengos, G., & Hatzimanikatis, V. (2017). |
|                 | redGEM: Systematic reduction and analysis of genome-scale metabolic     |
|                 | reconstructions for development of consistent core metabolic models.    |
|                 | PLoS computational biology, 13(7), e1005444.                            |
+-----------------+-------------------------------------------------------------------------+



