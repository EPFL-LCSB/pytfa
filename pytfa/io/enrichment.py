"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Tools to import or export enrichment to and from pytfa models


"""
import json

import numpy
import pandas as pd

class MyEncoder(json.JSONEncoder):
    """
    We define an encoder that takes care of the serialization of numpy types,
    which are not handled by json by default
    """
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

def write_lexicon(tmodel, filepath):
    """
    Writes a csv file in the format :
                    seed_id
    13BDglcn_c    cpd11791
    13dpg_c       cpd00203
    2pg_c         cpd00482
    3pg_c         cpd00169
    4abut_c       cpd00281

    Useful for exporting an annotation

    :type tmodel: pytfa.core.ThermoModel
    :param tmodel:
    :param filepath:
    :return:
    """

    lexicon = pd.DataFrame.from_dict({x.id:x.annotation    \
                                   for x in tmodel.metabolites},
                                  orient = 'index')
    lexicon.to_csv(filepath)
    return lexicon


def annotate_from_lexicon(model,lexicon):
    """
    Converts a lexicon into annotation for the metabolites

    :type model: cobra.Model
    :param model:
    :param lexicon:
    :return:
    """
    annotations = lexicon.to_dict(orient = 'index')

    for this_metabolite in model.metabolites:
        try:
            this_metabolite.annotation = annotations[this_metabolite.id]
        except KeyError:
            model.logger.warning(this_metabolite.id +   \
                                 '  not found in annotations')


def read_lexicon(filepath):
    return pd.read_csv(filepath, index_col = 0)


def write_compartment_data(tmodel, filepath):
    """

    :type tmodel: pytfa.core.ThermoModel
    :param tmodel:
    :return:
    """

    # compartment_data = pd.DataFrame.from_dict(tmodel.compartments,
    #                                           orient = 'index')
    # compartment_data.to_csv(filepath)
    # return compartment_data
    filepath = filepath + '.json' if not filepath.endswith('.json') \
            else filepath
    with open(filepath, 'w') as outfile:
        json.dump(tmodel.compartments, outfile, cls=MyEncoder)


def read_compartment_data(filepath):
    #return pd.read_csv(filepath, index_col = 0)
    filepath = filepath + '.json' if not filepath.endswith('.json') \
            else filepath
    with open(filepath) as json_data:
        compartment_data = json.load(json_data)
        return compartment_data


def apply_compartment_data(tmodel,compartment_data):
    # tmodel.compartments = compartment_data.to_dict(orient='index')
    tmodel.compartments = compartment_data


