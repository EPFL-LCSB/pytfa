"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Input/Output tools to vizualize results


"""
import pandas as pd

from ..utils.str import varnames2ids


def export_variable_for_escher(tmodel,variable_type,data,filename):
    """
    Exports all the variables of a given type into a csv file, indexed by
    variable.id. This format is readable by escher if the variable_type is a
    subclass of :pytfa:`pytfa.optim.variables.ReactionVariable` or
    :pytfa:`pytfa.optim.variables.MetaboliteVariable`

    :param pytfa.core.ThermoModel tmodel:
    :param ReactionVariable|MetaboliteVariable variable_type:
    :param pandas.Series data: indexed by variable name
    :param string filename:
    :return:
    """
    var_dict = {x.name: x.id   \
                     for x in tmodel.get_variables_of_type(variable_type)}

    if isinstance(data,pd.DataFrame):
        # it is a analysis
        var_values = data[:, var_dict.keys()]
    elif isinstance(data,pd.Series):
        var_values = data[var_dict.keys()]

    var_values.index = varnames2ids(tmodel, var_values.index)
    var_values.to_csv(filename)


def get_reaction_data(tmodel, data):
    """
    Exports values indexed by reaction ids. Reconciles Forward and Backwards
    variables.
    """
    dict_values = {}
    for reaction in tmodel.reactions:
        fwd_val = data[reaction.forward_variable.name]
        bwd_val = data[reaction.reverse_variable.name]
        dict_values[reaction.id] = fwd_val - bwd_val

    var_values = pd.Series(dict_values)
    return  var_values


def export_reactions_for_escher(tmodel,data,filename):
    """
    Exports values indexed by reaction ids. Reconciles Forward and Backwards
    variables. Writes it in a csv file. This format is readable by escher

    :param pytfa.core.ThermoModel tmodel:
    :param ReactionVariable|MetaboliteVariable variable_type:
    :param pandas.Series data: indexed by variable name
    :param string filename:
    :return:
    """
    var_values = get_reaction_data(tmodel, data)
    var_values.to_csv(filename,header=True)