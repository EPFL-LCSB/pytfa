import numpy as np
# from cobra.flux_analysis import flux_variability_analysis
from pytfa.analysis.variability import variability_analysis

def remove_blocked_reactions(model):
    epsilon = model.solver.configuration.tolerances.feasibility
    # fva = flux_variability_analysis(model)
    fva = variability_analysis(model, kind='reactions')
    # Blocked reactions have max and min at 0
    df = fva[ (fva.max(axis=1).abs()<1*epsilon)
            & (fva.min(axis=1).abs()<1*epsilon)]
    rid_to_rm = df.index

    model.remove_reactions(rid_to_rm)

    return df


def round(value, epsilon):
    n = int(-1*np.log10(epsilon))
    return np.round(value,n)

def trim_epsilon_mets(met_dict, epsilon, model=None):

    round_dict = {x:round(v,epsilon) for x,v in met_dict.items()}
    met_dict.update(round_dict)

    if model is None:
        rm_list = [x for x,v in met_dict.items() if abs(v) <= epsilon]
    else:
        rm_list = [x for x,v in met_dict.items() if abs(v) <= epsilon and x not in model.metabolites]

    [met_dict.pop(x) for x in rm_list]

    return met_dict

def set_medium(model, medium_dict, inplace):
    if inplace:
        new = model
    else:
        new = model.copy()

    if medium_dict is None or not medium_dict:
        return new

    for rxn_id, lb in medium_dict.items():
        new.reactions.get_by_id(rxn_id).lower_bound = lb

    return new