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


def trim_epsilon_mets(reaction, epsilon):
    rm_dict = {x:-v for x,v in reaction.metabolites.items() if abs(v)<=epsilon}
    reaction.add_metabolites(rm_dict)

    n = int(-1*np.log10(epsilon))
    round_dict = {x:-v+np.round(v,n) for x,v in reaction.metabolites.items()}
    reaction.add_metabolites(round_dict)

def set_medium(model, medium_dict, inplace):
    if inplace:
        new = model
    else:
        new = model.copy()

    for rxn_id, lb in medium_dict.items():
     new.reactions.get_by_id(rxn_id).lower_bound = lb

    return new