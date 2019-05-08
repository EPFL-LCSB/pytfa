from cobra.flux_analysis import flux_variability_analysis

def remove_blocked_reactions(model):
    epsilon = model.solver.configuration.tolerances.feasibility
    fva = flux_variability_analysis(model)
    # Blocked reactions have max and min at 0
    rid_to_rm = fva[  (fva.max(axis=1).abs()<1*epsilon)
                    & (fva.min(axis=1).abs()<1*epsilon)].index

    model.remove_reactions(rid_to_rm)


def trim_epsilon_mets(reaction, epsilon):
    rm_dict = {x:-v for x,v in reaction.metabolites.items() if abs(v)<=epsilon}
    reaction.add_metabolites(rm_dict)