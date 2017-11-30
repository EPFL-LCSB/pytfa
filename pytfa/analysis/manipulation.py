def apply_reaction_variability(tmodel, va, inplace = True):
    """
    Applies the VA results as bounds for the reactions of a cobra_model
    :param inplace:
    :param tmodel:
    :param va:
    :return:
    """
    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()

    for this_reaction in _tmodel.reactions:
        if this_reaction.id not in va.index:
            continue

        this_reaction.lower_bound = va.loc[this_reaction.id,'minimum']
        this_reaction.upper_bound = va.loc[this_reaction.id,'maximum']

    return _tmodel


def apply_generic_variability(tmodel,va, inplace = True):
    """
    Reactions a dealt with cobra, but the other variables added use pytfa's
    interface: the class GenericVariable. We use a different method to apply
    variability directly in the solver

    :param tmodel:
    :param va:
    :param inplace:
    :return:
    """
    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()


    for varname in va.index:
        the_min,the_max = va.loc[varname,['minimum','maximum']]
        _tmodel._var_dict[varname].variable.lb = the_min
        _tmodel._var_dict[varname].variable.ub = the_max

    return _tmodel


def apply_directionality(tmodel,solution, inplace = True):
    """
    Takes a flux solution and transfers its reaction directionality as
    constraints for the cobra_model

    :param inplace:
    :param tmodel:
    :param solution:
    :return:
    """

    if inplace:
        _tmodel = tmodel
    else:
        _tmodel = tmodel.copy()

    for this_reaction in _tmodel.reactions:

        backward_use = _tmodel.backward_use_variable.get_by_id(this_reaction.id)
        forward_use = _tmodel.forward_use_variable.get_by_id(this_reaction.id)

        backward_use.variable.lb = round(solution.x_dict[backward_use.name])
        backward_use.variable.ub = round(solution.x_dict[backward_use.name])

        forward_use.variable.lb  = round(solution.x_dict[forward_use.name])
        forward_use.variable.ub  = round(solution.x_dict[forward_use.name])

    return _tmodel