# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Model class
"""

from abc import ABC, abstractmethod
from collections import defaultdict

import pandas as pd
from numpy import empty
import optlang
from optlang.exceptions import SolverError
from cobra import DictList, Model
from cobra.core.solution import Solution

from ..utils.str import camel2underscores
from ..optim.variables import GenericVariable, ReactionVariable, MetaboliteVariable
from ..optim.constraints import ReactionConstraint, MetaboliteConstraint
from ..optim.utils import get_primal, get_all_subclasses

import time

def timeit(method):
    """
    Adapted from Andreas Jung's `blog
    <https://www.zopyx.com/andreas-jung/contents/a-python-decorator-for-measuring-the-execution-time-of-methods>`_

    :param method: The method to time
    :return:
    """


    def timed(self, *args, **kw):
        ts = time.time()
        result = method(self, *args, **kw)
        te = time.time()

        message = '%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts)

        try:
            self.logger.debug(message)
        except AttributeError:
            print(message)
        return result

    return timed

class LCSBModel(ABC):

    # @abstractmethod
    def __init__(self, model, name, sloppy=False):

        """
        Very much model specific
        """

        Model.__init__(self, model.copy(), name)

        self._cons_queue = list()
        self._var_queue = list()

        self._var_dict = dict()
        self._cons_dict = dict()

        self.sloppy=sloppy


    @abstractmethod
    def copy(self):
        """
        Needs to be reimplemented, as our objects have complicated hierarchy
        :return:
        """

    def print_info(self):
        """
        Print information and counts for the cobra_model
        :return:
        """

        n_metabolites = len(self.metabolites)
        n_reactions = len(self.reactions)
        n_constraints = len(self.constraints)
        n_variables = len(self.variables)

        info = pd.DataFrame(columns=['value'])
        info.loc['name'] = self.name
        info.loc['description'] = self.description
        info.loc['num constraints'] = n_constraints
        info.loc['num variables'] = n_variables
        info.loc['num metabolites'] = n_metabolites
        info.loc['num reactions'] = n_reactions
        info.index.name = 'key'

        print(info)

    def add_variable(self, kind, hook, queue=False, **kwargs):
        """ Add a new variable to a COBRApy cobra_model.

        :param kind:
        :param string,cobra.Reaction hook: Either a string representing the name
            of the variable to add to the cobra_model, or a reaction object if the
            kind allows it

        :returns: The created variable
        :rtype: optlang.interface.Variable

        """

        # Initialisation links to the cobra_model
        var = kind(hook,
                   # lb=lower_bound if lower_bound != float('-inf') else None,
                   # ub=upper_bound if upper_bound != float('inf') else None,
                   queue=queue,
                   **kwargs)

        self._var_dict[var.name] = var
        self.logger.debug('Added variable: {}'.format(var.name))
        # self.add_cons_vars(var.variable)

        return var

    def add_constraint(self, kind, hook, expr, queue=False,**kwargs):
        """ Add a new constraint to a COBRApy cobra_model

        :param kind:
        :param string,cobra.Reaction hook: Either a string representing the name
            of the variable to add to the cobra_model, or a reaction object if the
            kind allows it
        :param sympy.thermo.expr.Expr expr: The expression of the constraint

        :returns: The created constraint
        :rtype: optlang.interface.Constraint

        """

        if isinstance(expr, GenericVariable):
            # make sure we actually pass the optlang variable
            expr = expr.variable

        # Initialisation links to the cobra_model
        cons = kind(hook, expr, # problem = self.problem,
                    # lb=lower_bound if lower_bound != float('-inf') else None,
                    # ub=upper_bound if upper_bound != float('inf') else None,
                    queue=queue,
                    **kwargs)
        self._cons_dict[cons.name] = cons
        self.logger.debug('Added constraint: {}'.format(cons.name))
        # self.add_cons_vars(cons.constraint)

        return cons

    def remove_reactions(self, reactions, remove_orphans=False):
        # Remove the constraints and variables associated to these reactions
        all_cons_subclasses = get_all_subclasses(ReactionConstraint)
        all_var_subclasses = get_all_subclasses(ReactionVariable)

        self._remove_associated_consvar(all_cons_subclasses, all_var_subclasses,
                                        reactions)

        Model.remove_reactions(self,reactions,remove_orphans)

    def remove_metabolites(self, metabolite_list, destructive=False):
        # Remove the constraints and variables associated to these reactions
        all_cons_subclasses = get_all_subclasses(MetaboliteConstraint)
        all_var_subclasses = get_all_subclasses(MetaboliteVariable)

        self._remove_associated_consvar(all_cons_subclasses, all_var_subclasses,
                                        metabolite_list)

        Model.remove_metabolites(self, metabolite_list, destructive)

    def _remove_associated_consvar(self, all_cons_subclasses, all_var_subclasses,
                                   collection):
        """
        Removes both the constraints and variables associated to an element,
        as long as it was used as a hook in the cons/var declaration.
        For example, upon removing a reaction, also removes its associated
        deltaG variables and coupling constraints
        """

        if not hasattr(collection, '__iter__'):
            collection = [collection]

        strfy = lambda x:x if isinstance(x, str) else x.id

        for cons_type in all_cons_subclasses:
            for element in collection:
                try:
                    cons = self._cons_kinds[cons_type.__name__].get_by_id(strfy(element))
                    self.remove_constraint(cons)
                except KeyError as e:
                    pass
        for var_type in all_var_subclasses:
            for element in collection:
                try:
                    var = self._var_kinds[var_type.__name__].get_by_id(strfy(element))
                    self.remove_variable(var)
                except KeyError as e:
                    pass


    def remove_variable(self, var):
        """
        Removes a variable

        :param var:
        :return:
        """
        # Get the pytfa var object if an optlang variable is passed
        if isinstance(var,optlang.Variable):
            var = self._var_dict[var.name]

        self._var_dict.pop(var.name)
        self.remove_cons_vars(var.variable)
        self.logger.debug('Removed variable {}'.format(var.name))

    def remove_constraint(self, cons):
        """
        Removes a constraint

        :param cons:
        :return:
        """
        # Get the pytfa var object if an optlang variable is passed
        if isinstance(cons,optlang.Constraint):
            cons = self._cons_dict[cons.name]

        self._cons_dict.pop(cons.name)
        self.remove_cons_vars(cons.constraint)
        self.logger.debug('Removed constraint {}'.format(cons.name))

    def _push_queue(self):
        """
        updates the constraints and variables of the model with what's in the
        queue
        :return:
        """

        self.add_cons_vars(self._var_queue, sloppy=self.sloppy)
        self.add_cons_vars(self._cons_queue, sloppy = self.sloppy)

        if len(self._var_queue) > 0:
            self.regenerate_variables()
        if len(self._cons_queue) > 0:
            self.regenerate_constraints()

        self._var_queue = list()
        self._cons_queue = list()


    def regenerate_variables(self):
        """
        Generates references to the cobra_model's constraints in self._var_dict
        as tab-searchable attributes of the thermo cobra_model
        :return:
        """

        # Let us not forget to remove fields that might be empty by now
        if hasattr(self, '_var_kinds'):
            for k in self._var_kinds:
                attrname = camel2underscores(k)
                try:
                    delattr(self, attrname)
                except AttributeError:
                    pass # The attribute may not have been set up yet

        _var_kinds = defaultdict(DictList)
        for k, v in self._var_dict.items():
            _var_kinds[v.__class__.__name__].append(v)

        for k in _var_kinds:
            attrname = camel2underscores(k)
            setattr(self, attrname, _var_kinds[k])

        self._var_kinds = _var_kinds

    def regenerate_constraints(self):
        """
        Generates references to the cobra_model's constraints in self._cons_dict
        as tab-searchable attributes of the thermo cobra_model
        :return:
        """

        # Let us not forget to remove fields that migh be empty by now
        if hasattr(self, '_cons_kinds'):
            for k in self._cons_kinds:
                attrname = camel2underscores(k)
                try:
                    delattr(self, attrname)
                except AttributeError:
                    pass # The attribute may not have been set up yet

        _cons_kinds = defaultdict(DictList)

        for k, v in self._cons_dict.items():
            _cons_kinds[v.__class__.__name__].append(v)

        for k in _cons_kinds:
            attrname = camel2underscores(k)
            setattr(self, attrname, _cons_kinds[k])

        self._cons_kinds = _cons_kinds

    def repair(self):
        """
        Updates references to variables and constraints
        :return:
        """
        # self.add_cons_vars([x.constraint for x in self._cons_dict.values()])
        # self.add_cons_vars([x.variable for x in self._var_dict.values()])
        self._push_queue()
        Model.repair(self)
        self.regenerate_constraints()
        self.regenerate_variables()

    def get_primal(self, vartype, index_by_reactions=False):
        """
        Returns the primal value of the cobra_model for variables of a given type

        :param index_by_reactions:
        :param vartype: Class of variable. Ex: pytfa.optim.variables.ThermoDisplacement
        :return:
        """
        return get_primal(self, vartype, index_by_reactions)

    def get_solution(self):
        """
        Overrides the cobra.thermo.solution method, to also get the supplementary
        variables we added to the cobra_model

        *   :code:`solution.fluxes` in `cobrapy` is a transformed version of the solver
            output, as it actually calculates the _net_ flux of each reaction by
            substracting the reverse variable value to the forward variable value.
            This should be used anytime one needs the actual flux value

        *   :code:`solution.raw` is a clear copy of the solver output. From there one
            can access the value at solution for all the variables of the problem.
            However, looking for a reaction ID in there will only give the
            _forward_ flux. This should be used for any other variable than fluxes.

        *   :code:`solution.values` yields variables multiplied by their scaling factor
            (1 by default). Useful if you operated scaling on your equations for
            numerical reasons. This does _not_ include fluxes

        :return:
        """
        objective_value = self.solver.objective.value
        status = self.solver.status
        variables = pd.Series(data=self.solver.primal_values)

        fluxes = empty(len(self.reactions))
        rxn_index = list()
        var_primals = self.solver.primal_values

        for (i, rxn) in enumerate(self.reactions):
            rxn_index.append(rxn.id)
            fluxes[i] = var_primals[rxn.id] - var_primals[rxn.reverse_id]

        fluxes = pd.Series(index=rxn_index, data=fluxes, name="fluxes")

        solution = Solution(objective_value=objective_value, status=status,
                            fluxes=fluxes)

        self.solution = solution

        self.solution.raw = variables

        self.\
            solution.values = pd.DataFrame.from_dict({k:v.unscaled
                                                 for k,v in self._var_dict.items()},
                                                 orient = 'index')

        return solution

    def optimize(self, objective_sense=None, **kwargs):
        """
        Call the Model.optimize function (which really is but an interface to the
        solver's. Catches SolverError in the case of no solutions. Passes down
        supplementary keyword arguments (see cobra.thermo.Model.optimize)
        :type objective_sense: 'min' or 'max'
        """

        if objective_sense:
            self.objective.direction = objective_sense

        try:
            # self._hidden_optimize_call(kwargs)
            Model.optimize(self, **kwargs)
            solution = self.get_solution()
            self.solution = solution
            return solution
        except SolverError as SE:
            status = self.solver.status
            self.logger.error(SE)
            self.logger.warning('Solver status: {}'.format(status))
            raise (SE)

    # @timeit
    # def _hidden_optimize_call(self, kwargs):
    #     return Model.optimize(self, **kwargs)

    @timeit
    def slim_optimize(self, *args, **kwargs):
        return Model.slim_optimize(self, *args, **kwargs)

    def get_constraints_of_type(self, constraint_type):
        """
        Convenience function that takes as input a constraint class and returns
        all its instances within the cobra_model

        :param constraint_type:
        :return:
        """
        if isinstance(constraint_type,str):
            constraint_key = constraint_type
        else:
            #it is a class
            constraint_key = constraint_type.__name__
        return self._cons_kinds[constraint_key]

    def get_variables_of_type(self, variable_type):
        """
        Convenience function that takes as input a variable class and returns
        all its instances within the cobra_model

        :param variable_type:
        :return:
        """
        if isinstance(variable_type,str):
            variable_key = variable_type
        else:
            #it is a class
            variable_key = variable_type.__name__
        return self._var_kinds[variable_key]
