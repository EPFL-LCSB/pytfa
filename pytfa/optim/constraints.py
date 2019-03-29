# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Constraints declarations

"""

from ..utils.str import camel2underscores

###################################################
###                 CONSTRAINTS                 ###
###################################################


class GenericConstraint:
    """
    Class to represent a generic constraint. The purpose is that the interface
     is instantiated on initialization, to follow the type of interface used
     by the problem, and avoid incompatibilities in optlang

    Attributes:

        :id: Used for DictList comprehension. Usually points back at a
        enzyme or reaction id for ease of linking. Should be unique given
        a constraint type.
        :name: Should be a concatenation of the id and a prefix that is
        specific to the variable type. will be used to address the constraint at
        the solver level, and hence should be unique in the whole cobra_model
        :expr: the expression of the constraint (sympy.Expression subtype)
        :cobra_model: the cobra_model hook.
        :constraint: links directly to the cobra_model representation of tbe constraint
    """
    prefix = ''


    @property
    def __attrname__(self):
        """
        Name the attribute the instances will have
        Example: GenericConstraint -> generic_constraint
        :return:
        """
        return camel2underscores(self.__class__.__name__)

    def __init__(self, expr, id_='', model=None, hook = None, queue=False, **kwargs):
        """

        :param id_: will be used to identify the variable
            (name will be a concat of this and a prefix)
        :param model: the cobra.Model object
        :param queue: whether or not to queue the variable for update object
        :param kwargs: stuff you want to pass to the variable constructor
        """
        self.hook = hook
        self._id = id_
        self._model = model
        self.kwargs = kwargs
        self._name = self.make_name()
        self.get_interface(expr, queue)
        self.prefix = ''

    def get_interface(self, expr, queue):
        """
        Called upon completion of __init__, initializes the value of self.var,
        which is returned upon call, and stores the actual interfaced variable.

        :return: instance of Variable from the problem
        """
        if not self.name in self.model.constraints:
            constraint = self.model.problem.Constraint(expression = expr,
                                                       name = self.name,
                                                       **self.kwargs)
            if not queue:
                self.model.add_cons_vars(constraint, sloppy=self.model.sloppy)
            else:
                self.model._cons_queue.append(constraint)
        else:
            self.constraint = self.model.constraints.get(self.name)


    def make_name(self):
        """
        Needs to be overridden by the subclass, concats the id with a
         prefix

        :return: None
        """
        return self.prefix + self.id

    def change_expr(self, new_expr, sloppy=False):

        lb = self.constraint.lb
        ub = self.constraint.ub
        name = self.name

        # Remove former constraint to override it
        self.model.solver.remove(name)
        new_cons = self.model.solver.interface.Constraint(name = name,
                                                   expression = new_expr,
                                                   ub = ub,
                                                   lb = lb)
        # Add the new variant
        self.model.solver.add(new_cons, sloppy=sloppy)

    @property
    def expr(self):
        return self.constraint.expression

    @expr.setter
    def expr(self,value):
        self.constraint.expression = value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def id(self):
        """
        for cobra.thermo.DictList compatibility
        :return:
        """
        return self._id

    @property
    def constraint(self):
        return self.model.constraints[self.name]

    @constraint.setter
    def constraint(self, value):
        self.model.constraints[self.name] = value

    @property
    def model(self):
        return self._model

    def __repr__(self):
        return self.name + ': ' + self.constraint.expression.__repr__()


class ReactionConstraint(GenericConstraint):
    """
    Class to represent a variable attached to a reaction
    """

    def __init__(self, reaction, expr, **kwargs):
        model = reaction.model

        GenericConstraint.__init__(self,
                                   expr=expr,
                                   model=model,
                                   hook=reaction,
                                   **kwargs)


    @property
    def reaction(self):
        return self.hook

    @property
    def id(self):
        return self.reaction.id

    @property
    def model(self):
        return self.reaction.model

class MetaboliteConstraint(GenericConstraint):
    """
    Class to represent a variable attached to a enzyme
    """

    def __init__(self, metabolite, expr, **kwargs):
        model = metabolite.model

        GenericConstraint.__init__(self,
                                   expr=expr,
                                   model=model,
                                   hook=metabolite
                                   **kwargs)

    @property
    def metabolite(self):
        return self.hook

    @property
    def id(self):
        return self.metabolite.id

    @property
    def model(self):
        return self.metabolite.model

class NegativeDeltaG(ReactionConstraint):
    """
    Class to represent thermodynamics constraints.

    G: - DGR_rxn + DGoRerr_Rxn + RT * StoichCoefProd1 * LC_prod1
       + RT * StoichCoefProd2 * LC_prod2
       + RT * StoichCoefSub1 * LC_subs1
       + RT * StoichCoefSub2 * LC_subs2
       - ...
     = 0
    """

    prefix = 'G_'

class ForwardDeltaGCoupling(ReactionConstraint):
    """
    Class to represent thermodynamics coupling: DeltaG of reactions has to be
    DGR < 0 for the reaction to proceed forwards
    Looks like:
    FU_rxn: 1000 FU_rxn + DGR_rxn < 1000
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'FU_'

class BackwardDeltaGCoupling(ReactionConstraint):
    """
    Class to represent thermodynamics coupling: DeltaG of reactions has to be
    DGR > 0 for the reaction to proceed backwards
    Looks like:
    BU_rxn: 1000 BU_rxn - DGR_rxn < 1000
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'BU_'

class ForwardDirectionCoupling(ReactionConstraint):
    """
    Class to represent a forward directionality coupling with thermodynamics on
    reaction variables
    Looks like :
    UF_rxn: F_rxn - M FU_rxn < 0
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'UF_'


class BackwardDirectionCoupling(ReactionConstraint):
    """
    Class to represent a backward directionality coupling with thermodynamics on
    reaction variables
    Looks like :
    UR_rxn: R_rxn - M RU_rxn < 0
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'UR_'

class SimultaneousUse(ReactionConstraint):
    """
    Class to represent a simultaneous use constraint on reaction variables
    Looks like:
    SU_rxn: FU_rxn + BU_rxn <= 1
    """

    prefix = 'SU_'

class DisplacementCoupling(ReactionConstraint):
    """
    Class to represent the coupling to the thermodynamic displacement
    Looks like:
    Ln(Gamma) - (1/RT)*DGR_rxn = 0
    """

    prefix = 'DC_'

class ForbiddenProfile(GenericConstraint):
    """
    Class to represent a forbidden net flux directionality profile
    Looks like:
    FU_rxn_1 + BU_rxn_2 + ... + FU_rxn_n <= n-1
    """

    def __init__(self, model, expr, id_, **kwargs):

        GenericConstraint.__init__(self,
                                   id_=id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)

    prefix = 'FP_'
