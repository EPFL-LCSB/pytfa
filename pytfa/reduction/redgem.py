# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Genome-scale model reduction

.. moduleauthor:: pyTFA team

Model reduction implementation of

Ataman M, Hernandez Gardiol DF, Fengos G, Hatzimanikatis V (2017)
redGEM: Systematic reduction and analysis of genome-scale metabolic
reconstructions for development of consistent thermo metabolic models.
PLoS Computational Biology 13(7): e1005444.
https://doi.org/10.1371/journal.pcbi.1005444
"""


from ..thermo.tmodel import ThermoModel

def reduce(tmodel, reduction_parameters):
    reduction = ReducedModel(tmodel, reduction_parameters)

    #------------------- To complete -----------------

    return reduction

class ReducedModel(ThermoModel):
    """
    A class representing a thermo model that has been reduced.
    Such a class is needed as reduced models will have lumps, and meta information
    on the parameters of reduction.
    """

    def __init__(self, tmodel, reduction_parameters):

        # Import the tmodel
        self.import_tmodel(tmodel)

        self.reduction_parameters = reduction_parameters


    def import_tmodel(self,tmodel):
        """
        Imports the attributes of tmodel into this instance of ReducedModel
        :param tmodel:
        :return:
        """

        # Ugly hack
        tcopy = tmodel.copy()
        self.__dict__ = tcopy.__dict__
