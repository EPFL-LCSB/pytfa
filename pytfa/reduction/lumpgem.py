# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Genome-scale model reduction

.. moduleauthor:: pyTFA team

Model reduction implementation of

Ataman M, Hatzimanikatis V (2017).
lumpGEM: Systematic generation of subnetworks and elementally balanced lumped
reactions for the biosynthesis of target metabolites.
PLoS Computational Biology 13(7): e1005513.
https://doi.org/10.1371/journal.pcbi.1005513
"""

from cobra import Reaction

class LumpedReaction(Reaction):

    def __init__(self,**kwargs):
        Reaction.__init__(**kwargs)
