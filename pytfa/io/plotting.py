# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Plotting results

"""
import numpy as np
import pandas as pd

from bokeh.models import ColumnDataSource, DataRange1d, LinearAxis, \
    Grid,CategoricalTicker, FuncTickFormatter
from bokeh.models.glyphs import HBar
from bokeh.io import show
from bokeh.plotting import figure

from ..optim.variables import ThermoDisplacement

def plot_fva_tva_comparison(fva,tva):

    all_va = pd.merge(fva,tva, left_index = True, right_index=True)

    all_va['y'] = range(len(all_va))

    source = ColumnDataSource(all_va)

    xdr = DataRange1d()
    ydr = DataRange1d()

    plot = figure(title=None, x_range=xdr, y_range=ydr, plot_width=600,
        plot_height=1000, h_symmetry=False, v_symmetry=False, min_border=0)

    glyph_fva = HBar(y="y", right="maximum_x", left="minimum_x", height=0.5,
                  fill_color="#b3de69", fill_alpha = 0.5)
    glyph_tva = HBar(y="y", right="maximum_y", left="minimum_y", height=0.8,
                  fill_color="#2e86c1", fill_alpha = 0.5)
    plot.add_glyph(source, glyph_tva)
    plot.add_glyph(source, glyph_fva)


    # Fix ticks
    label_dict = {}
    for i, s in enumerate(tva.index):
        label_dict[i] = s
    plot.yaxis.formatter = FuncTickFormatter(code="""
                                            var labels = %s;
                                            return labels[tick];
                                        """ % label_dict)

    plot.yaxis.ticker = [x for x in range(len(tva))]

    return plot

def plot_thermo_displacement_histogram(tmodel,solution = None):
    """
    Plot a histogram of the thermodynamic displacement. if no solution is
    provided, will look at the cobra_model's own solution

    :param tmodel:
    :param solution:
    :return:
    """

    if solution is None:
        values = tmodel.get_primal(ThermoDisplacement)
    else:
        thermo_displacements = tmodel.get_variables_of_type(ThermoDisplacement)
        varnames = [x.name for x in thermo_displacements]

        values = solution.x_dict[varnames]

    p1 = plot_histogram(values)


    p1.legend.location = "center_right"
    p1.legend.background_fill_color = "darkgrey"
    p1.xaxis.axis_label = 'ln(Γ)'
    p1.yaxis.axis_label = 'count'
    p1.title = "Histogram: Thermodynamic displacement"

    return p1

def plot_histogram(values, **kwargs):
    """
    Convenience function. Plots a histogram of flat 1D data.

    :param values:
    :return:
    """


    hist, edges = np.histogram(values, **kwargs)

    p1 = figure(tools="save")

    p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color="#036564", line_color="#033649")


    return p1