# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

JSON serialization
"""

import json

import numpy

from .dict import model_from_dict, model_to_dict


class MyEncoder(json.JSONEncoder):
    """
    We define an encoder that takes care of the serialization of numpy types,
    which are not handled by json by default
    """
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)


def check_json_extension(filepath):
    if not filepath.endswith('.json'):
        filepath += '.json'
    return filepath

def save_json_model(model, filepath):

    filepath = check_json_extension(filepath)
    obj = model_to_dict(model)

    with open(filepath, 'w') as fid:
        json.dump(obj, fid, cls=MyEncoder)


def load_json_model(filepath):

    filepath = check_json_extension(filepath)
    with open(filepath, 'r') as fid:
        obj = json.load(fid)

    model = model_from_dict(obj)
    return model
