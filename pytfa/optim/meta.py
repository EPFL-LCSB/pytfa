# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Metaclass declarations to force the definition of prefixes in GenericVariable
and GeneriConstraint subclasses

Based on SethMMorton's answer on StackOverflow
https://stackoverflow.com/questions/45248243/most-pythonic-way-to-declare-an-abstract-class-property
https://stackoverflow.com/a/45250114
"""
from abc import ABCMeta

class RequirePrefixMeta(type):
    """Metaclass that enforces child classes define prefix."""

    def __init__(cls, name, bases, attrs):
        # Skip the check if there are no parent classes,
        # which allows base classes to not define prefix.
        if not bases:
            return
        if attrs.get('prefix', NotImplemented) is NotImplemented:
            # Choose your favorite exception.
            raise NotImplementedError('You forgot to define a prefix for your class ;)')

ABCRequirePrefixMeta = type('ABCRequirePrefixMeta',
                            (ABCMeta, RequirePrefixMeta), {})
