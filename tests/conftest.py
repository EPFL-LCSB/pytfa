"""Custom collection of test files.

.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

"""

import sys


if sys.version_info < (3, 6):
    collect_ignore_glob = ["test_equilibrator.py"]
