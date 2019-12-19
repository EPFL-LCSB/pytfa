"""Custom collection of test files.

.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

"""

import pytest
import sys


if sys.version_info < (3, 6):
    pytest.xfail("equilibrator-cache requires Python version >= 3.6")
