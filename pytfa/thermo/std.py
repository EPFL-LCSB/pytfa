# -*- coding: utf-8 -*-
"""Standard constants definitions
"""

from math import log

TEMPERATURE_0 = 298.15  # K
MIN_PH = 3
MAX_PH = 9

# Used to calculate the transformed Gibbs energy of formation of specie with
# given pH and ionic strength using formula given by Goldberg and Tewari, 1991
# equation 4.4-10 in Alberty's book
DEBYE_HUCKEL_B_0 = 1.6
DEBYE_HUCKEL_A = 1.17582 / log(10)

A_LOT = 5000
A_LITTLE = 0.5
A_COUPLE = 2.5 # A couple is usually considered to be 2 or 3
MANY = 100