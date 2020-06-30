"""
Named constants, or objects with values that are constant, used by other
modules.
"""

#  Universal gas constant [J mol^-1 K^-1]
GAS_CONST = 8.314
PR_ATM = 101325
KCAL_JL = 4184
HT_JL = 2625.5 * 1000
CAL_JL = 4.18
#  The cooling rate in *C/min to use if a cool down period is specified
COOL_RATE = -13
#  Tolerances for ODE solver
ABSOLUTE_TOLERANCE = 1E-11
RELATIVE_TOLERANCE = 1E-9

#  MW Dictionary values are:
# [0] MW
# [1] phase (g,s,lt,t,char,H20,CO,CO2)
# The rest of the columns are information that is only needed for tar
# components, so it may be blank in species of other phases
# [2] phenolic family (phenol/syringol)
# [3] number of [C,H,O] in the molecule
# [4] number of each type of C-bond/functional group in each molecule
# [carbonyl, aromatic C-O, aromatic C-C, aromatic C-H, aliphatic C-O,
# Methoxyl-aromatic, aliphatic C-C]
# .
# You can make changes to this dictionary to
# adjust which fraction each species
# is attributed to when lumping
# (commonly changed might be methanol, ETOH,
# ALD3, so those are located at the top).
# .
# The phenolic family entry indicates whether a tar species belongs to the
# phenol family or the syringol family (has methoxyl groups).
