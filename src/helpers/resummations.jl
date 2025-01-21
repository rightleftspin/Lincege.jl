"""
Resummation functions that take in thermodynamic properties
and return the resummed versions at various levels of correction.

NOTE: Please make sure to use at least two resummation techniques as
they will diverge at a further point than the bare NLCE summation.
This divergence between the resummation should be treated as the new
point to which the data is reasonable.
"""

using LinearAlgebra



