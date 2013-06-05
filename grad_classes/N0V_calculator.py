#!/usr/bin/env python
from math import pi

# Define constants
m    = 9.10938188e-31   # kg
hbar = 4.13566733e-15   # eV s
x    = 3.56#e-03         # eV
n    = 1.0e+26          # N/V, electron density in m^-3

# Convert to Joules
x    = x * 1.6021765e-19
hbar = hbar * 1.6021765e-19

# Calculate N(0)*V
N0V = ( x * m ) / ( hbar**2 * (3*pi**2*n)**(2/3.) ) * 0.88      # dimensionless
Ef  = ( hbar**2/(2*m) ) * (3*pi**2*n)**(2/3.) / 1.6021765e-19   # in eV

print '\n N(0)*V =', N0V,'\n'
print ' Ef =', Ef,'eV\n'

