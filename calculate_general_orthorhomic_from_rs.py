#!/usr/bin/env python

# backlib.py

import math

# constants

a0 = 5.291772108E-09  # Bohr radius in cm
pi = math.pi

# inputs

Rs = input('For what value of Rs would you like the box?   ')
print 'ax/ax:   1.0'
ay_over_ax = input('ay/ax:   ')
az_over_ax = input('az/ax:   ')
bx_over_ax = input('bx/ax:   ')
by_over_ax = input('by/ax:   ')
bz_over_ax = input('bz/ax:   ')
cx_over_ax = input('cx/ax:   ')
cy_over_ax = input('cy/ax:   ')
cz_over_ax = input('cz/ax:   ')
N = input('How many particles are there in your box?   ')

volume = (1.0/3)*(4*N*pi)*(Rs*a0)**3
vol_term = (by_over_ax*cz_over_ax-bz_over_ax*cy_over_ax) + \
            ay_over_ax*(bz_over_ax*cy_over_ax-bx_over_ax*cz_over_ax) + \
            az_over_ax*(bx_over_ax*cy_over_ax-by_over_ax*cx_over_ax)
ax = ( volume / vol_term )**(1.0/3.0)
ay = ax * ay_over_ax
az = ax * az_over_ax
bx = ax * bx_over_ax
by = ax * by_over_ax
bz = ax * bz_over_ax
cx = ax * cx_over_ax
cy = ax * cy_over_ax
cz = ax * cz_over_ax


print '\n N =\t',N
print 'rs =\t',Rs
print 'V =\t',volume / 1.0E-8**3, 'angstroms^3'
print 'ax =\t',ax / 1.0E-8, 'angstroms'
print 'ay =\t',ay / 1.0E-8, 'angstroms'
print 'az =\t',az / 1.0E-8, 'angstroms'
print 'bx =\t',bx / 1.0E-8, 'angstroms'
print 'by =\t',by / 1.0E-8, 'angstroms'
print 'bz =\t',bz / 1.0E-8, 'angstroms'
print 'cx =\t',cx / 1.0E-8, 'angstroms'
print 'cy =\t',cy / 1.0E-8, 'angstroms'
print 'cz =\t',cz / 1.0E-8, 'angstroms'
