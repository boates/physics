#!/usr/bin/env python

# backlib.py

import math

# constants

a0 = 5.291772108E-09  # Bohr radius in cm
pi = math.pi

# inputs

Rs = input('For what value of Rs would you like the box?   ')
ay_over_ax = input('ay/ax:   ')
az_over_ax = input('az/ax:   ')
N = input('How many particles are there in your box?   ')

volume = (1.0/3)*(4*N*pi)*(Rs*a0)**3
ax = ( volume / (ay_over_ax*az_over_ax) )**(1.0/3.0) 
ay = ax * ay_over_ax
az = ax * az_over_ax


print '\n N =\t',N
print 'rs =\t',Rs
print 'V =\t',volume / 1.0E-8**3, 'angstroms^3'
print 'ax =\t',ax / 1.0E-8, 'angstroms'
print 'ay =\t',ay / 1.0E-8, 'angstroms'
print 'az =\t',az / 1.0E-8, 'angstroms'
