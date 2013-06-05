#!/usr/bin/env python

import math

# constants

kb = 1.3806503e-23
Ha = 4.3597482e-18 # this is from gaussian
eV = 27.2117
# inputs

temperature = input('Temperature (K)? ')

energy = (temperature*kb)/Ha 

print '\n Energy is ' + repr(energy) + ' Ha.'
print '\n Energy is ' + repr(energy*eV) + ' eV.'
