#!/usr/bin/env python

# backlib.py

import math

# constants

kb = 1.3806503E-23
Ha_J = 4.35791998E-18

# inputs

numberOfParticles = input('Number of particles? ')

# ECLASSIC = input('Value of ECLASSIC, $5 (Ha)? ')
# EKS = input('Value of EKS, $4 (Ha)? ')
# KEIONS = ECLASSIC - EKS

KEIONS = input('ekin_ion? ')
temperature = (2.0/3.0)*(1/kb)*KEIONS*Ha_J/numberOfParticles

print '\n Temperature is ' + repr(temperature) + ' K.'
        
