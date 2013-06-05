#!/usr/bin/env python
"""
Author: Brian Boates
Date: May 5, 2009

Contains useful functions encountered in analysis codes
"""
from math import pi

def pbc_round(input_value):
    i = int(input_value)
    if (abs(input_value-i) >= 0.5):
        if (input_value > 0): i+= 1
        if (input_value < 0): i-= 1
    return i
            
def volume_from_rs(rs,Nel):
    """
    Return cell volume in angstroms^3.
    rs: density parameter for system
    Nel: number of valence electrons in system
    """
    a0 = 0.5291772   # Bohr radius (angstroms/bohr)
    volume = (4.0*pi/3.0)*Nel * (rs*a0)**3

    return volume
