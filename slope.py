#!/usr/bin/env python

import sys

T1 = float(sys.argv[1])
T2 = float(sys.argv[2])
P1 = float(sys.argv[3])
P2 = float(sys.argv[4])

H1 = float(input('H1: '))
H2 = float(input('H2: '))

Tslope = (H2-H1) / (T2-T1)
Pslope = (H2-H1) / (P2-P1)

print
print 'Tslope =', Tslope, 'H per K'
print 'Pslope =', Pslope, 'H per GPa'
print
print 'T_H =', T1 - H1/Tslope, 'K'
print 'P_H =', P1 - H1/Pslope, 'GPa'
print 
