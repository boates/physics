#!/usr/bin/env python

import numpy, math
from math import e

a0 = 60.8279910435
a1 = 0.0925493792456
a2 = -0.000214695155605

T0 = 11788.6
P0 = 326.56

Pmin = P0-10
Pmax = 650
Pstep = 1.0

P = numpy.arange( Pmin, Pmax+1, Pstep)
T = T0 * (1+(P-P0)/a0)**a1 * e**(-a2*(P-P0))

out = open('kechin_extend.dat','w')
for i in range(len(P)):
    out.write(str(P[i])+' '+str(T[i])+'\n')
out.close()
