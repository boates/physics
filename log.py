#!/usr/bin/env python

import math, numpy

t = range(1,1001)
x = [0.001, 0.005, 0.010, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09] + [i/10.0 for i in t]

log2  = open('log2.dat' ,'w')
log5  = open('log5.dat' ,'w')
log10 = open('log10.dat','w')

for v in x:

    log2.write(str(v) +" "+str( math.log(v,2)  )+"\n")
    log5.write(str(v) +" "+str( math.log(v,5)  )+"\n")
    log10.write(str(v)+" "+str( math.log(v,10) )+"\n")
    
log2.close()
log5.close()
log10.close()
