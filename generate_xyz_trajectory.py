#!/usr/bin/env python
"""
Create an xyz file for one particle
"""
import os, sys, commands
import Numeric, math

def main():

    pi = math.pi
    t = Numeric.arange(10000,typecode=Numeric.Float)*0.01
    x = Numeric.sin(2*pi*1.1*t)
    y = Numeric.sin(2*pi*2.0*t)
    z = Numeric.sin(2*pi*2.7*t)
    
    out = open('ATOM.xyz','w')
    for i in range(len(t)):
        out.write('1\n'+str(i)+'\n')
        out.write('N '+str(x[i])+' '+str(y[i])+' '+str(z[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
