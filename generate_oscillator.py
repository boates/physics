#!/usr/bin/env python
"""
Create an xyz file of a oscillating Nitrogen molecule.
"""
import Numeric

def main():

    x = Numeric.arange(2000, typecode=Numeric.Float)*0.1
    y = (Numeric.sin(x)+2.)/2.

    out = open('molecule.xyz','w')
    for i in range(len(y)):
        out.write('2\n'+str(i)+'\n')
        out.write('N '+str(y[i]/2.)+' 0.0 0.0\n')
        out.write('N '+str(-y[i]/2.)+' 0.0 0.0\n')
    out.close()

    out = open('bond_length.dat','w')
    for i in range(len(y)):
        out.write(str(x[i])+'   '+str(y[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
