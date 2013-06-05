#!/usr/bin/env python
"""
Convert data files from rs to volume per atom
"""
import os, sys, commands, glob
import Numeric

def main():

    try:
        f = open(sys.argv[1],'r')
#        line = f.readline()
        lines = f.readlines()
#        if '#' != line[0]:
#            lines = line + lines
        natom = int(sys.argv[2])
        nval = int(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+'  fname  natom  nval\n'
        sys.exit(0)

    rs, x = [], []
    for line in lines:
        row = line.split()
        rs.append(float(row[0]))
        x.append(float(row[1]))

    rs = Numeric.array(rs)

    volume = (4.0/3.0)*(natom*nval*Numeric.pi)*(rs*0.5291772108)**3
    volume = list(volume/natom)

    out = open(sys.argv[1]+'.V','w')
    for i in range(len(volume)):
        out.write(str(volume[i])+'  '+str(x[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
