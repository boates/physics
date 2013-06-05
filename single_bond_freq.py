#!/usr/bin/env python
# test for bond freq

import pylab

def main():

    f = open('TRAJEC.xyz','r')
    t = []
    dr = []
    lines = f.readlines()
    f.close()

    dr = []
    out = open('f_bond.dat','w')

    for i in range(len(lines)/4):
        natom = int(lines.pop(0))
        out.write(str(int(lines.pop(0)))+'    ')
        typat1,x1,y1,z1 = lines.pop(0).split()
        typat2,x2,y2,z2 = lines.pop(0).split()

        r = ( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )**0.5
        out.write(str(abs(r))+'\n')
        dr.append(abs(r))

        for i in range(natom-2): lines.pop(0)

    out.close()

    pxx, f = pylab.psd(dr,NFFT=32768)
    pylab.savefig('psd.png')

    out2 = open('psd.dat','w')
    for i in range(len(pxx)):
        out2.write(str(pxx[i])+'    '+str(f[i])+'\n')

    out2.close()


if __name__ == '__main__':
    main()
