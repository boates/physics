#!/usr/bin/env python
"""
make lines for phase diagram based on clausius-clapeyron slopes.
"""
import os, sys, commands, glob
global dP
dP = 2  # GPa

def main():

    try:
        f = open('clausius_clapeyron_slopes.dat','r')
        cc = f.readlines()
        if cc[0][0] == '#':
            cc.pop(0)
        f.close()
        g = open('TRANSITION_64.dat','r')
        tr = g.readlines()
        if tr[0][0] == '#':
            tr.pop(0)
        g.close()
    except:
        print '\n make sure specifically named dat files are present. \n'
        sys.exit(0)

    Tc, slopes, Tt, P = [], [], [], []
    for line in cc:
        Tc.append(float(line.split()[0]))
        slopes.append(float(line.split()[1]))
    for line in tr:
        Tt.append(float(line.split()[1]))
        P.append(float(line.split()[0]))

#    T = [t for t in Tc if t in Tt]
#
#    for i in range(len(slopes)):
#        if Tc[i] not in T:
#            Tc.pop(i)
#            slopes.pop(i)
#    for i in range(len(P)):
#        if Tt[i]

    Tt.pop(-1)
    P.pop(-1)

    for i in range(len(slopes)):
        f = open('cc_'+str(int(Tt[i]))+'.dat','w')
        f.write(str(P[i]-dP/2.0)+' '+str(Tt[i]-slopes[i]/(dP/2.0))+'\n')
        f.write(str(P[i]+dP/2.0)+' '+str(Tt[i]+slopes[i]/(dP/2.0))+'\n')
        f.close()


if __name__ == '__main__':
    main()
