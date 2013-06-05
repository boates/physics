#!/usr/bin/env python
"""
Nonlinear leasr-squares curve fitting for
melting curves using Kechin fit
"""
import os, sys, commands, glob
import Numeric, math
from scipy import optimize
from scipy import linalg

def kechin(P,T,P0,T0,verbose=False):

    def resid(a,P,T,P0,T0):
        return T - T0*(1.0 + (P-P0)/a[0])**a[1] * math.e**(-a[2]*(P-P0))

#    a0 = Numeric.array([0.0001,10.0,0.0001])
    a0 = Numeric.array([0.1,1.0,0.1])
    a, mesg = optimize.leastsq(resid,a0,args=(P,T,P0,T0),maxfev=100000)

    if verbose:
        print '\n Kechin fitting parameters are:'
        print ' a='+str(a[0])+', b='+str(a[1])+', c='+str(a[2])+', T0='+str(T0)+', P0='+str(P0)+'\n'

    return a

def main():

    # Retrieve data from user and format
    try:
        f = open(sys.argv[1])
        lines = f.readlines()
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' P_melt-T_melt.dat \n'
        sys.exit(0)

    P, T = [], []
    for line in lines:
        P.append(float(line.split()[0]))
        T.append(float(line.split()[1]))
    P = Numeric.array(P)
    T = Numeric.array(T)

    # Define P0 and T0
    P0, T0 = P[0], T[0]

    # Get the fitting parameters
    a = kechin(P,T,P0,T0,verbose=True)

    # Create a new data set from the fit
    P_fit = Numeric.arange(P[0],P[-1],1.0)
    T_fit = T0*(1.0 + (P_fit-P0)/a[0])**a[1] * math.e**(-a[2]*(P_fit-P0))

    # Write the fit to file
    out = open('kechin.fit','w')
    for i in range(len(P_fit)):
        out.write(str(P_fit[i])+'  '+str(T_fit[i])+'\n')
    out.close()
    

if __name__ == '__main__':
    main()
