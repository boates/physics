#!/usr/bin/env python
"""
Nonlinear leasr-squares curve fitting for EOS, using either
a Birch-Murnaghan fit of 2nd or 3rd order.
"""
import os, sys, commands, glob
import Numeric, math
from scipy import optimize
from scipy import linalg

def BM_2nd_order(rs,P,verbose=False):

    def resid_2nd(a,rs,P):

        return P - (3*a[0]/2.0)*( (a[1]/rs)**7. - (a[1]/rs)**5 )

    a0 = Numeric.array([-5.0,-1.0])

    a, mesg = optimize.leastsq(resid_2nd,a0,args=(rs,P),maxfev=100000)

    P_fit = (3*a[0]/2.0)*( (a[1]/rs)**7. - (a[1]/rs)**5 )

    if verbose:
        print a

    return P_fit

def BM_3rd_order(rs,P,verbose=False):

    V = rs**3*4./3.*math.pi*320.0

    def resid_3rd(a,rs,P):

        return P - (3*a[0]/2.0)*( (a[1]/V)**(7./3.) - (a[1]/V)**(5./3.) )*( 1+3/4.*(a[2]-4.)*((a[1]/V)**(2./3.)-1) )

    a0 = Numeric.array([1.0,2000.0,1.0])
    
    a, mesg = optimize.leastsq(resid_3rd,a0,args=(rs,P),maxfev=100000)

    P_fit = (3*a[0]/2.0)*( (a[1]/V)**(7./3.) - (a[1]/V)**(5./3.) )*( 1+3/4.*(a[2]-4.)*((a[1]/V)**(2./3.)-1) )

    if verbose:
        print '\nBirch-Murnaghan parameters are:', a[0],a[1],a[2],'\n'

    return P_fit

def main():

    # Retrieve data from user and format
    try:
        f = open(sys.argv[1])
        header = f.readline()
        if '#' != header.strip()[0]:
            f.close()
            f = open(sys.argv[1])
        lines = f.readlines()
        rs, P = [], []
        for line in lines:
            rs.append(float(line.split()[0]))
            P.append(float(line.split()[1]))
        rs = Numeric.array(rs)
        P = Numeric.array(P)
    except:
        print '\n usage: '+sys.argv[0]+' file_with_rs_and_P_as_1st_two_columns.dat \n'
        sys.exit(0)

#    P_new = BM_2nd_order(rs,P,verbose=True)
    P_new = BM_3rd_order(rs,P,verbose=True)

    # Write the fit to file
    out = open('BM_fit.dat','w')
    for i in range(len(rs)):
        out.write(str(rs[i])+'  '+str(P_new[i])+'\n')
    out.close()
    

if __name__ == '__main__':
    main()
