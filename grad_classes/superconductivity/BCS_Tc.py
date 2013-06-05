#!/usr/bin/env python
"""
Approximate maximum values for Tc from BCS theory for a
range of N(0)*V values and given Debye temperature.
"""
import os, sys, Numeric

def main():

    # Retrieve user input
    try:
        TD = float(sys.argv[1])
    except:
        print '\n usage: '+sys.argv[0]+' Debye_temperature(K)\n'
        sys.exit(0)

    # Create array of possible N(0)*V values
    N0V = Numeric.arange(0.005,1.001,0.001,typecode=Numeric.Float)

    # Calculate Tc from unapproximated expression
    Tc = 1.14 * TD / ( Numeric.e**(1/N0V) - 1 )
    Tc2 =  Numeric.e**(1/N0V) - 1

    # Calculate in the weak coupling limit
    WCL_Tc = 1.14 * TD * Numeric.e**(-1/N0V)
    WCL_Tc2 = Numeric.e**(-1/N0V)

    # Calculate the percentage difference between the two
    P = ( Tc2 - WCL_Tc2 ) / Tc2 *100

    # Write the results to file
    out = open('BCS.dat','w')
    out.write('# N(0)*V, %_difference, Tc(K), WCL_Tc(K), done using Debye_temperature of '+str(TD)+'\n')
    for i in range(len(P)):
        out.write(str(N0V[i])+' '+str(P[i])+' '+str(Tc[i])+' '+str(WCL_Tc[i])+'\n')
    out.close()
        

if __name__ == '__main__':
    main()
