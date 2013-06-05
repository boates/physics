#!/usr/bin/env python
"""
MC simulation from derivatives COMM 4202 w/ R. Nason
"""
import os, sys, commands
import random, numpy, pylab
global e
e = numpy.e

def main():

    # Define input parameters
    N   =  10000
    S   =  40.0
    K   =  40.0
    V   =  0.10
    T   =  0.50
    I   =  0.05
    t   =  0.0625
    KO  =  50.0
    KI  =  37.0

    # Calculate the prices from a random gaussian distribution
    P, random_list = [], []
    for i in range(N):
        P.append([S])
        for j in range(1,int(T/t)+1):
            RANDOM = random.normalvariate(0,V)
            random_list.append(RANDOM)
            PREV = P[i][j-1]
            P[i].append(PREV*( e**( t*(I - 0.5*V**2) + V*RANDOM*t**0.5 ) ) )
            if i == 0: print RANDOM, PREV, P[i][j]
        P[i].pop(0)
    P = numpy.array(P)

#    pylab.hist(random_list,1000)
#    pylab.show()

    # Initialize option arrays
    call = []
    put = []
    asian = []
    knock_out = []
    knock_in = []
    KIKO = []

    # Calculate the option P
    for i in range(len(P)):
        call.append(  max(0, P[i][-1] - K) )
        put.append(   max(0, K - P[i][-1]) )
        asian.append( max(0, numpy.average(P[i]) - K) )
    print call.count(0.0), put.count(0.0)

    # Determine the average option P
    callAvg, putAvg, asianAvg = 0.0, 0.0, 0.0
    for i in range(len(P)):
        callAvg  += call[i]  / len(call)
        putAvg   += put[i]   / len(put)
        asianAvg += asian[i] / len(asian)

    print 'call =', callAvg*e**(-I*T),'     put = ', putAvg*e**(-I*T),'    asian = ',asianAvg*e**(-I*T)

if __name__ == '__main__':
    main()
