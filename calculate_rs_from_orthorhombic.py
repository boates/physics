#!/usr/bin/env python

# backlib.py
import math

def main():

    ax = float(input('ax (angstroms): ')) / 0.5291772
    ay = float(input('ay (angstroms): ')) / 0.5291772
    az = float(input('az (angstroms): ')) / 0.5291772
    N = int(input('How many particles are there in your box? '))

    volume = ax*ay*az
    rs = ( (3*volume)/(4*math.pi*N) )**(1.0/3)

    print 'ax:',ax*0.5291772,'angstroms'
    print 'ay:',ay*0.5291772,'angstroms'
    print 'az:',az*0.5291772,'angstroms'
    print 'volume:',volume*0.5291772**3,'angstroms^3'
    print 'rs:',rs

if __name__ == '__main__':
    main()
