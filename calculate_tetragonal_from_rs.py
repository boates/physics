#!/usr/bin/env python

# backlib.py
import math

def main():
    # constants
    a0 = 5.291772108E-09  # Bohr radius in cm
    pi = math.pi

    # inputs
    Rs = input('For what value of Rs would you like the box? ')
    N = input('How many particles are there in your box? ')
    ratio = input('What is the c/a ratio? ')

    volume = (1.0/3)*(4*N*pi)*(Rs*a0)**3
    a = (volume*6.7483346e+24/ratio)**(1.0/3)
    c = a*ratio

    print 'rs:',Rs
    print 'ratio:',ratio
    print 'a:',a,'bohr'
    print 'c:',c,'bohr'
    print 'volume:',volume*6.7483346e+24,'bohr^3'
    print 'density:',N/(volume*6.7483346e+24),'particles/bohr^3'

if __name__ == '__main__':
    main()
