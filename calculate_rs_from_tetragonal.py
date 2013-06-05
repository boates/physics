#!/usr/bin/env python

# backlib.py
import math

def main():
    # inputs
    a = float(input('What is the value of a (bohr) for your tetragonal box? '))
    ratio = float(input('What is the c/a ratio? '))
    N = int(input('How many particles are there in your box? '))

    volume = a**3.0 * ratio
    length = volume**(1.0/3)
    rs = ( (3*volume)/(4*math.pi*N) )**(1.0/3)
    c = a*ratio

    print 'a:',a,'bohr'
    print 'c:',c,'bohr'
    print 'ratio:',ratio
    print 'volume:',volume,'bohr^3'
    print 'rs:',rs

if __name__ == '__main__':
    main()
