#!/usr/bin/env python

# Basically, just don't want to enter in tethered atom indices by hand
# print out the tethering lines of the FIELD file only, then insert (emacs ctrl-x i)
# it into your FIELD file assumes "harm"

import os, sys, commands

def main():

    k = float(input('Please give k in amu/ps^2\n'))
    i = int(input('Please give the index of the first atom you want to tether\n'))
    f = int(input('Please give the index of the last atom you want to tether\n'))

    out = open('teth.dat','w')
    for j in range(f-i+1):
        out.write('harm    '+str(j+i)+'    '+str(k)+'\n')

    out.close()

if __name__ == '__main__':
    main()

    
