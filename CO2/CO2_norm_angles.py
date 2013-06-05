#!/usr/bin/env python
"""
Rescale the angle probablities based on sin(180-angle)
"""
import os, sys, commands, math, Numeric

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
    except:
        print '\n usage: '+sys.argv[0]+' CO2_angles.dat\n'
        sys.exit(0)

    # Read in the data
    angles, probs = [], []
    for line in lines:
        row = line.split()
        angles.append(float(row[0]))
        probs.append(float(row[1]))

    # Do the rescaling
    for i in range(len(angles)):
        factor = ( (math.sin(angles[i]*(math.pi/180.0)))**2 )**0.5
        if factor < 1.0E-03: factor = 1.0E-03
        probs[i] /= factor

    norm = Numeric.sum( Numeric.array(probs) )

    # Write to output file
    out = open('CO2_angles.dat.norm','w')
    for i in range(len(angles)-1):
        out.write(str(angles[i])+'  '+str(probs[i]/norm)+'\n')
    out.close()


if __name__ == '__main__':
    main()
