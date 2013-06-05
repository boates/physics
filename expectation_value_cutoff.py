#!/usr/bin/env python
"""
Returns a weighted mean.
"""
import sys

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        header = f.readline()
        lines = f.readlines()
        column1 = int(sys.argv[2])
        column2 = int(sys.argv[3])
        cutoff = float(sys.argv[4])
    except:
        print '\n usage: '+sys.argv[0]+' file_with_header.dat column1 column2 cutoff\n'
        sys.exit(0)

    # Perform the calculation
    mean = 0.0
    norm = 0.0
    for line in lines:
        value = float(line.split()[column1 - 1])
        weight = float(line.split()[column2 - 1])
        if weight > cutoff:
            mean += weight * value
            norm += weight

    # Output the result
    print 'Expectation value:', mean / norm


if __name__ == '__main__':
    main()
