#!/usr/bin/env python
"""
XDATCAR_cat.py

Finds and sorts all ./XDATCAR.### files and removes
headers appropriately to concatenate them all together
"""
import os, sys, commands, glob

def main():

    # Locate all present XDATCAR.* files
    xdats = glob.glob('XDATCAR.*')
    if len(xdats) == 0:
        print '\n No XDATCAR.* files detected --- exiting...\n'
        sys.exit(0)

    # Sort XDATCAR filenames
    xdats.sort()

    # Begin the file with the first XDATCAR
    xdat0 = xdats.pop(0)
    os.system('cp '+xdat0+' XDATCAR.cat')

    # Remove the headers and cat the rest
    for xdat in xdats:
        lines = int(commands.getoutput('wc -l '+xdat).split()[0])
        os.system('tail -'+str(lines-5)+' '+xdat+' >> XDATCAR.cat')
    
        
if __name__ == '__main__':
    main()
