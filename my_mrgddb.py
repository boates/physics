#!/usr/bin/env python

import os, sys, commands, glob

def main():

    # Detect present DDB files
    raw = glob.glob('abinit/part*/OUT_DS*_DDB')
    dirs = [d for d in raw if 'DS11' not in d]
    dirs.sort()
    N = len(dirs)

    # Create merged file header
    os.system('grep -B1000000 blocks= '+dirs[0]+' | sed s/"blocks=    1"/"blocks=    '+str(N)+'"/g > mrgddb.out')

    # Write the DDB info to the merged file
    for d in dirs:
        os.system('echo "" >> mrgddb.out')
        os.system('grep -A1 2rd '+d+' >> mrgddb.out')
        os.system("grep -A1000000 2rd "+d+" | grep -v 2rd | grep -v qpt | awk '{print $3,$4,$1,$2,$5,$6}' >> mrgddb.out")

    # Write other line
    os.system('echo "" >> mrgddb.out')
    os.system('echo " List of bloks and their characteristics" >> mrgddb.out')

    # Write end of merged DDB file
    for d in dirs:
        os.system('echo "" >> mrgddb.out')
        os.system('grep -A1 2rd '+d+' >> mrgddb.out')
        
     
if __name__ == '__main__':
    main()
