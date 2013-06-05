#!/usr/bin/env python
"""
Report progress of CHAIN.cbn file being made (assuming made from TRAJEC.cbn file in same dir)
"""
import os, sys, commands, glob

def main():

    dirs = glob.glob('1.*')
    dirs.sort()
    for d in dirs:
        try:
            chain_lines = float(commands.getoutput("wc -l "+d+"/CHAIN.cbn | awk '{print $1}'"))
            trajec_lines = float(commands.getoutput("wc -l "+d+"/TRAJEC.cbn | awk '{print $1}'"))
            print d,':  ',round(chain_lines/trajec_lines*100,1),'%'
        except:
            print d,':  ','0.0 %'

if __name__ == '__main__':
    main()
