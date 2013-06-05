#!/usr/bin/env python
"""
Perform chain analysis
"""
import os, sys, commands, glob

def main():

    # Calculate chain_angles and chain_fractions using chain_angles.py, dissociation.x
    print 'Calculating chain_angles...'
    natom = int(commands.getoutput('head -1 TRAJEC.xyz').strip())
    os.system('mv -f TRAJEC.cbn FINAL-1000.cbn')
    os.system('tail -'+str((natom+2)*10000)+' wrapped.xyz > FINAL-10000.xyz')
    os.system('sed s/"FINAL-1000.xyz"/"FINAL-10000.xyz"/g cbn.in > cbn.in2')
    os.system('mv -f cbn.in2 cbn.in')
    os.system('nice -15 cbn_from_xyz.x < cbn.in > cbn.out')
    os.system('nice -15 chain_angles.py TRAJEC.cbn 1.50 90 90')
#    os.system('echo "CHAIN.cbn" > chain.in')
#    os.system('nice -15 chain_label.x < chain.in > chain_label.out')
    os.system('nice -15 dissociation_from_cbn.x < chain.in > dissociation.out')
    os.system('rm -f cbn.out chain.in chain_label.out dissociation.out')
    os.system('mv -f dissociation.dat chain_fraction.dat')
    os.system('running_avg.py chain_fraction.dat 2 50')
    os.system('mv -f avg.dat avg_chain_fraction.dat')
    print 'done'


if __name__ == '__main__':
    main()
