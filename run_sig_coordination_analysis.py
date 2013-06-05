#!/usr/bin/env python
"""
Run sig_coordination_analysis.py for all densities for all T
"""
import os, sys, commands, glob

def main():

    # Inform users and begin file search
    print '\n This code assumes you know what the hell you are doing, proceeding...\n'

    # Grab XXXXK directories and run get_conductivities.py
    temps = glob.glob('*K')
    cwd = os.getcwd()
    for temp in temps:
        os.chdir(temp)
        rss_rough = glob.glob('1.*')
        rss = [rs for rs in rss_rough if len(glob.glob(rs+'/sig*.dat')) != 0]
        rss.sort()
        T_dir = os.getcwd()
        for rs in rss:
            os.chdir(rs)
            os.system('sig_coordination_analysis.py '+rs)
            os.chdir(T_dir)
        os.chdir(cwd)

    # Make one huge file
#    os.system('tail -n-1 2000K/1.*/avg_cond_vs_coordination.dat solidified_2000K/1.*/avg_cond_vs_coordination.dat | grep "0\." > avg_cond_vs_coordination_2000K.dat')
    os.system('tail -n-1 2000K/1.*/avg_cond_vs_coordination.dat | grep "0\." > avg_cond_vs_coordination_2000K.dat')


if __name__ == '__main__':
    main()
