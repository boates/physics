#!/usr/bin/env python

# Script to fix the crappy fortran file writing format
# Customized towards stress.dat and statis.dat files

import commands, sys, os

def main():

    try:
        f1 = open('statis.dat','r')
        f2 = open('stress.dat','r')
    except:
        print '\nPlease make sure both statis.dat and stress.dat are present\n'
        sys.exit(0)

    sta_out = open('sta.out','w')
    str_out = open('str.out','w')

    N = int(commands.getoutput('wc -l statis.dat').split()[0]) / 2

    for i in range(N):
        sta_out.write(f1.readline().strip()+'  '+f1.readline().strip()+'\n')
        str_out.write(f2.readline().strip()+'  '+f2.readline().strip()+'\n')

    sta_out.close()
    str_out.close()

    os.system('rm stress.dat statis.dat')
    os.system('mv sta.out statis.dat')
    os.system('mv str.out stress.dat')

if __name__ == '__main__':
    main()
