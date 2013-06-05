#!/usr/bin/env python
"""
Script to retrieve CORERADIUSINBOHR and CUTOFFINANGSTROMS for INPUT_GROM
files for table writing
"""
import os,sys,commands,glob

def main():

    try:
        fin = sys.argv[1]
    except IndexError:
        print '\nPlease give Poten_NN.dat filename\n'
        sys.exit(0)

    os.system('cp /home/boates/bin/INPUT_GROM.template ./')

    rcut = float(commands.getoutput("tail -1 "+fin+" | awk '{print $1}'"))*0.529177
    rcut = str( round(int(rcut*100 - 1)/100.,2) )
    poten = commands.getoutput("awk '{print $2}' "+fin).split()
    for i in range(poten.count(poten[0])):
        poten.pop(0)
    core_poten = poten[4]
    core = str(round(float(commands.getoutput("grep "+core_poten+" "+fin+" | awk '{print $1}'")),2))
    os.system('sed -e s/CORERADIUSINBOHR/'+core[0:4]+'/g INPUT_GROM.template > INPUT_GROM.tmp')
    os.system('sed -e s/CUTOFFINANGSTROMS/'+rcut[0:4]+'/g INPUT_GROM.tmp > INPUT_GROM')
    os.system('rm -f INPUT_GROM.tmp INPUT_GROM.template')

    print '\nrcore =',core
    print 'rcut =',rcut,'\n'

if __name__ == '__main__':
    main()
