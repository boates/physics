#!/usr/bin/python
"""
force_match_vasp.py
Version: 2.0
Author: Brian Boates

Perform force-matching analysis on VASP output
to obtain pair-potentials for all species.

Assumptions:
 - Orthorhombic cell
"""
import os, sys, commands, glob
import optparse

def init_parser():

    usage   = "usage: %prog [options]"
    version = "%prog 1.0"
    parser  = optparse.OptionParser(usage=usage,version=version)

    parser.add_option('-f', '--file', type='string', dest='f', default='',
                      help='OUTCAR file to be used for force-matching \
                            [no default]')
    parser.add_option('-s', '--skip', type='int', dest='Nskip', default=0,
                      help='Number of steps to skip before force-matching \
                            [default=0]')
    parser.add_option('-b', '--block', type='int', dest='Nblock', default=100,
                      help='Number of steps in each block for force-matching \
                            [default=100]')
    parser.add_option('-z', '--stride', type='int', dest='stride', default=10,
                      help='Number of steps to skip within block for force-matching \
                            [default=10]')

    parser.add_option('-t', '--table', dest='table', default=True,
                      help='Flag to throw for table writing \
                            [default=True]')
    parser.add_option('-l', '--lambdas', dest='lambdas', default='0.00 0.01 0.02 0.04 0.06 \
                      0.08 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00',
                      help='List of lambdas to use to create scaled TABLE files \
                            Space separated string in qutations (see default) \
                            [default="0.00 0.01 0.02 0.04 0.06 0.08 0.10 0.15 0.20 \
                                      0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00"]')

    parser.add_option('-d', '--dlpoly', dest='dlpoly', default=True,
                      help='Flag to throw for creating DL_POLY input files: \
                            FIELD, CONFIG, & CONTROL \
                            [default=True]')
    parser.add_option('--densvar', dest='densvar', default='50',
                      help='CONTROL file parameter to allow for flexibility \
                            with DL_POLY neighbour lists, often needed > 0 \
                            [default=50]')
    parser.add_option('-v', '--levcfg', type='string', dest='levcfg', default='0',
                      help='DL_POLY levcfg parameter for CONFIG file \
                            (see DL_POLY documentation for all options) \
                            [default=0, i.e. only coordinates in CONFIG]')
    parser.add_option('-i', '--imcon', type='string', dest='imcon', default='1',
                      help='DL_POLY imcon parameter for CONFIG file \
                            (see DL_POLY documentation for all options) \
                            [default=1, i.e. cubic box]')
    parser.add_option('-k', '--keytrj', type='string', dest='keytrj', default='2',
                      help='DL_POLY keytrj for CONTROL file "traj" parameter\
                            (see DL_POLY documentation for all options) \
                            [default=2, i.e. write positions, velocities, & forces]')

    parser.add_option('--DLtemp', dest='DLtemp', default=False,
                      help='For CONTROL: Temperature for DL_POLY simulations \
                            [default=False, i.e. gets from OUTCAR file]')
    parser.add_option('--DLsteps', type='string', dest='DLsteps', default='20000',
                      help='For CONTROL: Number of steps for DL_POLY simulations \
                            [default=20000, i.e. ~10-20 picosecond]')
    parser.add_option('--DLtimestep', type='string', dest='DLtimestep', default='0.75',
                      help='For CONTROL: Timestep for DL_POLY simulations (in fs)\
                            [default=0.75]')
    parser.add_option('--DLstride', type='string', dest='DLstride', default='1',
                      help='CONTROL: Stride for writing to HISTORY & STATIS files \
                            [default=1, i.e. write to HISTORY at every step]')

    parser.add_option('-q', '--submit', dest='submit', default=False,
                      help='Use flag to submit DL_POLY calculations directly using \
                            this script, note: the executable path must be given \
                            (use the -e flag to pass the DL_POLY executable) \
                            [default=False]')
    parser.add_option('-e', '--exe', type='string', dest='exe', default='None',
                      help='path for DL_POLY exectuable for submission scripts \
                            (see DL_POLY documentation for all options) \
                            [default=None, must be provided]')

    parser.add_option('-n', '--nHISTORY', type='int', dest='nHISTORY', default=3000,
                      help='Number of steps to use at end of HISTORY file for \
                            lambda=1.00 calculations of the lambda!=1.00 ensembles \
                            (need enough for good thermodynamic averages) \
                            Note: stride can be controlled using --DLstride flag \
                            [default=3000]')
    parser.add_option('--DLclean', dest='DLclean', default=True,
                      help='Whether or not to delete DL_POLY output after \
                            obtaining the energies to preserve disk space \
                            [default=True]')
    parser.add_option('--checkTABLE', dest='checkTABLE', default=True,
                      help='Whether or not to check each TABLE file for comparison \
                            with original force-matched / numerical potentials \
                            [default=True]')

    return parser


def submit_script(exe):
    """
    Provide the submit script info for different machines
    based on the hostname

    Note: this function will have to be edited to allow for new machines
    
    Returns s and sub
    s: String to be used at the top of a submission script
    sub: qsub or msub command for given hostname
    """
    hostname = commands.getoutput('hostname')

    s = ''
    acenet = ['clhead', 'glooscap', 'fundy', 'mahone', 'placentia', 'brasdor']
    bonev = ['nova']
    llnl = ['hera', 'aztec', 'atlas', 'zeus']
    
    if [host for host in acenet if host in hostname]:
        s += '#!/bin/sh\n'
        s += '#$ -S /bin/sh\n'
        s += '#$ -cwd\n'
        s += '#$ -N dlpoly\n'
        s += '#$ -l h_vmem=1500M\n'
        s += '#$ -l h_rt=47:00:00\n\n'
        sub = 'qsub'
        if exe == 'None':
            exe = '/usr/local/dl_poly/bin/DLPOLY.Y'

    if [host for host in bonev if host in hostname]:
        s += '#$ -S /bin/bash\n'
        s += '#$ -cwd\n'
        s += '#$ -N dlpoly\n'
        s += '#$ -l h_vmem=2000M\n'
        s += '#$ -l h_rt=50:00:00\n\n'
        sub = 'qsub'
        if exe == 'None':
            exe = '/home/boates/bin/DLPOLY.Y'

    elif [host for host in llnl if host in hostname]:
        s += '#!/bin/sh\n'
        s += '#MSUB -A micphys\n'
        s += '#MSUB -l walltime=0:12:00:00\n'
        s += '#MSUB -l nodes=1\n'
        s += '#MSUB -l feature=aztec\n\n'
        sub = 'msub'

    else:
        print ' | \n | You have provided a valid executable - HOWEVER, the current'
        print ' | hostname is not recognized by this script and the jobs can not'
        print ' | be submitted, the currently acceptable hosts are LLNL machines:'
        print ' | hera, atlas, aztec, & zeus, and ACEnet machines: glooscap, fundy,'
        print ' | mahone, brasdor, & placentia'
        print ' | \n | ---> Thus, the DL_POLY jobs will not be submitted'
        print ' |      and submit scripts will not be created\n |'
        s = False
        sub = ''

    return s, sub, exe


def derivative(x,y):
    """
    Uses central difference method for all points less the
    first and final points where the forward and backward
    difference methods are used.
    """
    dydx = []
    for i in range(len(y)):
        if i == 0:
            dy = y[i+1] - y[i]
            dx = x[i+1] - x[i]
        elif i == len(y)-1:
            dy = y[i] - y[i-1]
            dx = x[i] - x[i-1]
        else:
            dy = y[i+1] - y[i-1]
            dx = x[i+1] - x[i-1]
        dydx.append( dy/dx )
    return dydx


def factorial(n):
    """
    Return the factorial of an integer n
    """
    f = 1
    for i in range(1,n+1):
        f *= i
    return f


def choose(m,n):
    """
    Return 'm choose n'
    """
    return factorial(m) / (factorial(n)*factorial(m-n))
    

def create_TABLE_file(s1, s2, rcore, rcut):
    """
    s1: species type 1
    s2: species type 2
    rcore: core radius in Bohr for INPUT_GROM file
    rcut: cut-off in angstroms for INPUT_GROM file
    """
    # INPUT_GROM file is the input file for f_table_writer.e
    out = open('INPUT_GROM','w')
    out.write('&begin\n')  # first line in INPUT_GROM file
    out.write('  %sites\n')  # specify the atomic species
    out.write('    '+s1+' '+s2+'\n')  # species types
    out.write('&\n')
    out.write('  %core\n')  # define core radius
    out.write('    '+s1+' '+s2+' '+rcore+'\n')  # core radius in Bohr
    out.write('&\n')
    out.write('  %rcut\n')  # define cut-off
    out.write('    rcut = '+rcut+'\n')
    out.write('&\n')
    out.write('  fit: table\n')  # fit for TABLE file
    out.write('&end\n')  # final line in INPUT_GROM
    out.close()

    # Do the TABLE writing
    os.system('f_table_writer.e >& tw.log')


def masses():
    """
    Creates an atomic mass dictionary
    """
    mass = {'H' : 1.00794,
            'He': 4.002602,
            'Li': 6.941,
            'Be': 9.012182,
            'B' : 10.811,
            'C' : 12.0107,
            'N' : 14.0067,
            'O' : 15.9994,
            'F' : 18.9984032,
            'Ne': 20.1797,
            'Na': 22.98976928,
            'Mg': 24.3050,
            'Al': 26.9815386,
            'Si': 28.0855,
            'P' : 30.973762,
            'S' : 32.065,
            'Cl': 35.453,
            'Ar': 39.948,
            'K' : 39.0983,
            'Ca': 40.078,
            'Sc': 44.955912,
            'Ti': 47.867,
            'V' : 50.9415,
            'Cr': 51.9961,
            'Mn': 54.938045,
            'Fe': 55.845,
            'Co': 58.933195,
            'Ni': 58.6934,
            'Cu': 63.546,
            'Zn': 65.38,
            'Ga': 69.723,
            'Ge': 72.64,
            'As': 74.92160,
            'Se': 78.96,
            'Br': 79.904,
            'Kr': 83.798,
            'Rb': 85.4678,
            'Sr': 87.62,
            'Y' : 88.90585,
            'Zr': 91.224,
            'Nb': 92.90638,
            'Mo': 95.96,
            'Tc': 98.0,
            'Ru': 101.07,
            'Rh': 102.90550,
            'Pd': 106.42,
            'Ag': 107.8682,
            'Cd': 112.411,
            'In': 114.818,
            'Sn': 118.710,
            'Sb': 121.760,
            'Te': 127.60,
            'I' : 126.90447,
            'Xe': 131.293,
            'Cs': 132.9054519,
            'Ba': 137.327,
            'La': 6.145,
            'Ce': 6.77,
            'Pr': 6.773,
            'Nd': 7.007,
            'Pm': 7.26,
            'Sm': 7.52,
            'Eu': 5.243,
            'Gd': 7.895,
            'Tb': 8.229,
            'Dy': 8.55,
            'Ho': 8.795,
            'Er': 9.066,
            'Tm': 9.321,
            'Yb': 6.965,
            'Lu': 174.9668,
            'Hf': 178.49,
            'Ta': 180.94788,
            'W' : 183.84,
            'Re': 186.207,
            'Os': 190.23,
            'Ir': 192.217,
            'Pt': 195.084,
            'Au': 196.966569,
            'Hg': 200.59,
            'Tl': 204.3833,
            'Pb': 207.2,
            'Bi': 208.98040,
            'Po': 210.0,
            'At': 210.0,
            'Rn': 222.0,
            'Fr': 223.0,
            'Ra': 226.0}

    return mass


def main():

    # Remove any files from a previous analysis
    os.system('rm -rf FORCE_OUT/ position_force TRAJECTORYn1 Input.in fmv.log w1.00')

    # Initiate option parser
    parser     = init_parser()
    options, args = parser.parse_args()
    f          = options.f
    Nskip      = options.Nskip
    Nblock     = options.Nblock
    stride     = options.stride
    table      = options.table
    lambdas    = options.lambdas.split()
    dlpoly     = options.dlpoly
    levcfg     = options.levcfg
    imcon      = options.imcon
    keytrj     = options.keytrj
    DLtemp     = options.DLtemp
    DLsteps    = options.DLsteps
    DLtimestep = str(float(options.DLtimestep)/1000.)
    DLstride   = options.DLstride
    submit     = options.submit
    exe        = options.exe
    nHISTORY   = options.nHISTORY
    DLclean    = options.DLclean
    checkTABLE = options.checkTABLE

    lambdas.sort()
    print
    # Retrieve OUTCAR filename if necessary
    if f == '':
        print " | You have not provided an OUTCAR file name with the -f flag"
        f = raw_input(" | Please provide OUTCAR filename:  ")
        print " |"

    # Exit if user's DLsteps, DLstride, & nHISTORY are inconsistent
    if float(DLsteps) / float(DLstride) < float(nHISTORY):
        print " | Values given for DLsteps, DLstride, & nHISTORY are inconsistent"
        print " | nHISTORY must be greater than DLsteps/DLstride, otherwise there won't"
        print " | be enough configurations in HISTORY for looping, currently you set"
        print " | DLsteps = "+DLsteps+", DLstride = "+DLstride+", & nHISTORY = "+str(nHISTORY)
        print " | \n | exiting...\n"
        sys.exit(1)

    #############################
    ### CREATING Input.in FOR ###
    ### FORCE-MATCHING CODE   ###
    print " | Creating Input.in for force-matching code"
    # Retrieve species list
    natoms = commands.getoutput("head -1000 "+f+" | grep 'ions per type'").split('=')[-1].split()
    Nspecies = len(natoms)
    typats = commands.getoutput("head -100 "+f+" | grep POTCAR | awk '{print $3}' | head -"+str(Nspecies)+"").split()
    species = []
    for i in range(Nspecies):
        species.append( [typats[i].split('_')[0], int(natoms[i])] )

    # Get lattice constants from OUTCAR [Ang]
    a, b, c = commands.getoutput("grep -A1 'length of vectors' "+f+" | tail -2 | grep -v length | awk '{print $1,$2,$3}'").split()
    a, b, c = float(a)/0.5291772, float(b)/0.5291772, float(c)/0.5291772

    # Get natom from OUTCAR
    natom = 0
    for i in range(Nspecies): natom += species[i][1]

    # Determine half the shortest lattice vector [Bohr]
    half = min([a, b, c]) / 2.0

    # Create the position_force file in [Bohr], [Ha/Bohr]
    os.system("grep -A"+str(natom+1)+" 'TOTAL-FORCE' "+f+" | grep -v TOTAL | grep -v '\-\-' | awk '{print $1/0.5291772,$2/0.5291772,$3/0.5291772,$4*0.5291772/27.211383,$5*0.5291772/27.211383,$6*0.5291772/27.211383}' > position_force")

    # Determine the number of steps
    Nsteps = int(commands.getoutput("grep pressure "+f+" | wc -l"))

    # Create the TRAJECTORYn1 file & remove position_force
    os.system("awk '{if (n % "+str(natom)+" == 0) {print n/"+str(natom)+" + 1;print $0} else print $0;n+=1}' position_force > TRAJECTORYn1")
    os.system("echo "+str(Nsteps + 1)+" >> TRAJECTORYn1")
    os.system('rm -f position_force')

    # Write the Input.in file
    out = open('Input.in','w')
    out.write('&begin Input\n\n')   # first line in Input.in
    out.write('number Types '+str(Nspecies)+'\n')   # number of different atomic species
    for i in range(Nspecies):                                         # Loop:
        out.write(" '"+species[i][0]+"'  "+str(species[i][1])+"\n")   #   species  natom
    out.write('\nCELL \ \n')                           # cell parameters flag
    out.write(str(a)+' '+str(b)+' '+str(c)+'\n\n')     # a, b, & c in Bohr
    out.write('path to potential \ \n')         # path for TRAJECTORYn1 file
    out.write('"./"\n')                         # current directory
    out.write('PERIODIC\n')   # use PBC
    out.write('THREE DIM\n')  # fit all three projections of reference forces
    out.write('ROTATE\n')     # apparently good option to use
    out.write('TEST FORCES NO\n\n')   # apparently good to use
    out.write('write only\n')    # write [all, only, ea]
    out.write('read H_xyz\n\n')  # has always been here
    out.write('Solver pda 1.0d-4 1.0d-4\n\n')  # subroutine to solve overdetermined eq'n and tolerance
    out.write('AVERAGING OVER TRAJECT\n\n')  # always use
    out.write('SKIP \ \n')        # number of steps to skip before force-match
    out.write(str(Nskip)+'\n\n')  # can be user-defined
    out.write('AVERAGE  IN\ \n')    # This number * steps in block = total number of configurations
    out.write(str(Nsteps)+'\n\n')   # total number of configs works as upper bound, code exits at eof anyway
    out.write('Steps in block \ \n')  # number of steps to use to build each overdetermined system
    out.write(str(Nblock)+'\n\n')     # can be user-defined
    out.write('Stride \ \n')       # number of steps to skip within trajecory block used to build
    out.write(str(stride)+'\n\n')  # overdetermined system, can be user-defined
    out.write('Radial mesh\n')  # Radial mesh for short ranged interaction
    for i in range(Nspecies):                # loop over all species pairs for interactions
        for j in range(i,Nspecies):          # length of interaction defined by alat / 2.
            out.write(species[i][0]+'-'+species[j][0]+'    0-'+str(half)+':0.05 \n')
    out.write('end\n\n')  # end Radial mesh
    out.write('&end Input \n')  # final line in Input.in
    out.close()
    print " |\n | done"

    ##########################
    ### RUN FORCE-MATCHING ###
    ###                    ###
    print ' |\n | Running force-matching code...\n |'
    os.system('m.e >& fmv.log')
    os.system('rm -f TRAJECTORYn1')
    print ' | done\n |'
    ##########################

#    print ' | Adding linear extrapolations to r=0 and V=0 for large r'
#    print ' | and replacing old Poten and Force files with them...\n |'
    ##################################
    ### MODIFY OBTAINED POTENTIALS ###
    ###                            ###
    # Extrapolate linearly to r=0 and r=r_max and create eV & eV/Angst files
#    cwd = os.getcwd()
#    os.chdir('FORCE_OUT')
#    os.mkdir('FINAL')

    # Loop over Poten_XX.dat files
#    potens = glob.glob('Poten_*.dat')
#    for p in potens:

        # Read in the Poten_XX.dat file and trim the beginning
#        pp = open(p,'r')
#        lines = pp.readlines()
#        pp.close()
        # Remove the "flat" region of pp for small r
#        V, r  = [], []
#        v0 = lines[0].split()[1]  # get plateau value
#        for line in lines:
#            ri, vi = line.split()
#            if vi != v0:
#                V.append(float(vi))
#                r.append(float(ri))
#        # Remove first few values, where potential begins plateauing
#        V = V[4:]
#        r = r[4:]

        # Linear extrapolation to r=0
#        dr = r[1] - r[0]
#        c1 = V[0]
#        c2 = (V[1] - V[0]) / dr
#        u = c1 - c2*r[0]
#        v = -c2
#        r_tmp = range(0,int( (r[0]-dr)*1000000. ),int(dr*1000000.))
#        r_head, V_head = [], []
#        for i in range(len(r_tmp)):
#            r_float = float(r_tmp[i]) / 1000000.
#            r_head.append( r_float )
#            V_head.append( u - v*r_float )
# numpy  r_head = numpy.arange(0.0,r[0]-dr,dr)
# numpy  V_head = a - b*r_head

#        r_tail, V_tail = [], []
        # Linear extarpolation at end to V=0
#        if abs(V[-1]) > 0.0001:
#            dr = r[-1] - r[-2]
#            c1 = V[-2]
#            c2 = (V[-1] - V[-2]) / dr
#            u = c1 - c2*r[-1]
#            v = -c2
#            r_tmp = range(int(r[-1]*1000000.),100*1000000,int(dr*1000000.))
#            r_tail, V_tail = [], []
#            for i in range(len(r_tmp)):
#                r_float = flaot(r_tmp[i]) / 1000000.
#                r_tail.append( r_float )
#                V_tail.append( u - v*r_float )
## numpy     r_tail = numpy.arange(r[-1],100,dr)
## numpy     V_tail = a - b*r_tail
#            for i in range(len(r_tail)):
#                print r_tail[i], V_tail[i]
#            # Find where V_tail crosses zero and end there
#            for i in range(len(V_tail)-1):
#                VV = V_tail[i]*V_tail[i+1]
#                if VV < 0:
#                    end = i+1
#            r_tail = r_tail[:end+1]
#            V_tail = V_tail[:end]+[0.0]
#        else:
#            r_tail = [r[-1]+dr]
#            V_tail = [0.0]
            
        # Concatenate all the pieces of r and V
#        r = list(r_head) + r + r_tail
#        V = list(V_head) + V + V_tail

        # Take derivatives for forces
#        F = derivative(r,V)

#        Vout1 = open('FINAL/'+p,'w')
#        Vout2 = open('FINAL/'+p.replace('dat','evang'),'w')
#        Fout1 = open('FINAL/'+p.replace('Poten','Force'),'w')
#        Fout2 = open('FINAL/'+p.replace('dat','evang').replace('Poten','Force'),'w')
#        for i in range(len(r)):
#            Vout1.write(str(r[i])+' '+str(V[i])+'\n')  # in Bohr and Hartree
#            Vout2.write(str(r[i]*0.5291772)+' '+str(V[i]*27.211383)+'\n')  # in Ang and eV
#            Fout1.write(str(r[i])+' '+str(-F[i])+'\n')  # in Bohr and Ha/Bohr
#            Fout2.write(str(r[i])+' '+str(-F[i]*27.211383/0.5291772)+'\n')  # in Ang and eV/Ang
#        Vout1.close()
#        Vout2.close()
#        Fout1.close()
#        Fout2.close()

    # Create Nod_Force_XX.dat files necessary for DL_POLY TABLE creation
#    forces = glob.glob('FINAL/Force_*.dat')
#    for force in forces:
#        os.system("awk '{print $1,$2,$2}' "+force+" > "+force.replace('Force','Nod_Force'))
#        os.system('cp '+force+' '+force.replace('Force','Nod_Force'))

    # Remove all original force-match files and
    # replace with the new ones just created
#    os.system('rm -r FORCE_OUT/RADIALS *.dat*.out')
#    os.system('mv FINAL/* ./')
#    os.system('rm -r FINAL')
#    os.chdir(cwd)
#    print ' | done'

    ##################################
    ###                            ###
    ###  FORCE MATHCED POTENTIALS  ###
    ###  ARE FINISHED              ###
    ###                            ###
    ###  NOW DO THE TABLE WRITING  ###
    ###  IF DESIRED                ###
    ###                            ###
    ##################################
    if table == True:

        print ' |\n | Doing the table writing for all lambdas...'

        # Create a template directory where lambda=1.00 (remove old first)
        os.system('rm -rf template_1.00')
        os.mkdir('template_1.00')
        # Move the FORCE_OUT directory to the template
        # directory, for convenience
        os.system('mv -f FORCE_OUT/ template_1.00/')

        # Write lambdas to a file for post-processing
        os.system('rm -f lambdas')
        L_out = open('lambdas','w')
        for L in lambdas:
            L_out.write(L+'\n')
        L_out.close()

        # Loop over lambdas and species pairs
        print ' | L =',
        for L in lambdas:
            os.system('rm -rf '+L)
            os.system('cp -r template_1.00 '+L)
            cwd = os.getcwd()
            os.chdir(L)
            print L,
            sys.stdout.flush()
            for i in range(Nspecies):
                for j in range(i,Nspecies):

                    # Define shorthand variables for species pair
                    s1, s2 = species[i][0], species[j][0]
                    # Annoying uppercase problems with force-matching output
                    s1u, s2u = species[i][0].upper(), species[j][0].upper()

                    # Get rcore and rcut for f_table_writer.e (must be in ANGSTROMS)
#                    rcore = '0.00'  # 0.00 because we linearly exrapolated the potentials to r=0
                    ff = open('FORCE_OUT/Force_'+s1u+s2u+'.dat','r')
                    lines = ff.readlines()
                    ff.close()
                    while float(lines[0].split()[1]) == 0.0: lines.pop(0)
                    lines.pop(0)
                    rcore = str(float(lines[0].split()[0]))         # in Bohr
                    rcut  = float(lines[-1].split()[0])*0.5291772   # in Ang
                    rcut  = str( round(int(rcut*100)/100.,2) - 0.01)[0:4]

#                    rcore = commands.getoutput("grep -v \"0\.0000\" FORCE_OUT/Force_"+s1u+s2u+".dat | grep -v grep | head -2 | tail -n-1 \
#                                                | awk '{print $1*1.}'")  # 2nd nonzero value in Force_XX.dat file in Bohr
#                    rcut = float(commands.getoutput("tail -n-1 FORCE_OUT/Poten_"+s1u+s2u+".dat \
#                                                     | awk '{print $1*0.5291772}'"))  # largest distance in Ang
#                    rcut = str( round(int(rcut*100)/100.,2) - 0.01)[0:4]

                    # Replace the Nod_Force_XX.dat files with lambda scaled ones
                    os.system("awk '{print $1,$2*"+L+",$3*"+L+"}' FORCE_OUT/Nod_Force_"+s1u+s2u+".dat > tmp")
                    os.system("mv tmp FORCE_OUT/Nod_Force_"+s1u+s2u+".dat")
                    os.system('mv FORCE_OUT/Nod_Force_'+s1u+s2u+'.dat FORCE_OUT/Nod_Force_XX.dat')

                    # Create the TABLE file
                    create_TABLE_file('X', 'X', rcore, rcut)

                    # Rename the Nod_Force file back to its proper name
                    os.system('mv FORCE_OUT/Nod_Force_XX.dat FORCE_OUT/Nod_Force_'+s1+s2+'.dat')

                    # Fix TABLE spacing and header
                    os.system('sed s/"E\-"/"XXX"/g OUT/TABLE | sed s/"\-"/" -"/g | sed s/"XXX"/"E-"/g > tab')
                    os.system('grep -v X tab > bottom')
                    tmp = open('OUT/TABLE','r')
                    header = tmp.readline().split('end')[-1].replace('|','')
                    tmp.close()
                    rvdw = header.split()[-2]
                    header += '       '+s1+'       '+s2+'\n'
                    head = open('header','w')
                    head.write('lambda='+L+'\n')
                    head.write(header)
                    head.close()
                    os.system('cat header bottom > OUT/TABLE')
                    os.system('rm -f tab header bottom')

                    # Create lambda scaled Poten files in eV & Ang for comparison with TABLE files
                    os.system("awk '{print $1*0.5291772,$2*27.211383*"+L+"}' FORCE_OUT/Poten_"+s1u+s2u+".dat \
                                                                           > FORCE_OUT/Poten_"+s1u+s2u+"_L"+L+".evang")
                
                    # Rename the OUT/ directory
#                    os.system('mv FORCE_OUT/Force_'+s1u+s2u+'.dat FORCE_OUT/Force_'+s1+s2+'.dat')
#                    os.system('mv FORCE_OUT/Poten_'+s1u+s2u+'.dat FORCE_OUT/Poten_'+s1+s2+'.dat')
#                    os.system('mv FORCE_OUT/Force_'+s1u+s2u+'.evang FORCE_OUT/Force_'+s1+s2+'.evang')
#                    os.system('mv FORCE_OUT/Poten_'+s1u+s2u+'.evang FORCE_OUT/Poten_'+s1+s2+'.evang')
                    os.system('mv -f OUT/force_XX.dat OUT/force_'+s1+s2+'.dat')
                    os.system('mv -f OUT/ener_XX.dat OUT/ener_'+s1+s2+'.dat')
                    os.system('rm -f OUT/FIELD_TAB.dat')
                    os.system('mv INPUT_GROM tw.log OUT/')
                    os.system('mv -f OUT/ TABLE_'+s1+s2)

            # Concatenate all TABLE files for use in DL_POLY, prepare FIELD.tab
            tables = glob.glob('*/TABLE')
            os.system('cp '+tables[0]+' '+tables[0]+'.tmp')        
            for tab in tables[1:]:
                tmp = commands.getoutput("head -2 "+tab+" | tail -n-1")
                os.system('grep -v lambda '+tab+' | grep -v "'+tmp+'" > '+tab+'.tmp')
            os.system('cat */TABLE.tmp > TABLE')
            os.system('rm -f */TABLE.tmp')
            os.system("grep -v E TABLE | grep -v lambda | awk '{print $1,\"   \",$2,\"    tab\"}' > FIELD.tab")
            os.system('mv FIELD.tab ../')

            os.chdir(cwd)

        # Remove the template directory
        os.system('rm -rf template_1.00')
        
        print '\n |\n | done'

        #####################################
        ###                               ###
        ###  TABLE FILES ARE ALL CREATED  ###
        ###                               ###
        ###  NOW CREATE DL_POLY INPUT     ###
        ###                               ###
        #####################################

        mass = masses()

        # Create DL_POLY input files if desired (FIELD, CONFIG, & CONTROL)
        # Same for all lambda!
        if dlpoly == True:

            print ' |\n | Creating FIELD, CONFIG, CONTROL, & clean files and w1.00/ dirs for DL_POLY'
            
            #########################
            ### Create FIELD file ###
            #########################
            os.system('rm -f FIELD')  # remove previous
            out = open('FIELD','w')
            out.write('FIELD file created by '+sys.argv[0]+'\n')  # header
            out.write('units eV\n')  # set DL_OPLY energy units to eV/mol
            out.write('molecules '+str(Nspecies)+'\n')  # number of atomic species
            # Loop over each species
            for i in range(Nspecies):
                out.write(species[i][0]+' atoms\n')  # atom type specific comment
                out.write('nummols '+str(species[i][1])+'\n')  # number of atoms of species type
                out.write('atoms 1\n')  # 1 atom in this "molecule" (so it is just an atom)
                out.write(species[i][0]+'        '+str(mass[species[i][0]])+'    0.000\n')  # species name, mass, & charge
                out.write('finish\n')
            out.write('vdw '+str(choose(Nspecies,2)+Nspecies)+'\n')  # number of pairs is 'N choose 2'
            tabs = open('FIELD.tab','r')  # retrieve atomic pairs from recently created FIELD.tab file
            out.write(tabs.read())  # write this info to the FIELD file
            out.write('close\n')  # final line in FIELD
            out.close()
            os.system('rm -f FIELD.tab')

            
            ##########################
            ### Create CONFIG file ###
            ##########################
            os.system('rm -f CONFIG')  # remove previous
            out = open('CONFIG','w')
            out.write('CONFIG of final configuration in '+f+'\n')  # header
            k = '        '
            out.write(k + levcfg + k + '  ' + imcon +'\n')
            # Convert a, b, & c to angstroms
            a, b, c = a*0.5291772, b*0.5291772, c*0.5291772
            out.write(k+' % .8f' % float(a)+k+' % .8f' % float(0)+k+' % .8f' % float(0)+'\n')
            out.write(k+' % .8f' % float(0)+k+' % .8f' % float(b)+k+' % .8f' % float(0)+'\n')
            out.write(k+' % .8f' % float(0)+k+' % .8f' % float(0)+k+' % .8f' % float(c)+'\n')
            
            # Get final configuration from given OUTCAR file
            coords = commands.getoutput("tail -10000 "+f+" | grep -A"+str(natom+1)+" TOTAL-FORCE | \
                                         tail -"+str(natom)+" | awk '{print $1,$2,$3}'").split('\n')
            # Write configuration to CONFIG file with very specific formatting
            count = 0
            for i in range(Nspecies):
                for j in range(species[i][1]):
                    count += 1
                    out.write(species[i][0] + k + k + str(count) +'\n')  # write species and index
                    # Loop over x, y, & z components of atomic coordinate
                    for r in coords[count-1].split():
                        if float(r) >= 0:
                            out.write('       '+'% .8f' % float(r))
                        elif float(r) < 0:
                            out.write('      '+'% .8f' % float(r))
                    out.write('\n')
            out.close()


            ###########################
            ### Create CONTROL file ###
            ###########################
            os.system('rm -f CONTROL')  # remove previous
            out = open('CONTROL','w')
            out.write('generated by '+sys.argv[0]+'\n\n')
            out.write('restart scale\n\n')
            out.write('integrator velocity verlet\n')
            # Get temperature from given OUTCAR file
            if DLtemp:
                temp = DLtemp
            else:
                temp = int(float(commands.getoutput("head -1000 "+f+" | grep TEBEG").split()[2].replace(';','')))
            out.write('temp  '+str(temp)+'\n')
            out.write('ensemble nvt hoover 0.05\n')
            out.write('no elec\n\n')
            out.write('steps  '+DLsteps+'\n')
            out.write('timestep  '+DLtimestep+'\n\n')
            out.write('scale  200\n')
            out.write('stack  100\n')
            out.write('print    1\n')
            out.write('stats    '+DLstride+'\n')
            out.write('rdf     10\n\n')
            dl_rcut = rvdw
            if float(rvdw) > half*0.5291772:                          # if rvdw is bigger than halfbox
                dl_rcut = str(round(int(half*0.5291772*100)/100.,2))  # set to just under half box
            out.write('cutoff      '+dl_rcut[:6]+'\n')
            out.write('rvdw        '+rvdw[:6]+'\n')
            out.write('delr width  0.5000\n\n')
            out.write('print rdf\n')
            out.write('traj  0 '+DLstride+' '+keytrj+'\n\n')
            out.write('densvar  50\n\n')
            out.write('job time    999999.9\n')
            out.write('close time      99.9\n\n')
            out.write('finish\n')
            out.close()


            # Create the template directory w1.00/ for total energy calculations using
            # unscaled interactions (lambda=1.00) on the scaled ensembles
            # CONTROL #
            os.mkdir('w1.00')
            os.mkdir('w1.00/template')
            os.system('grep -v restart CONTROL | grep -v temp | grep -v steps | grep -v traj | \
                       grep -v finish > tmp_control')
            out = open('w1.00_control_append','w')
            out.write('\nsteps  0\n')  # 0 steps for a single snapshot calculation
            out.write('traj 0 1 0\n')  # no need to print velocities and forces
            out.write('\nfinish\n')  # place the finish flag back at end of CONTROL file
            out.close()
            os.system('cat tmp_control w1.00_control_append > w1.00/template/CONTROL')
            os.system('rm -f tmp_control w1.00_control_append')

            # FIELD #
            os.system('cp FIELD w1.00/template/')  # same FIELD file as simulations

            # TABLE #
            one = [d for d in lambdas if d[0]=='1'][0]  # get unscaled (lambda=1.00) directory
            os.system('cp '+one+'/TABLE w1.00/template/')

            # clean #
            out = open('clean','w')
            out.write('#!/bin/bash\n\n')
            out.write('rm -f OUTPUT HISTORY STATIS REVCON REVIVE RDFDAT\n')
            out.close()
            os.system('chmod +x clean')

            # Copy FIELD, CONFIG, CONTROL, and clean files and w1.00/ directory to lambda directories
            for L in lambdas:
                os.system('cp FIELD   '+L)
                os.system('cp CONFIG  '+L)
                os.system('cp CONTROL '+L)
                os.system('cp clean   '+L)
#                if float(L) != 1.0:
                os.system('cp -r w1.00/ '+L)
            os.system('rm -rf FIELD CONFIG CONTROL clean w1.00/')

            print ' |\n | done'


            ######################################
            ###                                ###
            ###  CHECK TABLE FILES IF DESIRED  ###
            ###                                ###
            ######################################

            if checkTABLE:

                print ' | \n | Creating DL_POLY input for user requested TABLE checking'

                ### CREATE INPUT FILES FOR TABLE CHECKING ###

                for L in lambdas:

                    cwd = os.getcwd()
                    os.chdir(L)
                    os.system('mkdir checkTABLE/')
                    os.system('cp CONTROL FIELD TABLE checkTABLE/')
                    os.chdir('checkTABLE/')

                    for i in range(Nspecies):           # loop over all species
                        for j in range(i,Nspecies):     # pairs for TABLE tests
                            pair = species[i][0]+species[j][0]
                            os.mkdir(pair)
                            os.chdir(pair)

                            ### CONTROL FILE ###
                            os.system('cp ../../w1.00/template/CONTROL ./')

                            ### FIELD FILE ###
                            os.system('cp ../../FIELD ../')   # Get a "fresh" FIELD file to alter
                            N_in_pair = []
                            for k in range(Nspecies):
                                p = pair.split(species[k][0])   # Find out how many
                                if len(p) == 3:                 # of atom type are in 
                                    N_in_pair = 2               # current atom pair
                                else:                           #
                                    N_in_pair = p.count('')     #
#                                N_in_pair = list(pair).count(species[k][0])
                                os.system('sed s/"nummols '+str(species[k][1])+'"/"nummols '+str(N_in_pair)+'"/g ../FIELD > FIELD')
                                os.system('cp FIELD ../FIELD')

                            ### TABLE FILE ###
                            os.system('cp ../TABLE ./')

                            ### CONFIG FILE ###
                            out = open('CONFIG','w')
                            out.write('CONFIG file\n')
                            out.write('        0          1\n')
                            out.write('          20.0000000          0.00000000          0.00000000\n')
                            out.write('          0.00000000          20.0000000          0.00000000\n')
                            out.write('          0.00000000          0.00000000          20.0000000\n')
                            out.write(species[i][0]+'                1\n')
                            out.write('     0.00  0.0000  0.0000\n')
                            out.write(species[j][0]+'                2\n')
                            out.write('     XXXX  0.0000  0.0000\n')
                            out.close()

                            ### MOVE INPUT TO A TEMPLATE DIRECTORY ###
                            os.mkdir('template/')
                            os.system('mv CONTROL CONFIG FIELD TABLE template/')
                            
                            ### EXECUTION FILE ###
                            s, sub, exe = submit_script(exe)
                            out = open('run.py','w')
                            out.write('#!/usr/bin/python\n\n')
                            out.write('import os\n\n')
                            out.write('os.system("touch E.dat")\n\n')
                            out.write('for i in range(1,'+str(int(float(rcut)*100.))+'):\n')
                            out.write('    os.system("cp -r template/ tmp/")\n')
                            out.write('    os.chdir("tmp/")\n')
                            out.write('    os.system(\'sed s/"XXXX"/"\'+str(i/100.0)+\'"/g CONFIG > CONFIG.tmp\')\n')
                            out.write('    os.system("mv CONFIG.tmp CONFIG")\n')
                            out.write('    os.system("'+exe+'")\n')
                            out.write('    os.system("head -4 STATIS | tail -n-1 | awk \'{print $4/1.}\' > E")\n')
                            out.write('    os.system("echo "+str(i/100.)+" > r")\n')
                            out.write('    os.system("paste r E >> ../E.dat")\n')
                            out.write('    os.chdir("../")\n')
                            out.write('    os.system("rm -r tmp/")\n')
                            out.close()
                            os.system('chmod +x run.py')

                            os.chdir('../')
                    os.system('rm -f CONTROL FIELD TABLE')

                    os.chdir(cwd)

            print ' | \n | done'
            

            ##################################
            ###                            ###
            ###  ALL DL_POLY CALCULATIONS  ###
            ###  ARE NOW READY TO SUBMIT   ###
            ###                            ###
            ##################################

            #######################################
            ### CREATE dlpoly.s SUBMISSION FILE ###
            #######################################
                
            # Retrieve machine specific submit script header
#            s = False
            # If executable given AND given executable exists
#            if exe != 'None' and commands.getoutput('ls '+exe) == exe:
            s, sub, exe = submit_script(exe)
                    
            # If the user wants to submit DL_POLY jobs directly from this script
            if submit != False and exe == 'None':
                print ' | \n | You have selected to submit DL_POLY jobs directly from this'
                print ' | script, but did NOT provide a DL_POLY executable with -e flag'
                print ' | \n | ---> Thus, the DL_POLY jobs will not be submitted'
                print ' |      and submit scripts will not be created\n |'

            # If executable given AND given executable exists AND acceptable hostname
            if s:
                print ' |\n | Creating machine specific submission script (dlpoly.s)'
                # Create submit script
                out = open('dlpoly.s','w')
                out.write(s)

                # First, check the TABLE files if requested by the user
                if checkTABLE:
                    for i in range(Nspecies):
                        for j in range(i,Nspecies):
                            pair = species[i][0]+species[j][0]
                            out.write('cd checkTABLE/'+pair+'\n')
                            out.write('./run.py\n')
                            out.write('cd ../../\n\n')
                
                # Execute a DL_POLY simulation for given lambda
                out.write(exe+' >& LOG\n\n')  # capture DL_POLY notes/errors onto LOG file

                # Now copy HISTORY to the w1.00/ for lambda's != 1.00 and
                # compute lambda=1.00 energies of the lambda!=1.00 ensembles
                out.write('cp HISTORY w1.00/\n')
                out.write('cd w1.00\n')
                # To determine which timesteps from HISTORY file will be used
                out.write("tsteps=`grep timestep HISTORY | awk '{print $2}' | tail -n-"+str(nHISTORY)+"`\n")
                out.write("natom=`head -2 HISTORY | tail -n-1 | awk '{print $3}'`\n")
                out.write("levcfg=`head -2 HISTORY | tail -n-1 | awk '{print $1}'`\n")
                out.write("for t in $tsteps; do\n")  # loop over timesteps
                out.write("    echo $t\n")
                out.write("    cp -r template $t\n")  # create dir for snapshot based on template
                out.write("    head -2 HISTORY > CONFIG\n")                 # create CONFIG file
                out.write('    a=`echo "($natom*($levcfg+2)+4)*($t+1)+2" | bc`\n')  # need these numbers to parse
                out.write('    b=`echo "$natom*($levcfg+2)+3" | bc`\n')             # HISTORY for CONFIG
                out.write("    head -$a HISTORY | tail -n-$b >> CONFIG\n")  # from HISTORY
                out.write("    mv -f CONFIG $t\n")
                out.write("    cd $t\n")
                out.write("    "+exe+" >& LOG\n")  # run DL_POLY on snapshot
                out.write("    head -4 STATIS | tail -n-1 | awk '{print $3*1.0}' >> ../E\n")  # get energy
                out.write("    cd ../\n")
                if DLclean == True:  # delete DL_POLY files to preserve disk space if requested
                    out.write("    rm -rf $t\n")
                out.write("    cat E | awk '{sum+=$1} END {print sum/NR}' > E.avg\n\n")  # calculate average TOTAL E
                out.write("done\n\n")
                # Get average energy
                out.write("cat E | awk '{sum+=$1} END {print sum/NR}' > E.avg\n\n")  # calculate final average TOTAL E
                out.write("rm -f HISTORY\n")  # remove extra HISTORY file
                out.close()
                os.system('chmod +x dlpoly.s')

                # Copy dlpoly.s to lambda directories (special case for lambda=1.00)
                for L in lambdas:
                    os.system('cp dlpoly.s '+L)
#                    if float(L) == 1.0:
#                        os.system('grep -B10000 HISTORY dlpoly.s > '+L+'/dlpoly.s')
                os.system('rm -f dlpoly.s')

                ################################
                ### SUBMIT JOBS IF REQUESTED ###
                ################################
                if submit != False:
                    cwd = os.getcwd()
                    print ' | \n | Submitting DL_POLY jobs for L =',
                    for L in lambdas:
                        os.chdir(L)
                        os.system(sub+' dlpoly.s >& jobid')
                        os.chdir(cwd)
                        print L,
                        sys.stdout.flush()
                    print
                print ' |\n | done'
    print



if __name__ == '__main__':
    main()

