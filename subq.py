#!/usr/bin/env python
"""
qalter selected jobs to fill the sub.q on placentia2
"""
import os, sys, commands, glob

def main():

    # Number of free 16 & 8 cpu slots
    n16 = int(commands.getoutput('qstat -f -q sub.q | grep "0\/16" | grep -v S | wc -l'))
    n8  = int(commands.getoutput('qstat -f -q sub.q | grep "0\/8" | grep -v S | wc -l'))
    N = n16 + n8

    print 'n16 =', n16
    print 'n8  =', n8
    print 'n16 + n8 =', N

    if N == 0:
        print '\n No 16 or 8cpu nodes currently available in sub.q --- exiting...\n'
        sys.exit(0)

    j = int(input('How many jobs would you like to alter?\n '))

    print '\nWill first fill up 16cpu slots, followed by 8cpu slots until'
    print 'all given jobs are altered up to n16+n8 jobs\n'

    if j > N:
        print 'Requested number of jobs to alter > n16+n8 =',N
        print 'Will only do the first',N,' jobs given\n'

    # Retrieve PIDs from user
    yn = input('Are your PIDs stored in a file? (y=0/n=1)\n')
    if yn == 0:
        f = str( input('File with PIDs:\n') )
        PIDs = commands.getoutput("awk '{print $1}' "+f).split()[0:j]
    elif yn == 1:
        PIDs = []
        while len(PIDs) < j:
            PID = input('PID: ')
            PIDs.append(PID)
    else:
        print '\n Improper response --- proceeding assuming no file is present \n'

    # Fill the 16cpu slots
    for i in range(min(n16,j)):
        PID = str(PIDs.pop(0))
        current = commands.getoutput("qstat -j "+PID+" | grep 'hard resource_list' | awk '{print $3}'")
        new = current.replace('suspendable=false','susp=true')
        k = os.system('qalter -l '+new+' '+PID)
        if k == 0:
            print 'job:', PID, 'successfully altered'
        else:
            print 'PROBLEM altering job:', PID
        os.system('sleep 2s')

    # Fill the 8cpu slots
    for i in range(j - n16):
        PID = str(PIDs.pop(0))
        current = commands.getoutput("qstat -j "+PID+" | grep 'hard resource_list' | awk '{print $3}'")
        new = current.replace('suspendable=false','susp=true')
        k = os.system('qalter -l '+new+' '+PID)
        l = os.system('qalter -pe openmp 8 '+PID)
        if k+l == 0:
            print 'job:', PID, 'successfully altered'
        else:
            print 'PROBLEM altering job:', PID
        os.system('sleep 2s')


if __name__ == '__main__':
    main()
