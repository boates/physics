#!/usr/bin/env python
"""
Given a starting number works out the Collatz conjecture
 f(n) = { n/2   ,  n even
        { 3n+1  ,  n odd
"""
import sys

def collatz(n,steps=0):

    while n != 1:

        steps += 1

        if n % 2 == 0:
            n, steps = collatz( n/2, steps )
        elif n % 2 == 1:
            n, steps = collatz( 3*n+1, steps )

    return n, steps

def steps_checker(Nmin=1,Nmax=10000):
    """
    Calculate number of collatz steps for all integers
    from Nmin to Nmax and write to collatz_steps.dat
    """
    # Loop over desired range and keep track of steps
    nsteps = []
    interval = range(Nmin,Nmax)
    for n in interval:
        m, steps = collatz(n,steps=0)
        nsteps.append(steps)

    # Write to an output file
    out = open('collatz_steps.dat','w')
    for i in range(len(interval)):
        out.write(str(interval[i])+'    '+str(nsteps[i])+'\n')
    out.close()

def main():

    # Check for a single number test
    if len(sys.argv) > 1:
        try:
            m, steps = collatz(int(sys.argv[1]))
            if m == 1:
                print '\n',str(int(sys.argv[1])),'satisfies the Collatz conjecture after ' \
                                                                           +str(steps)+' steps\n'
        except IndexError:
            print '\nsys.argv[1] should be an integer\n'
            sys.exit(0)

    # Find the number of steps for desired range
    steps_checker()


if __name__ == '__main__':
    main()
