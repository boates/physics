#!/usr/bin/env python
"""
Based on a 0.8 fs timestep for Mg, give a scaled timestep
for another atomic mass
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        m_X = float(sys.argv[1])
    except:
        print '\n usage: '+sys.argv[0]+' mass_in_amu\n'
        sys.exit(0)

    m_Mg = 24.305
    t_Mg = 0.80

    t_X = t_Mg * ( m_X / m_Mg )**0.5

    print '\nFor a mass of',m_X,'amu, use a timestep of',round(t_X,3),'fs\n'


if __name__ == '__main__':
    main()
