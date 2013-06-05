#!/usr/bin/env python
"""
Hold all current jobs in the queue
requires use of status command
"""
import os, sys, commands

def main():

    jids = commands.getoutput("status | grep -A100000 Waiting | tail -n+4 | awk '{print $1}'").split()
    for jid in jids:
        os.system('qhold '+jid)


if __name__ == '__main__':
    main()
