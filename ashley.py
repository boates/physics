#!/usr/bin/env python
"""
Script for Ashley
"""

import os, sys, commands, glob

def main():

    x = int(input("How much do you love Brian Boates (out of 10):  "))

    if x > 10:
        print "\nI said, out of 10...\n"
    if 5 < x <= 10:
        print "\nCongratulations, you can give him a Wal-job\n"
    if x < 5:
        print "\n:(\n"

if __name__ == '__main__':
    main()
