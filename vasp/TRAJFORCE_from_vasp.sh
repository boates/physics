#!/bin/bash
# Extract TRAJFORCE data from vasp OUTCAR file and write to TRAJFORCE file for further processing

grep -A109 "POSITION                                       TOTAL-FORCE (eV/Angst)" OUTCAR | grep -v "\-\-" | grep -v POSITION > TRAJFORCE
