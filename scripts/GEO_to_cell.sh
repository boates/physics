#!/bin/bash

grep -A50 periodic $1 | awk '{print " ",$2,$3,$4}' | sed s/"  list Reduced coordinates"/"xred"/g | sed s/"  vectors of the"/"rprim"/g | grep -B50 data | grep -v data
