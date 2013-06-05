#!/bin/bash

# Run multiple continued DL_POLY runs for averaging purposes

run=1
increment=1
run_limit=10

while [ "$run" -le "$run_limit" ]
  do
    
  qsub -sync y dl.par
  mkdir run$run
  cp CONTROL CONFIG FIELD REV* run$run
  mv STATIS OUTPUT RDFDAT HISTORY dl_poly.* run$run
  $HOME/software/dl_poly/copy
  cp CONTROL.restart CONTROL
  run=`expr $run + $increment`

  done
