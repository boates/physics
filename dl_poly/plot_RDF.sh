#!/bin/bash

# Used by get_RDF.py to plot RDF results

cmd="gnuplot.scr"
cat > $cmd << END
set data style linespoints
set xlabel "distance"
set ylabel "g(r)"
load "plot_RDF.tmp"
pause -1 "Showing g(r).\nPress enter to quit"
quit
END

gnuplot -background white -xrm 'gnuplot*line1Color:blue' $cmd

\rm $cmd

#rm plot_RDF.tmp