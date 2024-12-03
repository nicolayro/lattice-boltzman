#!/usr/bin/env gnuplot

file=system("echo $A")
set term png size 600, 400
set output "imgs/".file.".png"

set view map
set cbrange[0:0.2]
set xrange[0:600]
set yrange[0:400]
plot "data/".file.".dat" binary array=600x400 format="%f" with image
