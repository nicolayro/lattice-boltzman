#!/usr/bin/env gnuplot

file=system("echo $A")
set term png size 600, 400
set output "imgs/".file.".png"
set view 0,0,1
set cbrange[0:0.0006]
set palette defined (0"black", 12"cyan", 16"white")
splot "data/".file.".dat" binary matrix with pm3d
