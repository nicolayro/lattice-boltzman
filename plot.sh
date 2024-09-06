#!/usr/bin/env gnuplot

files=system("find data/*.dat | awk -F '[/.\\n]' '{printf $2 \"\\n\"}'")
do for [file in files] {
    print file
    set term png
    set output "imgs/".file.".png"
    set view 0,0,1
    set cbrange[0:1]
    set palette defined (0"black", 12"cyan", 16"white")
    splot "data/".file.".dat" binary matrix with pm3d
}
