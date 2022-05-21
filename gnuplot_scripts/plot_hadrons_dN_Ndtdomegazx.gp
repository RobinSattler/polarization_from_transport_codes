# this script needs gnuplot > 5 and must be used as:
# gnuplot -c inputfile

set term png enh font "Helvetica, 20" size 900,750
inputfile = ARG1
outputfile = inputfile.".png"
set out outputfile

set xlabel "{/Symbol w}_{zx}"
set ylabel "t (fm)"
set cblabel "dN^2/(N d{/Symbol w}_{zx} dt)"

set size square 
set pm3d
unset surface  # don't need surfaces
set view map
#set contour
set key outside
set logscale zcb
set format cb "10^{%L}" 
set cbrange [1e-4:0.5]
set xrange [-0.375:0.375]
set cblabel offset 1
#set cntrparam cubicspline  # smooth out the lines
set pm3d interpolate 0,0# interpolate the color
set pal defined (1 '#00008f', 8 '#0000ff', 24 '#00ffff', 40 '#ffff00', 56 '#ff0000', 64 '#800000')
#splot inputfile u 1:2:(($3)+1.e-14) notitle with lines lt 1 lw 0
set yrange [0:45]
splot inputfile u 1:2:(($3)) notitle with lines lt 1 lw 0
