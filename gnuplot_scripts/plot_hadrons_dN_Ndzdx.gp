# this script needs gnuplot > 5 and must be used as:
# gnuplot -c inputfile

set term png enh font "Helvetica, 20" size 900,750
inputfile = ARG1
outputfile = inputfile.".png"
set out outputfile

set xlabel "z (fm)"
set ylabel "x (fm)"
set cblabel "dN^2/(N dx dz)"

set size square 
set pm3d
set bmargin at screen 0.15
set rmargin at screen 0.87
unset surface  # don't need surfaces
set view map
#set contour
set key outside
set logscale zcb
set format cb "10^{%L}" 
set cbrange [1e-6:0.008]
set cblabel offset 1
#set cntrparam cubicspline  # smooth out the lines
#set pm3d interpolate 0,0# interpolate the color
set pm3d 
set pal defined (1 '#00008f', 8 '#0000ff', 24 '#00ffff', 40 '#ffff00', 56 '#ff0000', 64 '#800000')
splot inputfile u 1:2:(($3)) notitle with lines lt 1 lw 0
