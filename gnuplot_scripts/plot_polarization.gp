set term pos eps col enh font "Helvetica, 22"

set format y "10^{%T}"

my_file="AuAu_ecm_7_7_b6_time_distr.txt"
comp_file="fig7a_art_vitiuk_bravina.dat"

set out "polarization_time_comparison.eps"
set key spacing 1.5
set logscale y
set mxtics 5
set xrange [0:30]
set ylabel "P_y"
set xlabel "t [fm]"
set label "0.1<p_T<3 GeV, |y|<1" at 14,0.1 tc "black"
plot my_file u 1:(-2*($5)) w p ps 1.2 pt 5 lc "blue" title "My data", comp_file w p ps 1.2 pt 7 lc "red" title "1910.06292"
unset label
