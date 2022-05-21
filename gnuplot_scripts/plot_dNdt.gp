set term pos eps col enh font "Helvetica, 22"

set format y "10^{%T}"

my_file="Lambda_time_distr_AuAuEcm7_b6_raplim1.dat"
comp_file="fig1_vitiuk.dat"

set out "Lambda_dN_over_dt_comparison.eps"
set key spacing 1.5
set logscale y
set mxtics 5
set ylabel "d<N>/dt"
set xlabel "t [fm]"
set label "0.1<p_T<3 GeV, |y|<1" at 40,0.00001 tc "black"
plot my_file u 1:3 w p ps 1.2 pt 5 lc "blue" title "My data", comp_file w p ps 1.2 pt 7 lc "red" title "1910.06292"
unset label
