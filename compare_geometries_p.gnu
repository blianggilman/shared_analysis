#creates a plot that compares the result of one geometry with another (default is set to new geometry by Ernst and Rey)
#
##some helpful resources for multiplotting in gnuplot:
#
#https://livebook.manning.com/book/gnuplot-in-action-second-edition/appendix-e/19
#http://gnuplot.sourceforge.net/demo/layout.html

#---- L E T ' S   B E G I N ----

reset
#c = "barrellradii-2"
#otherc = "" read in from shell script
#othercname also referenced in shell script 

set terminal pngcairo size 2500, 1500 font "Times" enhanced
set size 1,1
set origin 0,0
set output "plots/".c."/comparison_with_newgeo_.png"

set xtics font "Times, 15"
set ytics font "Times, 15"

set multiplot layout 3, 4 title '{/:Bold Ratio Plot of Momentum Resolution with ECCE Internal Note}' font "Times, 25"

set xlabel 'Track p [GeV/c]' font "Times, 15"
#set ylabel 'dp/p' font "Times, 15"

unset key

#read in dat file and plot

#1. eta = n30-n25
unset label
set label "-30 < {/Symbol h} < -25" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0: 0.035]

plot 'datafiles/'.c.'/data_mom_res_sim_etan30-n25.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan30-n25.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan30-n25.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#2. eta = n25-n20
unset label

set label "-25 < {/Symbol h} < -20" at graph 0.75, 0.75 center font "Times, 15"

set key top right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.035]


plot 'datafiles/'.c.'/data_mom_res_sim_etan25-n20.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan25-n20.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan25-n20.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#3. eta = n20-n15
unset label
set label "-20 < {/Symbol h} < -15" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_etan20-n15.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan20-n15.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan20-n15.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#4. eta = n15-n10
unset label
set label "-15 < {/Symbol h} < -10" at graph 0.75, 0.75 center font "Times, 15"

set key top right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.035]

plot 'datafiles/'.c.'/data_mom_res_sim_etan15-n10.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan15-n10.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan15-n10.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#5. eta = n10-n05
unset label
set label "-10 < {/Symbol h} < -5" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_etan10-n05.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan10-n05.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan10-n05.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#6. eta = n05-00
unset label
set label "-5 < {/Symbol h} < 0" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_etan05-00.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_etan05-00.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_etan05-00.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#7. eta = 00 - 05
unset label
set label "0 < {/Symbol h} < 5" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_eta00-05.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta00-05.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta00-05.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#8. eta = 05 - 10
unset label
set label "5 < {/Symbol h} < 10" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_eta05-10.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta05-10.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta05-10.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#9. eta = 10 - 15
unset label
set label "10 < {/Symbol h} < 15" at graph 0.75, 0.75 center font "Times, 15"

set key top right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.025]

plot 'datafiles/'.c.'/data_mom_res_sim_eta10-15.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta10-15.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta10-15.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#10. eta = 15 - 20
unset label
set label "15 < {/Symbol h} < 20" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.015]

plot 'datafiles/'.c.'/data_mom_res_sim_eta15-20.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta15-20.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta15-20.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#11. eta = 20 - 25
unset label
set label "20 < {/Symbol h} < 25" at graph 0.75, 0.75 center font "Times, 15"

set key top right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.025]

plot 'datafiles/'.c.'/data_mom_res_sim_eta20-25.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta20-25.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta20-25.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

#12. eta = 25 - 30
unset label
set label "25 < {/Symbol h} < 30" at graph 0.75, 0.25 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

set xrange [0 : 22]
set yrange [0 : 0.035]

plot 'datafiles/'.c.'/data_mom_res_sim_eta25-30.dat' w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#09AD00' title 'New Trial', 'datafiles/data_mom_res_ECCE_intnote_eta25-30.dat' w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'ECCE Internal Note', 'datafiles/'.otherc.'/data_mom_res_sim_eta25-30.dat' w l lw 3 lt 1 lc rgb '#4DBEEE' title ''.othercname.''

unset multiplot
unset output

