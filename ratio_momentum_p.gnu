#creates a ratio plot of results from simulation with the ECCE internal note (for dp/p in bins of p, ie momentum_p.gnu)
#
#some helpful resources for multiplotting in gnuplot:
#
#https://livebook.manning.com/book/gnuplot-in-action-second-edition/appendix-e/19
#http://gnuplot.sourceforge.net/demo/layout.html

#---- L E T ' S   B E G I N ----
#don't forget we need to combine the desired files, and then we can plot their ratio

reset
#c = "barrellradii-2"

set terminal pngcairo size 2500, 1500 font "Times" enhanced
set size 1,1
set origin 0,0
set output "plots/".c."/ratio_momentum_p_.png"

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
set yrange [0: 2.5] 

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan30-n25.dat datafiles/data_mom_res_ECCE_intnote_etan30-n25.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio with ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan30-n25.dat datafiles/'.c.'/data_mom_res_PWG_req_etan30-n25.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#2. eta = n25-n20
unset label
set label "-25 < {/Symbol h} < -20" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan25-n20.dat datafiles/data_mom_res_ECCE_intnote_etan25-n20.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan25-n20.dat datafiles/'.c.'/data_mom_res_PWG_req_etan25-n20.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#3. eta = n20-n15
unset label
set label "-20 < {/Symbol h} < -15" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan20-n15.dat datafiles/data_mom_res_ECCE_intnote_etan20-n15.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan20-n15.dat datafiles/'.c.'/data_mom_res_PWG_req_etan20-n15.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#4. eta = n15-n10
unset label
set label "-15 < {/Symbol h} < -10" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan15-n10.dat datafiles/data_mom_res_ECCE_intnote_etan15-n10.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan15-n10.dat datafiles/'.c.'/data_mom_res_PWG_req_etan15-n10.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#5. eta = n10-n05
unset label
set label "-10 < {/Symbol h} < -5" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan10-n05.dat datafiles/data_mom_res_ECCE_intnote_etan10-n05.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan10-n05.dat datafiles/'.c.'/data_mom_res_PWG_req_etan10-n05.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#6. eta = n05-00
unset label
set label "-5 < {/Symbol h} < 0" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_etan05-00.dat datafiles/data_mom_res_ECCE_intnote_etan05-00.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_etan05-00.dat datafiles/'.c.'/data_mom_res_PWG_req_etan05-00.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#7. eta = 00-05
unset label
set label "0 < {/Symbol h} < 5" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta00-05.dat datafiles/data_mom_res_ECCE_intnote_eta00-05.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta00-05.dat datafiles/'.c.'/data_mom_res_PWG_req_eta00-05.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#8. eta = 05-10
unset label
set label "5 < {/Symbol h} < 10" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta05-10.dat datafiles/data_mom_res_ECCE_intnote_eta05-10.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta05-10.dat datafiles/'.c.'/data_mom_res_PWG_req_eta05-10.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#9. eta = 10-15
unset label
set label "10 < {/Symbol h} < 15" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta10-15.dat datafiles/data_mom_res_ECCE_intnote_eta10-15.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta10-15.dat datafiles/'.c.'/data_mom_res_PWG_req_eta10-15.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#10. eta = 15-20
unset label
set label "15 < {/Symbol h} < 20" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta15-20.dat datafiles/data_mom_res_ECCE_intnote_eta15-20.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta15-20.dat datafiles/'.c.'/data_mom_res_PWG_req_eta15-20.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#11. eta = 20-25
unset label
set label "20 < {/Symbol h} < 25" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta20-25.dat datafiles/data_mom_res_ECCE_intnote_eta20-25.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta20-25.dat datafiles/'.c.'/data_mom_res_PWG_req_eta20-25.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

#12. eta = 25-30
unset label
set label "25 < {/Symbol h} < 30" at graph 0.7, 0.2 center font "Times, 15"

set key bottom right box opaque spacing 1 width 1 font "Times, 12"

plot '< paste datafiles/'.c.'/data_mom_res_sim_eta25-30.dat datafiles/data_mom_res_ECCE_intnote_eta25-30.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 5 ps 1.5 lc rgb '#0099AD' title 'Ratio of Data and ECCE Internal Note', '< paste datafiles/'.c.'/data_mom_res_sim_eta25-30.dat datafiles/'.c.'/data_mom_res_PWG_req_eta25-30.dat' using 1:($2/$4) w linespoints lw 3 lt 1 pt 7 ps 1.5 lc rgb '#4DBEEE' title 'Ratio with PWG Requirements'

unset multiplot
unset output
