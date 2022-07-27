#Plots momentum resolution dp/p in bins of momentum p, and compares with PWG requirements AND ECCE internal note
#
#Set overall margins for the combined set of plots and size them
# to generate a requested inter-plot spacing
#
#Finally, some helpful resources for multiplotting in gnuplot:
#
#https://livebook.manning.com/book/gnuplot-in-action-second-edition/appendix-e/19
#http://gnuplot.sourceforge.net/demo/layout.html

set multiplot layout 3, 4 title "Momentum Resolution in p bins, compared with PWG requirements and ECCE internal note" font ",14"

set tmargin 2
unset key

#read in dat file and plot

plot 'ENTER_HERE.dat' using 1:2 ti 'ENTER_HERE.dat'

 
