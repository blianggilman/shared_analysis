# creates plot for fits over an eta bin
#
#some helpful resources for multiplotting in gnuplot:
#
#https://livebook.manning.com/book/gnuplot-in-action-second-edition/appendix-e/19
#http://gnuplot.sourceforge.net/demo/layout.html

set mulitplot layout 5,2 title "Fits" font ",14"

set tmargin 2
set title "Plot 1"
set key right top Left title 'Legend' box 3

plot "ENTER.dat" using 1:2 ti "ENTER.dat"
#also plot std deviation


