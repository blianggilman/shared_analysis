# This shell script will run all the steps of the momentum analysis maybe?
# To set it up for the first time, you might have to run: chmod +x runMomResAnalysis.sh
# For regular setup, make sure you are in docker container w setup environment initialized

# ***** To run, do: ./runMomResAnalysis.sh [arg1] [arg2] *****

# [arg1] is the job number, ex. [arg1] = 60218804
# [arg2] is the configuration you are running (will go into filenames), ex. [arg2] = baseline, or [arg1] = fwddisk1_to_10m

# create directory to store plots with fits
cd plots/
if [ -d $2 ]; then
    echo "$2 is a directory that already exists"
else
    mkdir $2
fi
cd ..

cd datafiles/
if [ -d $2 ]; then
    echo "$2 is a directory that already exists"
else mkdir $2
fi
cd ..

# replace configuration name and input file, then run root file
sed -i "s/baseline/$2/" insert.h
sed -i "s/00000000/$1/" insert.h
root plotMomentumResolution.C

sed -i "s/$2/baseline/" insert.h
sed -i "s/$1/00000000/" insert.h

# run gnuplot code to generate plots here
#!/bin/bash

gnuplot -persist <<-EOFMarker
	c = "$2"
	otherc = "new_geo" #need to change this!!!
	othercname = "Ernst + Rey's geometry" #change this too!!!
	load "momentum_p.gnu"
	load "ratio_momentum_p.gnu"
	load "compare_geometries_p.gnu
EOFMarker

#done!
