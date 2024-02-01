#!/usr/bin/env bash

processing () {

	NAME=${1?No file provided}
	OUT=${2?No output name provided}

	echo Processing $NAME for ParaView, to be called $OUT.

	HEADER='1isteps gridX gridY gridZ velX velY velZ pressure coverageFactor'

	# Delete all empty lines
	sed -i '/^$/d' $1

	# Separate time series into individual files
	awk -v var=$OUT -v header="$HEADER" -F ' ' \
	'
    	BEGIN{N=-3; last_step=-1}
    	{
	        if(last_step != $1)
        	{
            	N++;
            	print header > var N".txt";
            	print >> var N".txt";
            	last_step=$1;
        	}
        	else
        	{
            	print >> var N".txt"
        	}
    	}
	' $NAME

 	rm $OUT\-1.txt $OUT\-2.txt
}

rm results_0/outlet*
rm results_0/inlet*
rm results_0/whole*
rm results_0/planeX*

./hemeXtract -X results_0/Extracted/outlet.dat > outletFirst.txt
./hemeXtract -X results_0/Extracted/inlet.dat >  inletFirst.txt
./hemeXtract -X results_0/Extracted/whole.dat >  wholeFirst.txt
./hemeXtract -X results_0/Extracted/planeX.dat >  planeXFirst.txt

processing "outletFirst.txt" "outletFirstPV"
processing "inletFirst.txt"  "inletFirstPV"
processing "wholeFirst.txt"  "wholeFirstPV"
processing "planeXFirst.txt"  "planeXFirstPV"

mv *First* results_0

