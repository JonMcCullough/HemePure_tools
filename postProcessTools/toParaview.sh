#!/usr/bin/env bash

processing () {

	NAME=${1?No file provided}
	OUT=${2?No output name provided}

	echo Processing $NAME for ParaView, to be called $OUT.

	sed -i '1,2d' $1
	sed -i '/^$/d' $1

	#awk -F ' ' 'BEGIN{N=-1; last_step=-1}{if(last_step != $1) {N++; print>$OUT"N".txt"; last_step=$1;} else {print>>$OUT"N".txt"}}' $NAME
	awk -v out=$2 -F ' ' 'BEGIN{N=-2; last_step=-1}{if(last_step != $1) {N++; print>out N".txt"; last_step=$1;} else {print>>out N".txt"}}' $NAME

	sed -i '1isteps gridX gridY gridZ velX velY velZ pressure coverageFactor' $OUT*
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

