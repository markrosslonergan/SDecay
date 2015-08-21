#!/bin/bash

SMASS=$1
ZMASS=$2

FOLDER=$SMASS"_"$ZMASS

mkdir -p data/$FOLDER

./inflight -m $SMASS -Z $ZMASS > data/all.dat

mv data/*.dat data/$FOLDER/.

cp plots/parula.pal data/$FOLDER
cp plots/panel_plots.plt data/$FOLDER
cp plots/plot.sh data/$FOLDER

cd data/$FOLDER

./plot.sh

mv panel_plot.pdf $SMASS"_"$ZMASS".pdf"
cp $SMASS"_"$ZMASS".pdf" ../../plots
