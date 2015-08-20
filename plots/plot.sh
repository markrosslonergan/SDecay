#!/bin/bash


gnuplot panel_plots.plt
epstopdf --autorotate=All panel_plot.eps

