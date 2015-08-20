set output 'panel_plot.eps'
set terminal postscript enhanced color linewidth 2 font 20 

set style fill solid 0.3

set size 4,2

set multiplot layout 4,2

set size 1,1 
set xlabel "E (GeV)"
set ylabel "Freq." offset 1
set xrange [0:2]

set origin 0,1

plot 'all.dat' u ((floor($1/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w boxes ls 1 lc 1 lw 2 ti "Total energy", 'all.dat' u ((floor($4/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w p ls 2 lc 3 lw 2 ps 1.5 ti "Initial sterile energy"

set size 1,1 
set xlabel "Cos{/Symbol q}"
set ylabel "Freq." offset 1
set xrange [-1:1]

set origin 1,1

plot 'all.dat' u ((floor(cos($2*3.1415/180.0)/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w boxes ls 1 lc 1 lw 2 ti "Average direction", 'all.dat' u ((floor(cos($5*3.1415/180.0)/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w p ls 2 lc 3 lw 2 ps 1.5 ti "Initial sterile direction"

set size 1,1 
set xlabel "E (GeV)"
set ylabel "Freq." offset 1
set xrange [0:2]

set origin 0,0

plot 'all.dat' u ((floor($6/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w boxes ls 1 lc 1 lw 2 ti "Highest electron energy", 'all.dat' u ((floor($8/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w p ls 2 lc 3 lw 2 ps 1.5 ti "Lowest electron energy"

set size 1,1 
set xlabel "Cos{/Symbol q}"
set ylabel "Freq." offset 1
set xrange [-1:1]

set origin 1,0

plot 'all.dat' u ((floor(cos($7*3.1415/180.0)/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w boxes ls 1 lc 1 lw 2 ti "Highest-energy electron direction", 'all.dat' u ((floor(cos($9*3.1415/180.0)/0.1) + 0.5)*0.1):(1.0/200000) smooth freq w p ls 2 lc 3 lw 2 ps 1.5 ti "Lowest-energy electron direction"

#unset multiplot


#set output 'panel_plot2.eps'
#set terminal postscript enhanced color linewidth 2 font 20

#set style fill solid 0.3

#set size 2,2
#set multiplot layout 2,2

set size 1,1 
set xlabel "Angular separation (deg.)"
set ylabel "Freq." offset 1
set xrange [0:180]

#set origin 0,1
set origin 2,1

plot 'all.dat' u ((floor($3/10) + 0.5)*10):(1.0/200000) smooth freq w boxes ls 1 lc 1 lw 2 ti "True separation", 'all.dat' u ((floor($10/10) + 0.5)*10):(1.0/200000) smooth freq w p ls 2 lc 3 lw 2 ps 1.5 ti "Foreshortened separation"

set size 1,1
set xlabel "Total Energy (GeV)"
set ylabel "Electron energy ratio (low/high)" offset 1
set xrange [0:2]
set yrange [0:1]
set tics nomirror out scale 0.75
set border 3 front
set cbtics scale 0
load 'parula.pal'

#set origin 1,1
set origin 2,0

plot 'Esum_EnergyRatio.dat' u ($2*0.01):($1*0.01):3 matrix w image notitle

set size 1,1
set xlabel "Total Energy (GeV)"
set ylabel "Foreshortened angular separation" offset 1
set xrange [0:2]
set yrange [0:180]
set tics nomirror out scale 0.75
set border 3 front
set cbtics scale 0
load 'parula.pal'

#set origin 0,0
set origin 3,1

plot 'Esum_FSangularsep.dat' u ($2*0.01):($1*1.0):3 matrix w image notitle

set size 1,1
set xlabel "True angular separation"
set ylabel "Foreshortened angular separation" offset 1
set xrange [0:180]
set yrange [0:180]
set tics nomirror out scale 0.75
set border 3 front
set cbtics scale 0
load 'parula.pal'

#set origin 1,0
set origin 3,0

plot 'Angularsep_FSangularsep.dat' u ($2*1.0):($1*1.0):3 matrix w image notitle


unset multiplot
