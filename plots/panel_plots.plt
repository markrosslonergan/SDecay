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

set palette negative defined ( \
    0 '#D53E4F',\
    1 '#F46D43',\
    2 '#FDAE61',\
    3 '#FEE08B',\
    4 '#E6F598',\
    5 '#ABDDA4',\
    6 '#66C2A5',\
    7 '#3288BD' )

f(x) = x>2 ? x : 1/x
g(x) = x<10 ? 0 : 1

set style fill solid 0.8 noborder

set size 1,1
set xlabel "Total Energy (GeV)"
set ylabel "Electron energy ratio (low/high)" offset 1
set xrange [0:2]
set yrange [0:1]

set origin 2,0

plot '<sort -n -k3 Esum_EnergyRatio.dat' u 1:2:(g($3)*0.01):(f($3)) w circle lc palette notitle

set size 1,1
set xlabel "Total Energy (GeV)"
set ylabel "Foreshortened angular separation" offset 1
set xrange [0:2]
set yrange [0:180]

#set origin 0,0
set origin 3,1

plot '<sort -n -k3 Esum_FSangularsep.dat' u 1:2:(g($3)*0.01):(f($3)) w circle lc palette notitle

set size 1,1
set xlabel "True angular separation"
set ylabel "Foreshortened angular separation" offset 1
set xrange [0:180]
set yrange [0:180]

#set origin 1,0
set origin 3,0

plot '<sort -n -r -k3 Angularsep_FSangularsep.dat' u 1:2:(g($3)):(f($3)) w circle lc palette notitle


unset multiplot
