#!/usr/bin/gnuplot -persist
#
reset
set term postscript enhanced "Times New Roman"  26 color solid
set termoption dash
set output "Kasolid10_p1m1.eps"
set lmargin 5
set bmargin 3.5
set rmargin 6
set tmargin 1
set xtics 0,1
set ytics 
set mytics 2 
set ytics nomirror
set y2tics nomirror 
set y2tics 0,0.1
set my2tics 2
set y2tics nomirror  
set title ""
set title  offset character 0, 0, 0 font "" norotate
set key title "Method 1 (10 days)"
set key right
set xlabel "Distance to the surface (mm)" 
set ylabel "Solid concentration (mol/L)" offset 2,0
set y2label "Porosity"  offset -2,0
set xrange [ 0 : 10 ] noreverse nowriteback  # (currently[0.151196:0.304308] )
set yrange [ * : 25 ] noreverse nowriteback    # (currently [0.411707:0.459901] )
set y2range [ 0 : 1 ] noreverse nowriteback   # (currently [0.411707:0.459901] )
set border  lw 2

set mxtics 2
set label "Attack direction" font ",22" at first 0.4, first 23
set arrow from 0.4,22 to 3,22  lt 1 lc 0 lw 2

plot 'Kap1m1.t2' u (100-$1*100):21 w linespoints lt 0 lc 0 pt 4 pi 3. lw 2 axis x1y1 title "Portlandite",\
     'Kap1m1.t2' u (100-$1*100):22 w linespoints lt 2 lc 0 pt 9 pi 3. lw 2 axis x1y1 title "Gypsum",\
     'Kap1m1.t2' u (100-$1*100):23 w linespoints lt 4 lc 0 pt 7 pi 3. lw 2 axis x1y1 title "C-S-H",\
     'Kap1m1.t2' u (100-$1*100):24 w l lt 1 lc 0  lw 2 axis x1y2 title "Porosity"
reset        
#    EOF
