set term postscript eps enhanced color 20

set termoption dash

# Size
set size square 0.75,1.


# Legends
set key title ""
set key top right #at first 90,first 12.5

# x-axis
set xlabel "Time (day)" font  ",24"
set xtics 0,50
set xrange [ 0 : 350 ] noreverse nowriteback

# y-axis
set ylabel "Surface pH" font  ",24"
set ytics 
set ytics nomirror  
set yrange [ 7 : 13 ] noreverse nowriteback

# Styles
set pointsize 2
set style line 10  lw 2 lt 1 pi -10 ps 2 pt 1 pi 500 lc rgb "black"
set style line 20  lw 2 lt 1 pi -10 ps 2 pt 8 lc rgb "black"
set style line 30  lw 2 lt 1 pi -10 ps 2 pt 6 lc rgb "black"
set style line 11  lw 2 lt 0 pi -10 ps 2 pt 4 lc rgb "black"
set style line 21  lw 2 lt 0 pi -10 ps 2 pt 8 lc rgb "black"
set style line 31  lw 2 lt 0 pi -10 ps 2 pt 6 lc rgb "black"


# Plot 1
set output "phneu.eps"
plot 'Yuan1.p3' us ($1/86400):2  w linespoints ls 10 title "10 ppm"

exit

plot 'h2s_10ppm.p4' u ($1/86400):2   w linespoints lt 1 lc 0 pt 1 pi 500 lw 2 axis x1y1 title "10 ppm",\
     'h2s_15ppm.p4' u ($1/86400):2  w linespoints lt 1 lc 0 pt 3 pi 500 lw 2 axis x1y1 title "15 ppm",\
     'h2s_35ppm.p4' u ($1/86400):2  w linespoints lt 1 lc 0 pt 4 pi 500 lw 2 axis x1y1 title "35 ppm",\
     'h2s_100ppm.p4' u ($1/86400):2  w linespoints lt 1 lc 0 pt 6 pi 500 lw 2 axis x1y1 title "100 ppm",\
     'h2s_200ppm.p4' u ($1/86400):2   w linespoints lt 1 lc 0 pt 8 pi 500 lw 2 axis x1y1 title "200 ppm"
reset          
#    EOF
