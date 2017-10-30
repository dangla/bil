set term postscript eps enhanced "Helvetica" 20
set output 'n.eps'
set size square 0.75,1.
set xlabel 'z' "Helvetica,24"
set ylabel 'n' "Helvetica,24"
plot 'm3.t1' us 1:5 w l 1\
    ,'m3.t1' us 1:6 w l 1