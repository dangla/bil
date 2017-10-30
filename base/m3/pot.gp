set term postscript eps enhanced "Helvetica" 20
set output 'pot.eps'
set size square 0.75,1.
set xlabel 'z' "Helvetica,24"
set ylabel '{/Symbol f}' "Helvetica,24"
set nokey
plot 'm3.t0' us 1:4 w l 1\
    ,'m3.t1' us 1:4 w l 1\
    ,'m3.t2' us 1:4 w l 1\
    ,'m3.t3' us 1:4 w l 1