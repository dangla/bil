set term postscript eps enhanced "Helvetica" 20
set output 'sl.eps'
set size square 0.75,1.
set xlabel 'r (m)' "Helvetica,20"
set ylabel 'S_l' "Helvetica,20"
set ytics 0,0.05,1
plot 't1.3.t1' us 1:30 title 't = 1.5e6 s' w l 1\
    ,'t1.3.t2' us 1:30 title 't = 50.e6 s' w l 1\
    ,'t1.3.t3' us 1:30 title 't = 300.e6 s' w l 1