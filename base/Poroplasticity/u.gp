set term postscript eps enhanced "helvetica"
set output 'u.eps'
set size square 7.5/10.,1.
set xlabel 'r (m)' "helvetica,20"
set ylabel 'u (m)' "helvetica,20"
set key bottom right
set nologscale  
plot 't1.3.t1' us 1:5 title 't = 1.5e6 s' w l 1\
    ,'t1.3.t2' us 1:5 title 't = 50.e6 s' w l 1\
    ,'t1.3.t3' us 1:5 title 't = 300.e6 s' w l 1