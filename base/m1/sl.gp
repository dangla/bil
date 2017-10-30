set term postscript eps enhanced "Helvetica" 20
set output 'sl.eps'
set size square 0.75,1.
set xlabel 'r (m)' "Helvetica,24"
set ylabel 'S_l' "Helvetica,24"
set xrange [0:4]
set yrange [0.75:1]
set ytics 0,0.05,1
set key samplen 1 spacing 4
plot 'm1.t0' us 1:6 title '0' w l 1\
    ,'m1.t1' us 1:6 title '1 an ' w l 1\
    ,'m1.t2' us 1:6 title '2 ans' w l 1\
    ,'m1.t3' us 1:6 title '4 ans' w l 1\
    ,'m1.t4' us 1:6 title '6 ans' w l 1\
    ,'m1.t5' us 1:6 title '8 ans' w l 1\
    ,'m1.t6' us 1:6 title '10 ans' w l 1\
    ,'m1.t7' us 1:6 title '20 ans' w l 1\
    ,'m1.t8' us 1:6 title '40 ans' w l 1\
    ,'m1.t9' us 1:6 title '50 ans' w l 1\
    ,'m1.t10' us 1:6 title '100 ans' w l 1