set term postscript eps enhanced "Helvetica" 20
set output 'sl.eps'
set size square 0.75,1.
set ylabel 'z (m)' "Helvetica,24"
set xlabel 'S_l' "Helvetica,24"
set xtics 0,0.05,1
set key samplen 1 spacing 4
set nokey
plot 'm2.t0' us 8:1 title '0' w l 1\
    ,'m2.t1' us 8:1 title '1 an ' w l 1\
    ,'m2.t2' us 8:1 title '2 ans' w l 1\
    ,'m2.t3' us 8:1 title '4 ans' w l 1\
    ,'m2.t4' us 8:1 title '6 ans' w l 1\
    ,'m2.t5' us 8:1 title '8 ans' w l 1\
    ,'m2.t6' us 8:1 title '10 ans' w l 1\
    ,'m2.t7' us 8:1 title '20 ans' w l 1\
    ,'m2.t8' us 8:1 title '40 ans' w l 1\
    ,'m2.t9' us 8:1 title '50 ans' w l 1\
    ,'m2.t10' us 8:1 title '100 ans' w l 1\
    ,'m2.t11' us 8:1 title '100 ans' w l 1