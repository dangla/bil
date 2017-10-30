set size square 0.75,1.
set xlabel 'x (m)' "Helvetica,24"
set ylabel 'S_g' "Helvetica,24"
set yrange[0:1]
set ytics 0,0.2,1
set key samplen 1 spacing 4
n=1
plot 'ex51.t0' us 1:6 title '0' w l 1
pause n
plot 'ex51.t1' us 1:6 title '1 heure ' w l 1
pause n
plot 'ex51.t2' us 1:6 title '2 heures' w l 1
pause n
plot 'ex51.t3' us 1:6 title '3 heures' w l 1
pause n
plot 'ex51.t4' us 1:6 title '4 heures' w l 1
pause n
plot 'ex51.t5' us 1:6 title '5 heures' w l 1
pause n
plot 'ex51.t6' us 1:6 title '6 heures' w l 1
pause n
plot 'ex51.t7' us 1:6 title '7 heures' w l 1
pause n
plot 'ex51.t8' us 1:6 title '8 heures' w l 1
pause n
plot 'ex51.t9' us 1:6 title '9 heures' w l 1
pause n
plot 'ex51.t10' us 1:6 title '10 heures' w l 1