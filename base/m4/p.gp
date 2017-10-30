set term postscript eps enhanced "Helvetica" 20
set output 'p.eps'
set size square 0.75,1.
set xlabel 't (heures)' "Helvetica,24"
set ylabel 'pressions (MPa)' "Helvetica,24"
#set yrange[12:15]
#set xtics 0,0.2,1
set ytics 0,1,20
set key samplen 1 spacing 4
plot 'm4.p1' us ($1/3600):(($2+$3)*1.e-6) title 'p' w l 1\
    ,'m4.p1' us ($1/3600):($2*1.e-6) title 'p_e' w l 2\
    ,'m4.p1' us ($1/3600):($3*1.e-6) title 'p_s' w l 2
reset
set output 'u.eps'
set size square 0.8,1.
set xlabel 't (heures)' "Helvetica,24"
set ylabel 'deplacement vertical (mm)' "Helvetica,24"
#set ytics 0,0.1,1
set key right bottom samplen 1 spacing 4
plot 'm4.p2' us ($1/3600):($4*1.e3) notitle w l 1