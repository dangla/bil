set term postscript eps enhanced "newcenturyschlbk-italic" 20
set output 's.eps'
set size square 7.5/10.,1.
set ylabel 's' "newcenturyschlbk-italic,20"
set xlabel 't'
set key 0.2,-0.8
set nokey
plot 'm7.p1' us 1:3  w l