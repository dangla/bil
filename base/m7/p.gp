set term postscript eps enhanced "newcenturyschlbk-italic" 16
set output 'p.eps'
set size square 7.5/10.,1.
set xlabel 'p'
set ylabel 'z'
set label "t = 0.001" at 0.8,-0.03
set label "t = 0.01"  at 0.8,-0.14
set label "t = 0.1"   at 0.7,-0.4
set label "t = 0.25"  at 0.55,-0.5
set label "t = 0.5"   at 0.37,-0.6
set label "t = 1"     at 0.15,-0.78
set key 0.2,-0.8
set nokey
plot 'm7.t1' us 4:1 title 't = 0.001' w l\
    ,'m7.t2' us 4:1 title 't = 0.01' w l\
    ,'m7.t3' us 4:1 title 't = 0.1' w l\
    ,'m7.t4' us 4:1 title 't = 0.25' w l\
    ,'m7.t5' us 4:1 title 't = 0.5' w l\
    ,'m7.t6' us 4:1 title 't = 1.' w l
