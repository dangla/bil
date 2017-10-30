# Coordinates(1) Pore pressure(4) Displacements(5) Fluid mass flow(8) Stresses(11) Plastic strains(20) Cumulative plastic shear strain(29) Permeability(30)

# Input Data file
file = 'Poroplasticity'

# To produce eps file
set term postscript eps enhanced color 20


# Size
set size square 1.,1.




# x-axis
set xlabel 'Radius (m)' font ",24"
x0 = 0
x1 = 10
set xrange [x0:x1] noreverse nowriteback
set xtics 2
#set format x "%1.0e"


# Labels
#set arrow from 0.3,0.4 to 0.4,0.4
#set label "p_c" at 0.28,0.4 right font "newcenturyschlbk-italic,28"
#set label "1" at 0.9,0.15 right
#set label "1.5" at 0.9,0.24 right
#set label "2" at 0.9,0.32 right
#set label "3" at 0.9,0.45 right


# Styles
set pointsize 2
set style line 10  lw 4 lt 1 pi -10 ps 2 pt 4  lc rgb "black"
set style line 20  lw 4 lt 1 pi -10 ps 2 pt 8  lc rgb "black"
set style line 30  lw 4 lt 1 pi -10 ps 2 pt 6  lc rgb "black"
set style line 40  lw 4 lt 1 pi -10 ps 2 pt 10 lc rgb "black"
set style line 50  lw 4 lt 1 pi -10 ps 2 pt 12 lc rgb "black"
set style line 60  lw 4 lt 1 pi -10 ps 2 pt 14 lc rgb "black"

set style line 100 lw 1 lt -1 lc rgb "black"
set style line 200 lw 2 lt -1 lc rgb "black"
set style line 300 lw 3 lt -1 lc rgb "black"

set style fill solid 0.7 border lt -1
set style fill pattern 6 border lt -1
set style rectangle fillcolor rgbcolor "#BBBBBB" fs pattern 4 border lt -1
set contour both



# Data files
file0 = file.'.t0'
file1 = file.'.t1'
file2 = file.'.t2'
file3 = file.'.t3'
file4 = file.'.t4'
fileI = file4


# 1st Plot
# --------
set output 'RadialDisplacement.eps'

# Legends
set key bottom right reverse Left samplen 2 spacing 2

# Size
set size square 1.,1.

# y-axis
set ylabel 'Radial displacement (m)' font ",24"
#set yrange[0:1.2] noreverse nowriteback
#set ytics norotate 0,0.2,1.2


plot \
     file1 us 1:5 w l ls 10 title 't = 1.5e6 s' \
    ,file2 us 1:5 w l ls 20 title 't = 50.e6 s' \
    ,file3 us 1:5 w l ls 30 title 't = 300.e6 s'


# 2nd Plot
# --------
set output 'PorePressure.eps'

# Legends
set key bottom right

# Size
set size square 1.,1.

# y-axis
set ylabel 'Pore pressure (Pa)' font ",24"
#set yrange[0:1.2] noreverse nowriteback
#set ytics norotate 0,0.2,1.2


plot \
     file1 us 1:4 w l ls 10 title 't = 1.5e6 s' \
    ,file2 us 1:4 w l ls 20 title 't = 50.e6 s' \
    ,file3 us 1:4 w l ls 30 title 't = 300.e6 s'



# 3rd Plot
# --------
set output 'VolumetricPlasticStrain.eps'

# Legends
set key bottom right

# Size
set size square 1.,1.

# y-axis
set ylabel 'Volumetric plastic strain (-)' font ",24"

plot \
     file1 us 1:($20+$24+$28) w l ls 10 title 't = 1.5e6 s' \
    ,file2 us 1:($20+$24+$28) w l ls 20 title 't = 50.e6 s' \
    ,file3 us 1:($20+$24+$28) w l ls 30 title 't = 300.e6 s'
    
    

# 4th Plot
# --------
set output 'EffectiveHoopStress.eps'

# Legends
set key bottom right

# Size
set size square 1.,1.

# y-axis
set ylabel 'Effective hoop stress (Pa)' font ",24"

plot \
     file1 us 1:($19+0.8*$4) w l ls 10 title 't = 1.5e6 s' \
    ,file2 us 1:($19+0.8*$4) w l ls 20 title 't = 50.e6 s' \
    ,file3 us 1:($19+0.8*$4) w l ls 30 title 't = 300.e6 s'
