# Gunplot file. To produce the figure type
#
# gnuplot PhaseDiagramOfCaO-CO2-H2O.gp
#
#
# Options not affected by reset
#------------------------------


# To produce eps file
set term postscript eps enhanced color 20
# To produce animation
#set term gif animate delay 100

set termoption solid
set termoption lw 2
set termoption font ",20"



# Equilibrium constants (at 293 K)
# --------------------------------
LogK_CH     = -5.15
LogK_Calcite    = -8.45
LogK_Aragonite  = -8.30
LogK_Vaterite   = -7.87
LogK_h2o    = -14.17
LogK_co3    = -3.67
LogK_hco3   = -10.4
LogK_h2co3  =  2.77
Loga_h2o    =  0

Logq_co2_Cal    = LogK_Calcite   - LogK_CH + LogK_co3 + LogK_hco3 + LogK_h2co3 + Loga_h2o
Logq_co2_Vat    = LogK_Vaterite  - LogK_CH + LogK_co3 + LogK_hco3 + LogK_h2co3 + Loga_h2o
Logq_co2_Ara    = LogK_Aragonite - LogK_CH + LogK_co3 + LogK_hco3 + LogK_h2co3 + Loga_h2o

Logq_Cal(x) = (x < Logq_co2_Cal) ? 0 : Logq_co2_Cal - x
Logq_Vat(x) = (x < Logq_co2_Vat) ? 0 : Logq_co2_Vat - x
Logq_Ara(x) = (x < Logq_co2_Ara) ? 0 : Logq_co2_Ara - x


q_co2_Cal = 10**(Logq_co2_Cal)
q_co2_Vat = 10**(Logq_co2_Vat)
q_co2_Ara = 10**(Logq_co2_Ara)
q_Cal(x)  = 10**(Logq_Cal(x))




# Linetypes
set linetype 1 lw 1 linecolor rgb "red"
set linetype 2 lw 1 linecolor rgb "light-magenta"
set linetype 3 lw 1 linecolor rgb "purple"
set linetype 4 lw 1 linecolor rgb "steelblue"
set linetype 5 lw 1 linecolor rgb "aquamarine"
set linetype 6 lw 1 linecolor rgb "bisque"
set linetype 7 lw 1 linecolor rgb "bisque"
set linetype 8 lw 1 linecolor rgb "light-goldenrod"
set linetype 9 lw 1 linecolor rgb "light-goldenrod"
set linetype 10 lw 1 linecolor rgb "black"
#set linetype cycle 9


# Styles
#set pointsize 2
set style fill solid 0.5
set style line 1  lw 1 lt 1 pi -10 ps 2 pt 4 lc rgb "red"
set style line 2  lw 1 lt 2 pi -10 ps 2 pt 8 lc rgb "light-magenta"
set style line 3  lw 1 lt 3 pi -10 ps 2 pt 6 lc rgb "purple"
set style line 4  lw 1 lt 4 pi -10 ps 2 pt 4 lc rgb "steelblue"
set style line 5  lw 1 lt 5 pi -10 ps 2 pt 8 lc rgb "bisque"
set style line 6  lw 1 lt 6 pi -10 ps 2 pt 6 lc rgb "green"
set style line 6  lw 1 lt 7 pi -10 ps 2 pt 6 lc rgb "green"
set style line 7  lw 1 lt 8 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"
set style line 7  lw 1 lt 9 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"
set style line 11 lw 3 lt 0 lc rgb "black"




# Plot 1
#-------
# Size
set size square 1,1

set size square 0.75,1.1

set origin 0.,0.


# x-axis
#set xlabel '{/Symbol r}_{CO_2}/{/Symbol r}@_{CO_2}^{CH}' font ",26"
set xlabel 'Activity of {CO_2}' font ",26"
set xtics 1
#set mxtics 5
#set format x "10^{%T}"
#set logscale x
x_min = 1.e-2*q_co2_Cal
x_max = 1.e4*q_co2_Cal
x_min = -17
x_max = -10
set xrange[x_min:x_max]
set format x "10^{%+-3.0f}"


# y-axis
set ylabel 'Saturation index of portlandite' font ",26"
set ytics 1
#set mytics 5
#set format y "10^{%T}"
#set logscale y
y_min = -4
y_max = 1
set yrange[y_min:y_max]
set format y "10^{%+-3.0f}"



#set label "n_{CH} {/Symbol \263} 0" at 1.e-1,6.e-2
#set label "n_{CH} = 0" at 1.3,6.e-2
#set label "n_{CC} = 0" at 1.e-1,2.e-2
#set label "n_{CC} {/Symbol \263} 0" at 1.3,2.e-2
#set label "Portlandite" at Logq_co2_Cal,0.25527 right
#set label "Calcite" at (2.7 + Logq_co2_Cal),-2.3 center rotate by -50
set label "Portlandite" at -16,0.255 center font ",22"
set label "Calcite"     at -13,-1.3 center rotate by -54 font ",22"
set label "Vaterite"    at -12,-1.7 center rotate by -54 font ",22"

#set arrow from 1.35,0.74 to 2,1 nohead lw 1 lt 0
#set arrow from 1.84,0.539 to 3,1 nohead lw 1 lt 0
#set arrow from 2.5,0.4 to 4,1 nohead lw 1 lt 0
#set arrow from 3.26,0.3 to 4,0.6 nohead lw 1 lt 0

#set arrow from 4,0.6 to 2,0.5 head lw 1 lt 1
#set label "Q_{CC}/K_{CC} = 1" at 5,0.6 front

#set arrow from 1*q_co2_Cal,y_min to 1,1 nohead lw 2 lt 0
#set arrow from x_min,1 to 1,1 nohead lw 4 lt 1



# Legends
#set key lmargin  horizontal Left reverse samplen 2 spacing 1
#set key lmargin  horizontal Right samplen 2 spacing 1
set key top right reverse Left spacing 2 samplen 1

set output 'PhaseDiagramOfCaO-CO2-H2O.eps'
#set samples 500

plot \
     Logq_Vat(x)  with lines ls 11 lw 3 notitle \
    ,Logq_Cal(x)  with lines lt 10 lw 3 notitle \






exit

plot q_Cal(x)  with filledcurves below y1=y_max linecolor rgbcolor "red" notitle\
    ,q_Cal(x)  with lines lw 3 lt -1 notitle

exit

set parametric

r(q_Cal,q_cc) = q_cc/q_Cal
t_min = 1.e-6
t_max = 1.e2
set trange[t_min:t_max]

plot r(t,1),(r(t,1) < 1) ? 1 : t  with filledcurves below y1=y_max linecolor rgbcolor "red" notitle\
    ,r(t,1),(r(t,1) < 1) ? 1 : t  with lines lw 3 lt -1 notitle
