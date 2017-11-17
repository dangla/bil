# Gunplot file. To produce the figure type
#
# gnuplot PhaseDiagramOfCaO-SiO2-H2O.gp
#
#
# Options not affected by reset
#------------------------------


# To produce eps file
set term postscript eps enhanced color 20 size 10cm, 10cm
# To produce animation
#set term gif animate delay 100

set termoption solid
set termoption lw 2
set termoption font ",20"



# Linetypes
set linetype 1 lw 4 linecolor rgb "red"
set linetype 2 lw 4 linecolor rgb "light-magenta"
set linetype 3 lw 4 linecolor rgb "purple"
set linetype 4 lw 4 linecolor rgb "steelblue"
set linetype 5 lw 4 linecolor rgb "aquamarine"
set linetype 6 lw 4 linecolor rgb "bisque"
set linetype 7 lw 4 linecolor rgb "bisque"
set linetype 8 lw 4 linecolor rgb "light-goldenrod"
set linetype 9 lw 4 linecolor rgb "light-goldenrod"
#set linetype cycle 9




# Styles
set pointsize 2
set style line 10  lw 4 lt 1 pi -1 ps 2 pt 12 lc rgb "black"
set style line 11  lw 4 lt 1 pi -1 ps 2 pt 12 lc rgb "red"
set style line 12  lw 4 lt 1 pi -1 ps 2 pt 12 lc rgb "blue"
set style line 13  lw 4 lt 1 pi -1 ps 2 pt 12 lc rgb "green"
set style fill pattern 7 border lc rgb "black"






# Equilibrium constants
# ---------------------
logk_ch      = -5.14
logk_sh      = -2.76
logk_csh_Tob = -11.77
logk_csh_Jen = -17.36




# Functions
Langmuir(s,n) = (s)**n/(1 + (s)**n)
FITX1(s) = x1*Langmuir(s/s1,n1) + x2*Langmuir(s/s2,n2)
FITSSH1(s) = ((1 + (s/s1)**n1)**(-x1/n1)) * ((1 + (s/s2)**n2)**(-x2/n2))
RT = 2436
FITLOGK(s) = FITX1(s)*log(s) + log(FITSSH1(s))
FITG1(s) = RT*FITLOGK(s)

x1 = 0.88
x2 = 0.98

n1 = 0.88
s1 = 1.87e-6
n2 = 0.98
s2 = 6.9e-2

#fit FITX1(x) 'cal-greenberg-a1' us 2:1 via n1,s1,n2,s2
# After optimization
n1 = 0.88
s1 = 1.87e-6
n2 = 0.98
s2 = 6.9e-2

FITSSH_TOB(s) = 10**(logk_csh_Tob - logk_sh - 0.8333 * (logk_ch + log10(s)))
FITSSH_JEN(s) = 10**(logk_csh_Jen - logk_sh - 1.6666 * (logk_ch + log10(s)))
FITSSH_SH(s)  = 1




# Plot 1
# ------
# Size
set size square 1.,1.

# x-axis
set origin 0.,0.
#set xlabel '{/Symbol r}_{CO_2}/{/Symbol r}@_{CO_2}^{CH}' font ",24"
set xlabel 'Saturation index of portlandite' font ",24"
set logscale x
set format x "10^{%T}"
set xtics mirror 1.e-2
set xrange[1.e-8:1]

# y-axis
set ylabel 'Saturation index of silica' font ",24"
set logscale y
set format y "10^{%T}"
set yrange[1.e-6:1.e1]
set ytics mirror






# Legends
set key default
set key bottom left reverse Left spacing 1.5 samplen 2
#set key left at 1.e-5,8 reverse Left spacing 1.5 samplen 2
#set key left at 1.e-5,1.35 reverse Left spacing 1.5 samplen 1



set output 'PhaseDiagramOfCaO-SiO2-H2O.eps'


plot \
     FITSSH_SH(x)   w l lw 2 lt 4  title 'Tobermorite' \
    ,FITSSH_TOB(x)  w l lw 2 lt 2  title 'Silica' \
    ,FITSSH_JEN(x)  w l lw 2 lt 3  title 'Jennite' \
    ,FITSSH1(x)     w l lw 3 lt 1  title 'C-S-H' \


exit

# Plot 2
set output 'x-s_ch-greenberg.eps'
# y-axis
set ylabel 'C/S ratio' font ",24"
set nologscale y
set format y "%g"
set ytics 0.2
set ytics mirror
set yrange[0:1.8]
set key default
set key top left reverse Left spacing 1.5 samplen 2
plot \
     'cal-greenberg-a1' us 2:1 w p lw 4 pt 1 ps 2 title 'C/S (exp. Greenberg)' \
    , FITX1(x)  w l lw 6 lt 0 title 'fit' \
#    ,'cal-greenberg-a05' us 2:1 w p lw 4 pt 1 ps 2 title 'C/S (exp. Greenberg)' \
    
    
# Plot 3
set output 'ph-s_ch-greenberg.eps'
# y-axis
set ylabel 'pH' font ",20"
set nologscale y
set format y "%g"
set ytics 1
set ytics nomirror
set yrange[7:14]
set key default
set key top left reverse Left spacing 1.5 samplen 2
plot \
     'cal-greenberg-a1' us 2:17 w p lw 4 pt 1 ps 2 title 'pH (exp. Greenberg)'\
    ,'cal-greenberg-a1' us 2:18 w p lw 4 pt 1 ps 2 title 'pH (theory)' \
    ,'cal-greenberg-a05' us 2:17 w p lw 4 pt 1 ps 2 title 'pH (exp. Greenberg)'\
    ,'cal-greenberg-a05' us 2:18 w p lw 4 pt 1 ps 2 title 'pH (theory)' \
    
    
# Plot 4
set output 'c_i-s_ch-greenberg.eps'
# y-axis
set ylabel 'Concentration (mol/L)' font ",24"
set autoscale y
set nologscale y
set format y "%g"
set yrange[0:2.e-2]
set ytics nomirror 0,1.e-2
set key default
set key top left reverse Left spacing 1.5 samplen 2
plot \
     'cal-greenberg-a1' us 2:14 w p lw 4 pt 1 ps 2  title 'Charge+' \
    ,'cal-greenberg-a1' us 2:15 w p lw 4 pt 1 ps 2  title 'Charge-' \
    ,'cal-greenberg-a05' us 2:14 w p lw 4 pt 1 ps 2  title 'Charge+' \
    ,'cal-greenberg-a05' us 2:15 w p lw 4 pt 1 ps 2  title 'Charge-' \
#    ,'cal-greenberg-a1' us 2:19 w p lw 4 pt 1 ps 2  title 'Ionic Strength' \
#    ,'cal-greenberg-a05' us 2:19 w p lw 4 pt 1 ps 2  title 'Ionic Strength' \
    


# Plot 5
set output 'deltas_sht-s_ch-greenberg.eps'
# y-axis
set ylabel 'Saturation index of dissolved silica' font ",24"
set logscale y
set format y "10^{%T}"
set yrange[1.e-6:1.e1]
set ytics mirror 0.1
set key default
set key bottom left reverse Left spacing 1.5 samplen 2
set nokey
plot 1/0. notitle \
    , '+' using 1:(FITSSH1($1)*10**(-0.1*FITX1($1))):(FITSSH1($1)*10**(0.1*FITX1($1))) with filledcurves  ls 10 notitle \
    , FITSSH1(x)*10**(0.1*FITX1(x))   w l ls 10 notitle \
    , FITSSH1(x)*10**(-0.1*FITX1(x))  w l ls 10 notitle
    
    
# Plot 6: Gibbs energy
set output 'LippmannDiagram-greenberg.eps'
reset

set xlabel 'C/S fraction or S_{CH}' font ",24"
set xrange[0:1]
set logscale x
set format x "10^{%T}"
set xtics mirror 1.e-2
set xrange[1.e-10:1]
# y-axis
set ylabel 'Molar Gibbs Energy of CSH (kJ/mol)' font ",24"


set key bottom left reverse Left spacing 1.5 samplen 2

set sample 10000
plot '+' us ($1):(FITG1($1)*1.e-3) w l ls 10 lw 8 title 'solutus' \
    ,'+' us (FITX1($1)/1.8):(FITG1($1)*1.e-3) w l ls 11 lw 8 title 'solidus'


# Plot 7: Gibbs energy
set output 'Gibbs-greenberg.eps'
reset

set xlabel 'C/S ratio' font ",24"
set nologscale x
set xrange[0:2]
set xtics mirror 4.e-1
# y-axis
set ylabel 'Molar Gibbs Energy of CSH (kJ/mol)' font ",24"


set key bottom left reverse Left spacing 1.5 samplen 2

set parametric
set trange[-20:0]

plot FITX1(10**t),FITG1(10**t)*1.e-3 w l ls 10 lw 8 notitle



# Plot 8: Solubility product constant
set output 'LogK-x-greenberg.eps'

# x-axis
set origin 0.,0.
set xlabel 'C/S ratio' font ",24"
set nologscale x
set xrange[0:2]
set xtics mirror 4.e-1

# y-axis
#set ylabel 'Log(K) = log((a_{Ca^{2+}})^x * (a_{OH^{-}})^{2x} * (a_{H_4SiO_4}))' font ",24"
set ylabel 'Log(K) = log((a_{Ca^{2+}})^x * (a_{H^{+}})^{-2x} * (a_{H_4SiO_4}))' font ",24"


set key default
set key bottom left reverse Left spacing 1.5 samplen 2

set parametric
set trange[-20:0]

plot FITX1(10**t),(FITLOGK(10**t) - 5.2*FITX1(10**t) - 2.71 +28*FITX1(10**t)) w l ls 10 lw 8 notitle
