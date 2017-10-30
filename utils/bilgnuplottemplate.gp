# Units
# -----
m   = 1
m2  = (m*m)
m3  = (m*m*m)
cm  = (1.e-2*m)
cm3 = (cm*cm*cm)
nm  = (1.e-9*m)
mm  = (1.e-3*m)
Pa  = 1.
MPa = 1.e6
GPa = 1.e9
J   = Pa*m3
mol = 1
K   = 1
sec = 1
day = 86400.*sec


# Physical constants
# ------------------
R   = 8.314 * J / mol / K
T   = 293. * K
RT  = R * T * J / mol
V_w = 18. * cm3 / mol
gamma = 0.072 * J / m2



# To produce eps file
set term postscript eps enhanced color 20

# To change the options of the terminal
set termoption solid
set termoption lw 2
set termoption font ",20"


# Linetypes
set linetype 1 lw 3 linecolor rgb "red"
set linetype 2 lw 3 linecolor rgb "light-magenta"
set linetype 3 lw 3 linecolor rgb "purple"
set linetype 4 lw 3 linecolor rgb "steelblue"
set linetype 5 lw 3 linecolor rgb "aquamarine"
set linetype 6 lw 3 linecolor rgb "bisque"
set linetype 7 lw 3 linecolor rgb "bisque"
set linetype 8 lw 3 linecolor rgb "light-goldenrod"
set linetype 9 lw 3 linecolor rgb "light-goldenrod"
set linetype 10 lw 3 linecolor rgb "black"
set linetype 11 lw 3 linecolor rgb "black"





# Figure 1
# --------
set output 'Fig1.eps'

# Size
set size square 0.75,1.
set origin 0,0


# x-axis
set xlabel "x-axis {/Symbol b}@^{eq}_p" font ",24"
x0 = 1.e0
x1 = 1.e6
set xrange [x0:x1] noreverse nowriteback
set xtics out nomirror
set format x '10^{%T}'
#set autoscale x
set logscale x


# y-axis
set ylabel "y-axis (unit)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
#set ytics 0.2 nomirror
#set autoscale y
set ytics in nomirror


# Legends
set key top left reverse Left samplen 1 spacing 1.5 maxrows 5



b  = 1
V_c = 710 * cm3
a_p = 1.e-8 * mol / m3 / sec
dt = 1 * day
beta = 1.e5
E = 20 * GPa
strain0 = 1.e-4
strainf = 3.9e-3

Pc(beta_p) = RT / V_c * log(beta_p)
stress(S_c,beta_p) = b * S_c * Pc(beta_p)
strain(S_c,beta,beta_p,t) = (t <= 0) ? 0 : strain(S_c,beta,beta_p,t-dt) + dt * V_c / S_c * a_p * (1 - beta_p/beta)

Estrain(S_c,beta,beta_p,t) = E * strain(S_c,beta,beta_p,t)

plot \
     S_c = 0.1 \
    ,stress(S_c,x)/MPa w l lt 1 title sprintf("bS_cP_c for S_c = %g",S_c) \
    ,t = 5*dt \
    ,Estrain(S_c,beta,x,t)/MPa w l lt 2 title sprintf("E{/Symbol e} at t = %g days",t/day) \
    ,t = 16*dt \
    ,Estrain(S_c,beta,x,t)/MPa w l lt 3 title sprintf("E{/Symbol e} at t = %g days",t/day) \
    ,t = 100*dt \
    ,Estrain(S_c,beta,x,t)/MPa w l lt 4 title sprintf("E{/Symbol e} at t = %g days",t/day) \


reset


# Figure 2
set output 'Fig2.eps'

# Size
set size square 0.75,1.
set origin 0,0



# Input file
file  = 'anyfile'

# Files
file1  = file.'.t1'
file2  = file.'.t2'
file3  = file.'.t3'
fileI = file2



# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics 0.05 nomirror
set autoscale x

# y-axis
set ylabel "Ion concentration (mol/L)" font ",24"
y0 = 0
y1 = 2.1
#set yrange [y0:y1] noreverse nowriteback
set ytics nomirror
set autoscale y
set tics in


plot \
     fileI us ($1):($25) w l ls 1 title 'n_{CH}' \
    ,fileI us ($1):($26) w l ls 2  title 'n_{CSH_2}' \
    ,fileI us ($1):($27) w l ls 3 title 'n_{CSH}' \
    ,fileI us ($1):($47*1000) w l ls 4 title 'n_{AFt} 10^{3}' \
    ,fileI us ($1):($32) w l ls 5 title 'C/S'
    

reset

exit


     
