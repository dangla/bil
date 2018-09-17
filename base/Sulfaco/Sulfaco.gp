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





# Time(1) ph(2) c_oh(3) c_h(4) c_ca(5) c_caoh(6) c_h2sio4(7) c_h3sio4(8) c_h4sio4(9) c_cah2sio4(10) c_cah3sio4(11) c_h2so4(12) c_hso4(13) c_so4(14) c_caso4aq(15) c_cahso4(16) c_k(17) c_koh(18) zn_ca_s(19) zn_si_s(20) s_ch(21) s_csh2(22) n_ch(23) n_csh2(24) n_csh(25) porosite(26) potentiel_electrique(27) charge(28) V_CSH(29) C/S(30) W_Si(31) W_Ca(32) W_S(33) P_CSH2(34) Damage(35) c_al(36) c_alo4h4(37) zn_al_s(38) s_ah3(39) s_afm(40) s_aft(41) s_c3ah6(42) n_ah3(43) n_afm(44) n_aft(45) n_c3ah6(46) W_Al(47) W_q(48) N_Ca(49) N_Si(50) N_S(51) N_Al(52) N_K(53) N_Cl(54) Saturation degree of crystal(55) Pore entry radius(56) Equilibrium saturation index of AFt(57) Crystallization pressure(58) Strain(59)





# Linetypes
set linetype 1  lw 2 linecolor rgb "red"
set linetype 2  lw 2 linecolor rgb "light-magenta"
set linetype 3  lw 2 linecolor rgb "purple"
set linetype 4  lw 2 linecolor rgb "steelblue"
set linetype 5  lw 2 linecolor rgb "aquamarine"
set linetype 6  lw 2 linecolor rgb "bisque"
set linetype 7  lw 2 linecolor rgb "dark-violet"
set linetype 8  lw 2 linecolor rgb "light-goldenrod"
set linetype 9  lw 2 linecolor rgb "sea-green"
set linetype 10 lw 2 linecolor rgb "cyan"
set linetype 11 lw 2 linecolor rgb "dark-red"
set linetype 12 lw 2 linecolor rgb "blue"
set linetype 13 lw 2 linecolor rgb "dark-orange"
set linetype 14 lw 2 linecolor rgb "black"
set linetype 15 lw 2 linecolor rgb "goldenrod"
set linetype cycle 15



# Input files
file1  = 'Sulfaco'


# Files
file1p1  = file1.'.p1'
fileI    = file1p1



# To produce eps file
set term postscript eps enhanced color 20

# To change the options of the terminal
set termoption solid
set termoption lw 2
set termoption font ",20"





# Figure 1
# --------
set output 'Strain.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Strain (-)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics 0.001 nomirror
set tics in
set autoscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):($59) w l lt 1 title '1'  \


reset






# Figure 2
# --------
set output 'EttringiteSaturationIndex.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
set logscale x

# y-axis
set ylabel "Saturation index of ettringite (-)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
set logscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI  us ($1/day):($41) w l lt 1 title '1' \

reset







# Figure 3
# --------
set output 'EttringiteContent.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
set logscale x

# y-axis
set ylabel "Ettringite content (mol/L)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI  us ($1/day):($45) w l lt 1 title '1' \

reset







# Figure 4
# --------
set output 'CrystallizationPressure.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Crystallization pressure (Pa)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):($58) w l lt 1 title '1' \

reset







# Figure 5
# --------
set output 'pH.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "pH (-)" font ",24"
y0 = 0
y1 = 14
set yrange [y0:y1] noreverse nowriteback
set ytics 0.05 nomirror
set tics in
set autoscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):($2) w l lt 1 title '1' \

reset






# Figure 6
# --------
set output 'SaturationIndexes.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Saturation indexes (-)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
set logscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key center left reverse Left samplen 2 spacing 1.5 maxrows 2




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):($41/$57) w l lt 1 title '1: supersat AFt' \
    ,fileI   us ($1/day):($22)     w l lt 2 title '1: CSH_2' \



reset







# Figure 7
# --------
set output 'SolidContents.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Solid content (mol/L)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key center left reverse Left samplen 2 spacing 1.5 maxrows 3




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):($23) w l lt 7 title '1: CH' \
    ,fileI   us ($1/day):($45) w l lt 1 title '1: AFt' \
    ,fileI   us ($1/day):($46) w l lt 3 title '1: C_3AH_6' \


reset






# Figure 8
# --------
set output 'IonConcentrations.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Ion concentration (mol/L)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y
#set logscale y
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5

set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI   us ($1/day):(10**$14) w l lt 1 title '1: SO@_4^{-2}' \

reset




# Figure 9
# --------
set output 'StressStrain.eps'


# x-axis
set xlabel "Strain (-)" font ",24"
x0 = 0
x1 = 1.e-2
set xrange [x0:x1] noreverse nowriteback
set format x '%g'
set xtics 5.e-3 out nomirror
set autoscale x


# y-axis
set ylabel "Effective stress: S_cP_c (Pa)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
#set ytics 0.2 nomirror
set autoscale y
set ytics in nomirror


# Legends
set key top left reverse Left samplen 1 spacing 1.5 maxrows 5

# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI  us ($59):($55*$58) w l lt 1 title '1' \




reset




# Figure 9
# --------
set output 'Charge.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics  nomirror
set autoscale x
#set logscale x


# y-axis
set ylabel "Charge (M)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
#set ytics 0.2 nomirror
set autoscale y
set ytics in nomirror


# Legends
set key top left reverse Left samplen 1 spacing 1.5 maxrows 5

# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI  us ($1/day):($28) w l lt 1 title '1' \




reset




