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





#   x	  r	 sc	 beta	 BETAi	 stress	 s100
#Experimental and numerical study on cement paste degradation under ESA
# strain ratio a. 3mm-1.5g/L ; b.3mm-30g/L ; 
# x1a(1)    stress(MPa)(2) strain(3)  x2b(4) stress(MPa)(5) strain(6)


# Time(1) ph(2) c_ca(3) c_caoh(4) c_h2sio4(5) c_h3sio4(6) c_h4sio4(7) c_cah2sio4(8) c_cah3sio4(9) c_h2so4(10) c_hso4(11) c_so4(12) c_caso4aq(13) c_cahso4(14) c_al(15) c_alo4h4(16) c_k(17) c_koh(18) C_Ca(19) C_Si(20) C_S(21) C_Al(22) C_K(23) s_ch(24) s_csh2(25) s_ah3(26) s_afm(27) s_aft(28) s_c3ah6(29) n_ch(30) n_csh2(31) n_csh(32) n_ah3(33) n_afm(34) n_aft(35) n_c3ah6(36) porosite(37) potentiel_electrique(38) charge(39) V_CSH(40) C/S(41) W_Si(42) W_Ca(45) W_S(48) W_Al(51) W_q(54) P_CSH2(57) Damage(58) N_Ca(59) N_Si(60) N_S(61) N_Al(62) N_K(63) N_Cl(64) Saturation degree of crystal(65) Interface equilibrium saturation index of AFt(66) Pore wall equilibrium saturation index of AFt(67) Crystallization pressure(68) Strain tensor(69) Stress tensor(78) Displacement vector(87) Damage(90) Hardening variable(91) Yield function(92)

# Time(1) ph(2) c_ca(3) c_caoh(4) c_h2sio4(5) c_h3sio4(6) c_h4sio4(7) c_cah2sio4(8) c_cah3sio4(9) c_h2so4(10) c_hso4(11) c_so4(12) c_caso4aq(13) c_cahso4(14) c_al(15) c_alo4h4(16) c_k(17) c_koh(18) C_Ca(19) C_Si(20) C_S(21) C_Al(22) C_K(23) s_ch(24) s_csh2(25) s_ah3(26) s_afm(27) s_aft(28) s_c3ah6(29) n_ch(30) n_csh2(31) n_csh(32) n_ah3(33) n_afm(34) n_aft(35) n_c3ah6(36) n_so4^ads(37) porosite(38) potentiel_electrique(39) charge(40) V_CSH(41) C/S(42) W_Si(43) W_Ca(46) W_S(49) W_Al(52) W_q(55) P_CSH2(58) Damage(59) N_Ca(60) N_Si(61) N_S(62) N_Al(63) N_K(64) N_Cl(65) Saturation degree of crystal(66) Interface equilibrium saturation index of AFt(67) Pore wall equilibrium saturation index of AFt(68) Crystallization pressure(69) Strain tensor(70) Stress tensor(79) Displacement vector(88) Damage(91) Hardening variable(92) Yield function(93)



# Linetypes
set linetype 1  lw 3 linecolor rgb "red" 
set linetype 2  lw 3 linecolor rgb "green" 
set linetype 3  lw 3 linecolor rgb "blue"
set linetype 4  lw 3 linecolor rgb "purple"
set linetype 5  lw 3 linecolor rgb "cyan"
set linetype 6  lw 3 linecolor rgb "golden"
set linetype 7  lw 3 linecolor rgb "dark-violet"
set linetype 8  lw 3 linecolor rgb "light-goldenrod"
set linetype 9  lw 3 linecolor rgb "sea-green"
set linetype 10 lw 3 linecolor rgb "black"
set linetype 11 lw 3 linecolor rgb "dark-red"
set linetype 12 lw 3 linecolor rgb "steelblue"
set linetype 13 lw 3 linecolor rgb "dark-orange"
set linetype 14 lw 3 linecolor rgb "aquamarine"
set linetype 15 lw 3 linecolor rgb "goldenrod"
set linetype 16 lw 3 linecolor rgb "light-magenta"
set linetype 17 lw 3 linecolor rgb "dark-blue"
set linetype 18 lw 3 linecolor rgb "dark-green"
#set linetype 19 lw 3 linecolor rgb "dark-purple"
set linetype 20 lw 3 linecolor rgb "light-green"
set linetype 21 lw 3 linecolor rgb "light-red"
set linetype cycle 15



# Input files
basefile  = 'SulfacoESA3d'
file1  = 'Ma.txt'
file2  = basefile.'.p1'
#file2  = basefile.'30.p1'
file3  = basefile.'1.5.p1'
file4  = 'strain30.txt'
file5  = 'strain1.5.txt'
file6  = 'stress-strain30.sample'
file7  = 'stress-strain1.5.sample'




# To produce eps file
set term postscript eps enhanced color 20

# To change the options of the terminal
set termoption solid
set termoption lw 3
set termoption font ",20"





# Figure 1
# --------
set output basefile.'Stress.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 230
set xrange [x0:x1] noreverse writeback
set xtics  nomirror
set autoscale x
#set logscale x


# y-axis
set ylabel "Compressive axial stress (MPa)" font ",24"
y0 = 0
y1 = 1
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y


set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

# Rod stiffness
k = 9.257e9
# Poisson's ratio
nu = 0.2
# Young's modulus
E = 36.e9

# coef
coef = (1-2*nu)/(1 + E/k)
# unit
MPa = 1.e-6

plot \
     file2  us ($1/86400):(coef*$65*$68*MPa)  w l lt 1 title 'Simulation 3mm--30g/L' \
    ,file2  us ($1/86400):(-coef*$78*MPa)  w l lt 2 title 'Simulation 3mm--30g/L'  

reset
exit



# Figure 2
# --------
set output 'mastrain.eps'


# x-axis
set xlabel "Time (day)" font ",24"
x0 = 0
x1 = 230
set xrange [x0:x1] noreverse writeback
set xtics  nomirror
set autoscale x
#set logscale x


# y-axis
set ylabel "Expansion (%)" font ",24"
y0 = 0
y1 = 1
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y


set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     file5 us ($1):($2)  w l lt 3 title '3mm--1.5g/L'  \
     ,file4 us ($1):($2)  w l lt 1 title '3mm--30g/L'  

reset
exit



# Figure 3
# --------
set output 'Ma-stress-strain.eps'


# x-axis
set xlabel "Axial strain (10^{-4})" font ",24"
x0 = 0
x1 = 0.015
set xrange [x0:x1] noreverse writeback
set xtics  nomirror
set autoscale x
#set logscale x


# y-axis
set ylabel "Compressive axial stress (MPa)" font ",24"
y0 = 0
y1 = 1.3
set yrange [y0:y1] noreverse nowriteback
set ytics  nomirror
set tics in
set autoscale y


set grid


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5




# Size
set size square 0.75,1.
set origin 0,0

plot \
     file7 us ($4*100):($2)  w l lt 3 title '3mm--1.5g/L'  \
     ,file6 us ($4*100):($2)  w l lt 1 title '3mm--30g/L'  

exit




