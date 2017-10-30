# Coordinates(1) ph(4) c_oh(5) c_h(6) c_ca(7) c_caoh(8) c_h2sio4(9) c_h3sio4(10) c_h4sio4(11) c_cah2sio4(12) c_cah3sio4(13) c_h2so4(14) c_hso4(15) c_so4(16) c_caso4aq(17) c_cahso4(18) c_k(19) c_koh(20) zn_ca_s(21) zn_si_s(22) s_ch(23) s_csh2(24) n_ch(25) n_csh2(26) n_csh(27) porosite(28) potentiel_electrique(29) charge(30) V_CSH(31) C/S(32) W_Si(33) W_Ca(34) W_S(35) P_CSH2(36) Damage(37) c_al(38) c_alo4h4(39) zn_al_s(40) s_ah3(41) s_afm(42) s_aft(43) s_c3ah6(44) n_ah3(45) n_afm(46) n_aft(47) n_c3ah6(48) W_Al(49) W_q(50) N_Ca(51) N_Si(52) N_S(53) N_Al(54) N_K(55) N_Cl(56) Saturation degree of crystal(57) Pore entry radius(58) Equilibrium saturation index of AFt(59) Crystallization pressure(60)



# Input file
file  = 'Guda2'

# Files
file1  = file.'.t1'
file2  = file.'.t2'
file3  = file.'.t3'
fileI = file2



# To produce eps file
set term postscript eps enhanced color 20

# To change the options of the terminal
set termoption solid
set termoption lw 2
set termoption font ",20"





# Figure 1
# --------
set output 'SolidContents.eps'


# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics 0.05 nomirror

# y-axis
set ylabel "Solid content (mol/L)" font ",24"
y0 = 0
y1 = 5
set yrange [y0:y1] noreverse nowriteback
set ytics 1 nomirror
set tics in
#set y2label "Gipsum content (mol/L)" font ",24"
#set y2range [y0:] noreverse nowriteback
#set y2tics 5


# Legends
set key top left reverse Left samplen 2 spacing 1.5 maxrows 5


# Styles
set style line  1 lt 1 lw 5 lc rgb "black"
set style line  2 lt 1 lw 5 lc rgb "red"
set style line  3 lt 1 lw 5 lc rgb "orange"
set style line  4 lt 1 lw 5 lc rgb "blue"
set style line  5 lt 1 lw 5 lc rgb "web-green

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


# Size
set size square 0.75,1.
set origin 0,0

plot \
     fileI us ($1):($25) w l ls 1 title 'n_{CH}' \
    ,fileI us ($1):($26) w l ls 2  title 'n_{CSH_2}' \
    ,fileI us ($1):($27) w l ls 3 title 'n_{CSH}' \
    ,fileI us ($1):($47*1000) w l ls 4 title 'n_{AFt} 10^{3}' \
    ,fileI us ($1):($32) w l ls 5 title 'C/S'



reset


# Figure 2
set output 'SaturationIndexes.eps'


# Legends
set key center left reverse Left samplen 1 spacing 1.5 maxrows 5


# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 0.2
set xrange [x0:x1] noreverse nowriteback
set xtics 0.05 nomirror
set autoscale x

# y-axis
set ylabel "Saturation index" font ",24"
y0 = 0
y1 = 2.1
set yrange [y0:y1] noreverse nowriteback
set ytics nomirror
set autoscale y
set tics in



plot \
     fileI  us ($1):($23) w l lt 5 title '{/Symbol b}_{CH}' \
    ,fileI  us ($1):($24) w l lt 6 title '{/Symbol b}_{CSH2}' \
    ,fileI  us ($1):($41) w l lt 1 title '{/Symbol b}_{AH3}' \
    ,fileI  us ($1):($42) w l lt 2 title '{/Symbol b}_{AFm}' \
    ,fileI  us ($1):($43) w lp lt 3 title '{/Symbol b}_{AFt}' \
    ,fileI  us ($1):($44) w l lt 4 title '{/Symbol b}_{C3AH6}' \



reset


# Figure 3
set output 'IonConcentrations.eps'


# Legends
set key center left reverse Left samplen 1 spacing 1.5 maxrows 5


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
     fileI  us ($1):($16) w l lt 1 title 'log({/Symbol r}_{SO@_4^{2-}})' \
    ,fileI  us ($1):($15) w l lt 2 title 'log({/Symbol r}_{HSO@_4^{-}})' \
    ,fileI  us ($1):($14) w l lt 3 title 'log({/Symbol r}_{H_2SO_4})' \



