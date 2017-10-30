# Options not affected by reset
#------------------------------
`set term`       , `set output`  , `set loadpath`, 
`set fontpath`   , `set linetype`, `set encoding`, 
`set decimalsign`, `set locale`  , `set psdir`

# To produce eps file
set term postscript eps enhanced color 20
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
#set pointsize 2
set style line 1  lw 4 lt 1 pi -10 ps 2 pt 4 lc rgb "red"
set style line 2  lw 4 lt 2 pi -10 ps 2 pt 8 lc rgb "light-magenta"
set style line 3  lw 4 lt 3 pi -10 ps 2 pt 6 lc rgb "purple"
set style line 4  lw 4 lt 4 pi -10 ps 2 pt 4 lc rgb "steelblue"
set style line 5  lw 4 lt 5 pi -10 ps 2 pt 8 lc rgb "bisque"
set style line 6  lw 4 lt 6 pi -10 ps 2 pt 6 lc rgb "green"
set style line 6  lw 4 lt 7 pi -10 ps 2 pt 6 lc rgb "green"
set style line 7  lw 4 lt 8 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"
set style line 7  lw 4 lt 9 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"



# Main file
mainfile = 'Yuan3c_si'

# Coordinates(1) ph(4) c_oh(5) c_h(6) c_ca(7) c_caoh(8) c_h2sio4(9) c_h3sio4(10) c_h4sio4(11) c_cah2sio4(12) c_cah3sio4(13) c_h2so4(14) c_hso4(15) c_so4(16) c_caso4aq(17) c_cahso4(18) c_k(19) c_koh(20) zn_ca_s(21) zn_si_s(22) s_ch(23) s_csh2(24) n_ch(25) n_csh2(26) n_si_s(27) porosite(28) potentiel_electrique(29) charge(30) V_CSH(31) C/S(32) W_Si(33) W_Ca(34) W_S(35) P_CSH2(36) Damage(37) c_al(38) c_alo4h4(39) zn_al_s(40) s_ah3(41) s_afm(42) s_aft(43) s_c3ah6(44) n_ah3(45) n_afm(46) n_aft(47) n_c3ah6(48) W_Al(49) W_q(50) N_Ca(51) N_Si(52) N_S(53) N_Al(54) N_K(55) N_Cl(56)

# Time(1) ph(2) c_oh(3) c_h(4) c_ca(5) c_caoh(6) c_h2sio4(7) c_h3sio4(8) c_h4sio4(9) c_cah2sio4(10) c_cah3sio4(11) c_h2so4(12) c_hso4(13) c_so4(14) c_caso4aq(15) c_cahso4(16) c_k(17) c_koh(18) zn_ca_s(19) zn_si_s(20) s_ch(21) s_csh2(22) n_ch(23) n_csh2(24) n_si_s(25) porosite(26) potentiel_electrique(27) charge(28) V_CSH(29) C/S(30) W_Si(31) W_Ca(32) W_S(33) P_CSH2(34) Damage(35) c_al(36) c_alo4h4(37) zn_al_s(38) s_ah3(39) s_afm(40) s_aft(41) s_c3ah6(42) n_ah3(43) n_afm(44) n_aft(45) n_c3ah6(46) W_Al(47) W_q(48) N_Ca(49) N_Si(50) N_S(51) N_Al(52) N_K(53) N_Cl(54)


# Plot 1
#-------
set output 'SolidContents.eps'
# Size
set size square 1,1


# x-axis
set xlabel 'Width (dm)' font ",24" offset 2,-1 #rotate parallel
x0 = -4
x1 = 0
#set xrange [x0:x1] noreverse nowriteback
#set xtics nomirror out norotate x0,2,x1


# y-axis
set ylabel 'Solid content (mol/L)' font ",24" offset -1.5,-1.5 #rotate parallel
y0 = -8
y1 = +8
#set yrange[y0:y1] noreverse nowriteback
#set ytics  nomirror out scale 4. norotate offset 0,-0.5 y0,4,y1
#set logscale y
set ytics  nomirror
set y2tics  nomirror
set y2label 'Gipsum contents (mol/L)' font ",24" offset -1.5,-1.5 


# Legends
#set key lmargin  horizontal Left reverse samplen 2 spacing 1
#set key lmargin  horizontal Right samplen 2 spacing 1
set key top left

# Data files
file = mainfile.'.t2'

plot \
      file us 1:25 w l ls 1  title 'n_{ch}'\
    , file us 1:26 axis x1y2 w l ls 2  title 'n_{csh2}'\
    , file us 1:27 w l ls 3  title 'n_{si}(s)'\
    , file us 1:45 w l ls 4  title 'n_{ah3}'\
    , file us 1:46 w l ls 5  title 'n_{afm}'\
    , file us 1:47 w l ls 6  title 'n_{aft}'\
    , file us 1:48 w l ls 7  title 'n_{c3ah6}'\



# Plot 2
#-------
set output 'SaturationDegreeOfDissolvedSolids.eps'
# Size
set size square 1,1


# x-axis
set xlabel 'Width (dm)' font ",24" offset 2,-1 #rotate parallel
x0 = 0.98
x1 = 1.
#set xrange [x0:x1] noreverse nowriteback
#set xtics nomirror out norotate x0,2,x1
set autoscale x


# y-axis
set ylabel 'Saturation degree of dissolved solid (-)' font ",24" offset -1.5,-1.5 #rotate parallel
y0 = 0
y1 = 8
#set yrange[y0:y1] noreverse nowriteback
#set ytics  nomirror out scale 4. norotate offset 0,-0.5 y0,4,y1
set ytics  nomirror
unset y2tics
unset y2label #'Saturation degree of dissolved gipsum (-)' font ",24" offset -1.5,-1.5 #rotate parallel


# Legends
#set key lmargin  horizontal Left reverse samplen 2 spacing 1
#set key lmargin  horizontal Right samplen 2 spacing 1
set key top left

# Data files
file = mainfile.'.t2'


plot \
      file us 1:23 w l ls 1 title 's_{ch}'\
    , file us 1:24 w l ls 2 title 's_{csh2}'\
    , file us 1:41 w l ls 4 title 's_{ah3}'\
    , file us 1:42 w l ls 5 title 's_{afm}'\
    , file us 1:43 w l ls 6 title 's_{aft}'\
    , file us 1:44 w l ls 7 title 's_{c3ah6}'\

