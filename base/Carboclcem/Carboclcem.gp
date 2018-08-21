# Units
# -----
m   = 1
m2  = (m*m)
m3  = (m*m*m)
cm  = (1.e-2*m)
cm3 = (cm*cm*cm)
Liter = 1.e3*cm3
nm  = (1.e-9*m)
mm  = (1.e-3*m)
Pa  = 1.
MPa = 1.e6
GPa = 1.e9
J   = Pa*m3
mol = 1
K   = 1
sec = 1
hour = 3600 * sec
day = 86400 * sec
Mol = mol / Liter


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





# Input file
file  = 'Carboclcem'

# Files
filet0  = file.'.t1'
filet1  = file.'.t1'



# Outputs at a given time
#------------------------
# Coordinates(1) c_co2(4) ph(5) n_Si_s(6) porosity(7) n_CH(8) c_ca(9) c_co3(10) c_hco3(11) n_CC(12) c_h(13) c_oh(14) saturation(15) grad_psi(16) charge(17) c_na(18) c_naoh(19) c_nahco3(20) c_naco3(21) c_k(22) c_koh(23) c_caoh(24) c_cahco3(25) c_caco3aq(26) c_caoh2aq(27) p_l(28) c_h3sio4(29) n_Na(30) n_Ca(31) n_Si(32) n_Ca_s(33) c_cah2sio4(34) c_cah3sio4(35) CsurS(36) c_h4sio4(37) c_h2sio4(38) c_cl(39) I(40) x_csh(41) n_si_s(42) s_ch(43) s_sh(44) k_l(45) verma-pruess(46) dn_chsdt(47) dn1sdt(48) coeff_dnCH(49) v_csh(50) v_ch(51) v_cc(52)


# Outputs at a given point
# ------------------------
# Time(1) c_co2(2) ph(3) n_Si_s(4) porosity(5) n_CH(6) c_ca(7) c_co3(8) c_hco3(9) n_CC(10) c_h(11) c_oh(12) saturation(13) grad_psi(14) charge(15) c_na(16) c_naoh(17) c_nahco3(18) c_naco3(19) c_k(20) c_koh(21) c_caoh(22) c_cahco3(23) c_caco3aq(24) c_caoh2aq(25) p_l(26) c_h3sio4(27) n_Na(28) n_Ca(29) n_Si(30) n_Ca_s(31) c_cah2sio4(32) c_cah3sio4(33) CsurS(34) c_h4sio4(35) c_h2sio4(36) c_cl(37) I(38) x_csh(39) n_si_s(40) s_ch(41) s_sh(42) k_l(43) verma-pruess(44) dn_chsdt(45) dn1sdt(46) coeff_dnCH(47) v_csh(48) v_ch(49) v_cc(50)




# Figure 1
# --------
set output 'Concentrations.eps'

# Size
set size square 2.25,1.
set origin 0,0



# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 5
set xrange [x0:x1] noreverse nowriteback
#set xtics 2 nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Concentrations (mol/L)" font ",24"
y0 = 0
y1 = 1
set yrange [y0:y1] noreverse nowriteback
set ytics out 0.1 nomirror
set autoscale y
set y2range [0:12] noreverse nowriteback
set y2tics out 1
set autoscale y2
set tics out


# grid
set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5



set size square 0.75,1.
set origin 0,0

plot \
     filet0 us ($1):($4)  w l lt 1 title 'CO2' \
    ,filet0 us ($1):($9)  w l lt 2 title 'Ca' \
    ,filet0 us ($1):($18) w l lt 3 title 'Na' \
    ,filet0 us ($1):($39) w l lt 4 title 'Cl' \
    ,filet0 us ($1):($5) axis x1y2 w l lt 0 title 'pH' \

    

reset





# Figure 2
# --------
set output 'SolidContents.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 24
set xrange [x0:x1] noreverse nowriteback
set xtics 0.2 nomirror
set autoscale x
#set logscale x

# y-axis
set ylabel "Solid contents (mol/L)" font ",24"
y0 = -20
y1 = 0
set yrange [y0:y1] noreverse nowriteback
#set ytics out 5 nomirror
set autoscale y
set tics out


# grid
set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5
     

set size square 0.75,1.
set origin 0,0

plot \
     filet0 us ($1):($8)     w l lt 1 title 'CH' \
    ,filet0 us ($1):($12)    w l lt 2 title 'CC' \

    

    

reset





# Figure 3
# --------
set output 'LiquidSaturation.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Depth (dm)" font ",24"
x0 = 0
x1 = 24
set xrange [x0:x1] noreverse nowriteback
set xtics 0.2 nomirror
set autoscale x

# y-axis
set ylabel "Liquid saturation degree (-)" font ",24"
y0 = -20
y1 = 0
set yrange [y0:y1] noreverse nowriteback
#set ytics out 5 nomirror
set autoscale y
set tics out


# grid
set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5
     

set size square 0.75,1.
set origin 0,0

plot \
     filet0 us ($1):($15)    w l lt 1 title 'S_L'

    
    

reset
exit





# Figure 3
# --------
set output 'SaltConcentration.eps'

# Size
set size square 2.25,1.
set origin 0,0



# x-axis
set xlabel "Depth (cm)" font ",24"
x0 = 0
x1 = 5
set xrange [x0:x1] noreverse nowriteback
set xtics 1 nomirror
#set autoscale x

# y-axis
set ylabel "Salt concentration (mol/L)" font ",24"
y0 = 0
y1 = 3
#set yrange [y0:y1] noreverse nowriteback
set ytics out 1 nomirror
set autoscale y
set tics out


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5

set multiplot


set size square 0.75,1.
set origin 0,0

set title '0 M' font ",24"



plot \
     file0t1 us ($1/cm):($5/Mol) w l lt 1 title '6 hrs' \
    ,file0t4 us ($1/cm):($5/Mol) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 0.75,0

set title '1 M' font ",24"

plot \
     file1t1 us ($1/cm):($5/Mol) w l lt 1 title '6 hrs' \
    ,file1t4 us ($1/cm):($5/Mol) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 1.5,0

set title '3 M' font ",24"


plot \
     file3t1 us ($1/cm):($5/Mol) w l lt 1 title '6 hrs' \
    ,file3t4 us ($1/cm):($5/Mol) w l lt 4 title sprintf("%3.0f days",t4/day) \

unset multiplot

    

reset





# Figure 4
# --------
set output 'Strain.eps'

# Size
set size square 2.25,1.
set origin 0,0



# x-axis
set xlabel "Depth (cm)" font ",24"
x0 = 0
x1 = 5
set xrange [x0:x1] noreverse nowriteback
set xtics 1 nomirror
#set autoscale x

# y-axis
set ylabel "Strain" font ",24"
y0 = 0
y1 = 3
set yrange [y0:y1] noreverse nowriteback
set ytics out 1.e-4 nomirror
set autoscale y
set tics out


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5

set multiplot


set size square 0.75,1.
set origin 0,0

set title '0 M' font ",24"



plot \
     file0t1 us ($1/cm):($17) w l lt 1 title '6 hrs' \
    ,file0t4 us ($1/cm):($17) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 0.75,0

set title '1 M' font ",24"



plot \
     file1t1 us ($1/cm):($17) w l lt 1 title '6 hrs' \
    ,file1t4 us ($1/cm):($17) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 1.5,0

set title '3 M' font ",24"



plot \
     file3t1 us ($1/cm):($17) w l lt 1 title '6 hrs' \
    ,file3t4 us ($1/cm):($17) w l lt 4 title sprintf("%3.0f days",t4/day) \

unset multiplot

reset





# Figure 5
# --------
set output 'IcePressure.eps'

# Size
set size square 2.25,1.
set origin 0,0



# x-axis
set xlabel "Depth (cm)" font ",24"
x0 = 0
x1 = 5
set xrange [x0:x1] noreverse nowriteback
set xtics 1 nomirror
#set autoscale x

# y-axis
set ylabel "Ice pressure (MPa)" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 10 nomirror
set autoscale y
set tics out


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5

set multiplot


set size square 0.75,1.
set origin 0,0

set title '0 M' font ",24"



plot \
     file0t1 us ($1/cm):($14/MPa) w l lt 1 title '6 hrs' \
    ,file0t4 us ($1/cm):($14/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 0.75,0

set title '1 M' font ",24"



plot \
     file1t1 us ($1/cm):($14/MPa) w l lt 1 title '6 hrs' \
    ,file1t4 us ($1/cm):($14/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 1.5,0

set title '3 M' font ",24"



plot \
     file3t1 us ($1/cm):($14/MPa) w l lt 1 title '6 hrs' \
    ,file3t4 us ($1/cm):($14/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \
    
unset multiplot
    

reset





# Figure 6
# --------
set output 'LiquidPressure.eps'

# Size
set size square 2.25,1.
set origin 0,0



# x-axis
set xlabel "Depth (cm)" font ",24"
x0 = 0
x1 = 5
set xrange [x0:x1] noreverse nowriteback
set xtics 1 nomirror
#set autoscale x

# y-axis
set ylabel "Liquid pressure (MPa)" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 10 nomirror
set autoscale y
set tics out


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5

set multiplot


set size square 0.75,1.
set origin 0,0

set title '0 M' font ",24"



plot \
     file0t1 us ($1/cm):($4/MPa) w l lt 1 title '6 hrs' \
    ,file0t4 us ($1/cm):($4/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 0.75,0

set title '1 M' font ",24"



plot \
     file1t1 us ($1/cm):($4/MPa) w l lt 1 title '6 hrs' \
    ,file1t4 us ($1/cm):($4/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \


set size square 0.75,1.
set origin 1.5,0

set title '3 M' font ",24"



plot \
     file3t1 us ($1/cm):($4/MPa) w l lt 1 title '6 hrs' \
    ,file3t4 us ($1/cm):($4/MPa) w l lt 4 title sprintf("%3.0f days",t4/day) \
    
unset multiplot
    

reset





# Figure 7
# --------
set output 'LiquidPressureVsTime.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Time (day)" font ",24"
x0 = 1.e-2
x1 = 100
set xrange [x0:x1] noreverse nowriteback
set xtics out 10 nomirror
#set autoscale x
set logscale x

# y-axis
set ylabel "Liquid pressure (MPa)" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 10 nomirror
set autoscale y
set tics out

set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5


plot \
     file0p1 us ($1/day):($2/MPa) w l lt 1 title '0 M' \
    ,file1p1 us ($1/day):($2/MPa) w l lt 2 title '1 M' \
    ,file3p1 us ($1/day):($2/MPa) w l lt 3 title '3 M' \

    

reset





# Figure 8
# --------
set output 'SaltConcentrationVsTime.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Time (day)" font ",24"
x0 = 1.e-2
x1 = 100
set xrange [x0:x1] noreverse nowriteback
set xtics out 10 nomirror
#set autoscale x
set logscale x

# y-axis
set ylabel "Salt concentration (mol/L)" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 1 nomirror
set autoscale y
set tics out


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5


plot \
     file0p1 us ($1/day):($3/Mol) w l lt 1 title '0 M' \
    ,file1p1 us ($1/day):($3/Mol) w l lt 2 title '1 M' \
    ,file3p1 us ($1/day):($3/Mol) w l lt 3 title '3 M' \

    

reset





# Figure 9
# --------
set output 'StrainVsTime.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Time (day)" font ",24"
x0 = 1.e-2
x1 = 100
set xrange [x0:x1] noreverse nowriteback
set xtics out 10 nomirror
#set autoscale x
set logscale x

# y-axis
set ylabel "Strain" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 1.e-4 nomirror
set autoscale y
set tics out

set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5


plot \
     file0p1 us ($1/day):($15) w l lt 1 title '0 M' \
    ,file1p1 us ($1/day):($15) w l lt 2 title '1 M' \
    ,file3p1 us ($1/day):($15) w l lt 3 title '3 M' \

    

reset





# Figure 10
# --------
set output 'IceSaturationVsTime.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Time (day)" font ",24"
x0 = 1.e-2
x1 = 100
set xrange [x0:x1] noreverse nowriteback
set xtics out 10 nomirror
#set autoscale x
set logscale x

# y-axis
set ylabel "Ice saturation" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 1.e-1 nomirror
set autoscale y
set tics out

set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5


plot \
     file0p1 us ($1/day):(1-$5) w l lt 1 title '0 M' \
    ,file1p1 us ($1/day):(1-$5) w l lt 2 title '1 M' \
    ,file3p1 us ($1/day):(1-$5) w l lt 3 title '3 M' \

    

reset





# Figure 11
# --------
set output 'IcePressureVsTime.eps'

# Size
set size square 0.75,1.
set origin 0,0



# x-axis
set xlabel "Time (day)" font ",24"
x0 = 1.e-2
x1 = 100
set xrange [x0:x1] noreverse nowriteback
set xtics out 10 nomirror
#set autoscale x
set logscale x

# y-axis
set ylabel "Ice pressure (MPa)" font ",24"
y0 = 0
y1 = 100
set yrange [y0:y1] noreverse nowriteback
set ytics out 20 nomirror
set autoscale y
set tics out

set grid


# Legends
set key top right reverse Left samplen 1 spacing 1.5 maxrows 5


plot \
     file0p1 us ($1/day):($12/MPa) w l lt 1 title '0 M' \
    ,file1p1 us ($1/day):($12/MPa) w l lt 2 title '1 M' \
    ,file3p1 us ($1/day):($12/MPa) w l lt 3 title '3 M' \

    

reset

exit


     
