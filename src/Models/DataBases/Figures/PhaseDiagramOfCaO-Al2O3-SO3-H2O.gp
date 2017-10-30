# Gunplot file. To produce the figure type
#
# gnuplot PhaseDiagramOfCaO-Al2O3-SO3-H2O.gp
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



# Min functions
# -------------
Min(a,b) = ((a < b) ? a : b)
Min3(a,b,c) = Min(a,Min(b,c))
Min4(a,b,c,d) = Min(a,Min3(b,c,d))
Min5(a,b,c,d,e) = Min(a,Min4(b,c,d,e))

Minc(a,b,A,B) = ((a < b) ? A : B )
Minc3(a,b,c,A,B,C) = (a < Min(b,c) ? A : Minc(b,c,B,C))
Minc4(a,b,c,d,A,B,C,D) = (a < Min3(b,c,d) ? A : Minc3(b,c,d,B,C,D))
Minc5(a,b,c,d,e,A,B,C,D,E) = (a < Min4(b,c,d,e) ? A : Minc4(b,c,d,e,B,C,D,E))



# Max functions
# -------------
Max(a,b) = ((a > b) ? a : b)
Max3(a,b,c) = Max(a,Max(b,c))
Max4(a,b,c,d) = Max(a,Max3(b,c,d))
Max5(a,b,c,d,e) = Max(a,Max4(b,c,d,e))

Maxc(a,b,A,B) = ((a > b) ? A : B )
Maxc3(a,b,c,A,B,C) = (a < Max(b,c) ? A : Maxc(b,c,B,C))
Maxc4(a,b,c,d,A,B,C,D) = (a < Max3(b,c,d) ? A : Maxc3(b,c,d,B,C,D))
Maxc5(a,b,c,d,e,A,B,C,D,E) = (a < Max4(b,c,d,e) ? A : Maxc4(b,c,d,e,B,C,D,E))


# Equilibrium constants
# ---------------------
logk_ch    = -5.14
logk_csh2  = -4.58
logk_ah3   = -68
logk_afm   = -96
logk_aft   = -112
logk_c3ah6 = -89.75
logk_as3h6 = 0.31
logk_h2so4 = -32


# Log(S_i) as function of log(S_CH), log(S_AH3) and log(a_H2SO4 / k)
# ------------------------------------------------------------------

# Gypsum
LS_CSH2(ls_ch,la_h2so4)       = - logk_csh2 + (ls_ch + logk_ch) + (la_h2so4 - logk_h2so4)
  
# Aluminium-sulfur compounds
LS_AS3H6(ls_ah3,la_h2so4)     = - logk_as3h6 + (ls_ah3 + logk_ah3) + 3*(la_h2so4 - logk_h2so4)

# Monosulfoaluminate
LS_AFm(ls_ch,ls_ah3,la_h2so4) = - logk_afm + 4*(ls_ch + logk_ch) + (ls_ah3 + logk_ah3) + (la_h2so4 - logk_h2so4)

# Ettringite
LS_AFt(ls_ch,ls_ah3,la_h2so4) = - logk_aft + 6*(ls_ch + logk_ch) + (ls_ah3 + logk_ah3) + 3*(la_h2so4 - logk_h2so4)
  
# Hydrogarnets
LS_C3AH6(ls_ch,ls_ah3)        = - logk_c3ah6 + 3*(ls_ch + logk_ch) + (ls_ah3 + logk_ah3)


# The following surfaces, LS_AH3_i(), give the value of ls_ah3 
# for which log(S_i) = 0.
# The undersaturated domain (given by log(S_i) < 0) is therefore beneath these surfaces
LS_AH3_AFm(ls_ch,la_h2so4)   = - LS_AFm(ls_ch,0,la_h2so4)
LS_AH3_AFt(ls_ch,la_h2so4)   = - LS_AFt(ls_ch,0,la_h2so4)
LS_AH3_C3AH6(ls_ch)          = - LS_C3AH6(ls_ch,0)
LS_AH3_AS3H6(la_h2so4)       = - LS_AS3H6(0,la_h2so4)
# The following surfaces, LS_CH_i(), give the value of ls_ch 
# for which log(S_i) = 0.
# The undersaturated domain is therefore beneath these surfaces
LS_CH_CSH2(la_h2so4)         = - LS_CSH2(0,la_h2so4)


# We define the surface as the Min of the previous surfaces
Surface1(a,b) = Min5(LS_AH3_AFm(a,b),LS_AH3_AFt(a,b),LS_AH3_C3AH6(a),LS_AH3_AS3H6(b),0)
# The same but without AFt
Surface2(a,b) = Min4(LS_AH3_AFm(a,b),LS_AH3_C3AH6(a),LS_AH3_AS3H6(b),0)


# We cap the surface by the saturation degree of CSH2
Surface(a,b) = ((LS_CSH2(a,b) > 0) ? 1/0 : Surface1(a,b))
#Surface(a,b) = ((LS_CSH2(a,b) > 0) ? 1/0 : Surface2(a,b))

# The different parts of the surface
LS_AH3_AFm_(a,b)   = ((LS_AH3_AFm(a,b)   > Surface(a,b)) ? 1/0 : LS_AH3_AFm(a,b))
LS_AH3_AFt_(a,b)   = ((LS_AH3_AFt(a,b)   > Surface(a,b)) ? 1/0 : LS_AH3_AFt(a,b))
LS_AH3_C3AH6_(a,b) = ((LS_AH3_C3AH6(a)   > Surface(a,b)) ? 1/0 : LS_AH3_C3AH6(a))
LS_AH3_AS3H6_(a,b) = ((LS_AH3_AS3H6(b)   > Surface(a,b)) ? 1/0 : LS_AH3_AS3H6(b))
LS_AH3_AH3_(a,b)   = ((0                 > Surface(a,b)) ? 1/0 : 0              )



# Invariant points
# ----------------
# LS_CH    = X
# LA_H2SO4 = Y
# LS_AH3   = Z
# LS_CSH2  = logk_ch   - logk_csh2             +   X +       (Y - logk_h2so4)
# LS_AS3H6 = logk_ah3  - logk_as3h6                  + Z + 3*(Y - logk_h2so4) 
# LS_AFm   = 4*logk_ch - logk_afm   + logk_ah3 + 4*X + Z +   (Y - logk_h2so4)
# LS_AFt   = 6*logk_ch - logk_aft   + logk_ah3 + 6*X + Z + 3*(Y - logk_h2so4)
# LS_C3AH6 = 3*logk_ch - logk_c3ah6 + logk_ah3 + 3*X + Z

# P1: X = 0, LS_C3AH6 = 0, LS_AFt = 0
X1 = 0
Z1 = - (3*logk_ch - logk_c3ah6 + logk_ah3 + 3*X1)
Y1 = - (6*logk_ch - logk_aft   + logk_ah3 + 6*X1 + Z1)/3 + logk_h2so4

# P2: LS_C3AH6 = 0, Z = 0, LS_AFt = 0
Z2 = 0
X2 = -(3*logk_ch - logk_c3ah6 + logk_ah3 + Z2)/3
Y2 = - (6*logk_ch - logk_aft   + logk_ah3 + 6*X2 + Z2)/3 + logk_h2so4

# P3: LS_CSH2 = 0, Z = 0, LS_AFt = 0
#0 = logk_ch   - logk_csh2           + X + Y - logk_h2so4
#0 = 6*logk_ch - logk_aft   + logk_ah3 + 6*X + Z + 3*(Y - logk_h2so4)
#0 = 3*logk_ch - logk_aft   + logk_ah3 + 3*X + Z - 3*(- logk_csh2)
Z3 = 0
X3 = - (3*logk_ch - logk_aft   + logk_ah3 + Z3 - 3*(- logk_csh2))/3
Y3 = - (logk_ch   - logk_csh2  + X3) + logk_h2so4

# P4: X = 0 , LSCSH2 = 0 , LS_AFt = 0
X4 = 0
Y4 = - (logk_ch - logk_csh2  + X4) + logk_h2so4
Z4 = - (6*logk_ch - logk_aft   + logk_ah3 + 6*X4 + 3*(Y4 - logk_h2so4))

# P5: LS_CSH2 = 0, Z = 0, LS_AS3H6 = 0
Z5 = 0
Y5 = - (logk_ah3 - logk_as3h6 + Z5)/3 + logk_h2so4
X5 = - (logk_ch - logk_csh2  + Y5 - logk_h2so4)













# Styles
#set pointsize 2
#set style line 1  lw 1 lt 1 pi -10 ps 2 pt 4 lc rgb "red"
#set style line 2  lw 1 lt 2 pi -10 ps 2 pt 8 lc rgb "light-magenta"
#set style line 3  lw 1 lt 3 pi -10 ps 2 pt 6 lc rgb "purple"
#set style line 4  lw 1 lt 4 pi -10 ps 2 pt 4 lc rgb "steelblue"
#set style line 5  lw 1 lt 5 pi -10 ps 2 pt 8 lc rgb "bisque"
#set style line 6  lw 1 lt 6 pi -10 ps 2 pt 6 lc rgb "green"
#set style line 6  lw 1 lt 7 pi -10 ps 2 pt 6 lc rgb "green"
#set style line 7  lw 1 lt 8 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"
#set style line 7  lw 1 lt 9 pi -10 ps 2 pt 4 lc rgb "light-goldenrod"




# Plot 1
#-------
# Size
set size square 1.3,1


# x-axis
#set xlabel 'log({/Symbol b}_{CH})' font ",24" offset 2,0 rotate parallel
set xlabel 'Saturation index of {CH}' font ",24" offset 1,-1 rotate parallel
#x0 = -32
x0 = -6
x1 = 0
set xrange [x0:x1] noreverse nowriteback
set xtics nomirror out norotate x0,2,x1
set format x "10^{%+-3.0f}"


# y-axis
#set ylabel 'log(a_{H_2SO_4})' font ",24" offset -1.5,-1.5 #rotate parallel
set ylabel 'Activity of {H_2SO_4}' font ",24" offset 0,-2 rotate parallel
y0 = -8 + logk_h2so4
y1 = +8 + logk_h2so4
#y1 = +32 + logk_h2so4
set yrange[y0:y1] noreverse nowriteback
set ytics  nomirror out scale 4. norotate offset 0,-0.5 y0,4,y1
set format y "10^{%+-3.0f}"


# z-axis
#set zlabel 'log({/Symbol b}_{AH3})' font ",24" offset 0,0 rotate parallel
set zlabel 'Saturation index of {AH3}' font ",24" offset 0,0 rotate parallel
z0 = -16
z1 = 0
set zrange[z0:z1] noreverse nowriteback
set ztics nomirror out norotate z0,4,z1
set format z "10^{%+-3.0f}"

#set view 60, 30, 1, 1 # default
set view 70, 105, 1,1
#set view equal xy

set xyplane at z0

set samples 100
set isosamples 50,80

set contour surface
set cntrparam levels discrete -7 #incremental -7,40,-7
set clabel #font "bold,14"
unset clabel
unset contour

#set table 'StableSurface.dat'
#splot  Surface(x,y) with lines ls 100 notitle
#unset table

set hidden3d

# Legends
#set key lmargin  horizontal Left reverse samplen 2 spacing 1
#set key lmargin  horizontal Right samplen 2 spacing 1
#set key tmargin left  horizontal Right samplen 2 spacing 1 maxrows 2
set nokey
    

# Plotting the lines joining the invariant points
set arrow from X1,Y1,Z1 to X2,Y2,Z2 nohead lw 2 front
set arrow from X2,Y2,Z2 to X3,Y3,Z3 nohead lw 2 front
set arrow from X3,Y3,Z3 to X4,Y4,Z4 nohead lw 2 front
set arrow from X1,Y1,Z1 to X4,Y4,Z4 nohead lw 2 front
set arrow from X1,Y1,Z1 to X1,y0,Z1 nohead lw 2 front
set arrow from X2,Y2,Z2 to X2,y0,Z2 nohead lw 2 front
set arrow from X3,Y3,Z3 to x0,Y3-x0+X3,Z3 nohead lw 2 front
#set arrow from X3,Y3,Z3 to X5,Y5,Z5 nohead lw 2 front
set arrow from X4,Y4,Z4 to X4,Y4,z0 nohead lw 2 front
#set arrow from X5,Y5,Z5 to x0,Y5,Z5 nohead lw 2 front
#set arrow from X5,Y5,Z5 to X5 + z0/3,Y5 -z0/3,z0 nohead lw 2 front


set label "Portlandite" at 0,-37,-12 center front textcolor rgbcolor "black" font ",24"
set label "Ettringite" at -1,-32,-8 center front textcolor rgbcolor "black" font ",24"
#set label "Monosulfo" at -1,-31,-4 center front textcolor rgbcolor "black" font ",24"
set label "Gipsum" at -2,-26,-8 center front textcolor rgbcolor "black" font ",24"
set label "Hydrogarnet" at -1,-37,-3 center front textcolor rgbcolor "black" font ",24"
set label "Gibbsite" at -4,-36,0 center front textcolor rgbcolor "black" font ",24"


set output 'PhaseDiagramOfCaO-Al2O3-SO3-H2O.eps'
set isosamples 100,150
#set isosamples 80,100

splot \
      (LS_AH3_AFt_(x,y))    with lines lt 2 title 'AFt' \
    , (LS_AH3_AH3_(x,y))    with lines lt 5 title 'AH3' \
    , (LS_AH3_AS3H6_(x,y))  with lines lt 1 title 'AS3H6' \
    , (LS_AH3_C3AH6_(x,y))  with lines lt 4 title 'C_3AH_6' \
    , (LS_AH3_AFm_(x,y))    with lines lt 3 title 'AFm' \
    , '++' using (0):($2):((($2 < logk_csh2 - logk_ch + logk_h2so4) && ((z0*$1)/x0 <= Surface(0,$2))) ? (z0*$1)/x0 : 1/0) with lines lt 6 title 'CH' \
    , '++' using (LS_CH_CSH2($2)):($2):((($2 > logk_csh2 - logk_ch + logk_h2so4) && ((z0*$1)/x0 <= Surface(LS_CH_CSH2($2),$2))) ? (z0*$1)/x0 : 1/0) with lines lt 8 title 'CSH_2' \


reset

