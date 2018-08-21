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



# Linear function
# ---------------
Linear(T,A,B) = ((343 - T)*(A) + (T - 293)*(B))/50.


# Equilibrium constants
# ---------------------
logk_ch(T)    = Linear(T,-5.14,-5.80)   # CH    = Ca[2+] + 2OH[-]
logk_csh2(T)  = Linear(T,-4.58,-4.58)   # CSH2  = Ca[2+] + SO4[2-] + 2H2O
logk_ah3(T)   = Linear(T,-68,-68)       # AH3   = 2Al[3+] + 6OH[-]
logk_afm(T)   = Linear(T,-96,-94)       # AFm   = 4Ca[2+] + 2Al[3+] + SO4[2-] + 12OH[-] + 6H2O
logk_aft(T)   = Linear(T,-112,-105)     # AFt   = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 12OH[-] + 26H2O
logk_c3ah6(T) = Linear(T,-89.75,-73.46) # C3AH6 = 3Ca[2+] + 2Al[3+] + 12OH[-]
logk_as3h6(T) = Linear(T,0.31,0.31)     # AS3H6 = 2Al[3+] + 3SO4[2-] + 6H2O
logk_h2so4(T) = Linear(T,-32,-29)       # SO4[-2] + 2H2O = H2SO4 + 2 OH[-]



# Log(S_i) as function of log(S_CH), log(S_AH3) and log(a_H2SO4 / k)
# ------------------------------------------------------------------

# Gypsum
LS_CSH2(ls_ch,la_h2so4,T)       = - logk_csh2(T) + (ls_ch + logk_ch(T)) + (la_h2so4 - logk_h2so4(T))
  
# Aluminium-sulfur compounds
LS_AS3H6(ls_ah3,la_h2so4,T)     = - logk_as3h6(T) + (ls_ah3 + logk_ah3(T)) + 3*(la_h2so4 - logk_h2so4(T))

# Monosulfoaluminate
LS_AFm(ls_ch,ls_ah3,la_h2so4,T) = - logk_afm(T) + 4*(ls_ch + logk_ch(T)) + (ls_ah3 + logk_ah3(T)) + (la_h2so4 - logk_h2so4(T))

# Ettringite
LS_AFt(ls_ch,ls_ah3,la_h2so4,T) = - logk_aft(T) + 6*(ls_ch + logk_ch(T)) + (ls_ah3 + logk_ah3(T)) + 3*(la_h2so4 - logk_h2so4(T))
  
# Hydrogarnets
LS_C3AH6(ls_ch,ls_ah3,T)        = - logk_c3ah6(T) + 3*(ls_ch + logk_ch(T)) + (ls_ah3 + logk_ah3(T))


# The following surfaces, LS_AH3_i(), give the value of ls_ah3 
# for which log(S_i) = 0.
# The undersaturated domain (given by log(S_i) < 0) is therefore beneath these surfaces
LS_AH3_AFm(ls_ch,la_h2so4,T)   = - LS_AFm(ls_ch,0,la_h2so4,T)
LS_AH3_AFt(ls_ch,la_h2so4,T)   = - LS_AFt(ls_ch,0,la_h2so4,T)
LS_AH3_C3AH6(ls_ch,T)          = - LS_C3AH6(ls_ch,0,T)
LS_AH3_AS3H6(la_h2so4,T)       = - LS_AS3H6(0,la_h2so4,T)
# The following surfaces, LS_CH_i(), give the value of ls_ch 
# for which log(S_i) = 0.
# The undersaturated domain is therefore beneath these surfaces
LS_CH_CSH2(la_h2so4,T)         = - LS_CSH2(0,la_h2so4,T)


# We define the surface as the Min of the previous surfaces
Surface1(a,b,T) = Min5(LS_AH3_AFm(a,b,T),LS_AH3_AFt(a,b,T),LS_AH3_C3AH6(a,T),LS_AH3_AS3H6(b,T),0)
# The same but without AFt
Surface2(a,b,T) = Min4(LS_AH3_AFm(a,b,T),LS_AH3_C3AH6(a,T),LS_AH3_AS3H6(b,T),0)


# We cap the surface by the saturation degree of CSH2
Surface(a,b,T) = ((LS_CSH2(a,b,T) > 0) ? 1/0 : Surface1(a,b,T))
#Surface(a,b,T) = ((LS_CSH2(a,b,T) > 0) ? 1/0 : Surface2(a,b,T))

# The different parts of the surface
LS_AH3_AFm_(a,b,T)   = ((LS_AH3_AFm(a,b,T)   > Surface(a,b,T)) ? 1/0 : LS_AH3_AFm(a,b,T))
LS_AH3_AFt_(a,b,T)   = ((LS_AH3_AFt(a,b,T)   > Surface(a,b,T)) ? 1/0 : LS_AH3_AFt(a,b,T))
LS_AH3_C3AH6_(a,b,T) = ((LS_AH3_C3AH6(a,T)   > Surface(a,b,T)) ? 1/0 : LS_AH3_C3AH6(a,T))
LS_AH3_AS3H6_(a,b,T) = ((LS_AH3_AS3H6(b,T)   > Surface(a,b,T)) ? 1/0 : LS_AH3_AS3H6(b,T))
LS_AH3_AH3_(a,b,T)   = ((0                   > Surface(a,b,T)) ? 1/0 : 0                )



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

# Invariant points of Ettringite
# ------------------------------
# P1: X = 0, LS_C3AH6 = 0, LS_AFt = 0
X1(T) = 0
Z1(T) = - (3*logk_ch(T) - logk_c3ah6(T) + logk_ah3(T) + 3*X1(T))
Y1(T) = - (6*logk_ch(T) - logk_aft(T)   + logk_ah3(T) + 6*X1(T) + Z1(T))/3 + logk_h2so4(T)

# P2: LS_C3AH6 = 0, Z = 0, LS_AFt = 0
Z2(T) = 0
X2(T) = - (3*logk_ch(T) - logk_c3ah6(T) + logk_ah3(T) + Z2(T))/3
Y2(T) = - (6*logk_ch(T) - logk_aft(T)   + logk_ah3(T) + 6*X2(T) + Z2(T))/3 + logk_h2so4(T)

# P3: LS_CSH2 = 0, Z = 0, LS_AFt = 0
#0 = logk_ch   - logk_csh2           + X + Y - logk_h2so4
#0 = 6*logk_ch - logk_aft   + logk_ah3 + 6*X + Z + 3*(Y - logk_h2so4)
#0 = 3*logk_ch - logk_aft   + logk_ah3 + 3*X + Z - 3*(- logk_csh2)
Z3(T) = 0
X3(T) = - (3*logk_ch(T) - logk_aft(T)   + logk_ah3(T) + Z3(T) - 3*(- logk_csh2(T)))/3
Y3(T) = - (logk_ch(T)   - logk_csh2(T)  + X3(T)) + logk_h2so4(T)

# P4: X = 0 , LSCSH2 = 0 , LS_AFt = 0
X4(T) = 0
Y4(T) = - (logk_ch(T) - logk_csh2(T)  + X4(T)) + logk_h2so4(T)
Z4(T) = - (6*logk_ch(T) - logk_aft(T)   + logk_ah3(T) + 6*X4(T) + 3*(Y4(T) - logk_h2so4(T)))

# Invariant point of Gibbsite-Gipsum-AS3H6
# ----------------------------------------
# P5: LS_CSH2 = 0, Z = 0, LS_AS3H6 = 0
Z5(T) = 0
Y5(T) = - (logk_ah3(T) - logk_as3h6(T) + Z5(T))/3 + logk_h2so4(T)
X5(T) = - (logk_ch(T) - logk_csh2(T)  + Y5(T) - logk_h2so4(T))



# Invariant points of Monosulfoaluminate
# --------------------------------------
# Q1: X = 0, LS_C3AH6 = 0, LS_AFm = 0
Xm1(T) = 0
Zm1(T) = - (3*logk_ch(T) - logk_c3ah6(T) + logk_ah3(T) + 3*Xm1(T))
Ym1(T) = - (4*logk_ch(T) - logk_afm(T)   + logk_ah3(T) + 4*Xm1(T) + Zm1(T)) + logk_h2so4(T)

# Q2: LS_C3AH6 = 0, Z = 0, LS_AFm = 0
Zm2(T) = 0
Xm2(T) = - (3*logk_ch(T) - logk_c3ah6(T) + logk_ah3(T) + Zm2(T))/3
Ym2(T) = - (4*logk_ch(T) - logk_afm(T)   + logk_ah3(T) + 4*Xm2(T) + Zm2(T)) + logk_h2so4(T)

# Q3: LS_CSH2 = 0, Z = 0, LS_AFm = 0
Zm3(T) = 0
Xm3(T) = - (3*logk_ch(T) - logk_afm(T)   + logk_ah3(T) + Zm3(T) - (- logk_csh2(T)))/3
Ym3(T) = - (logk_ch(T)   - logk_csh2(T)  + Xm3(T)) + logk_h2so4(T)

# Q4: X = 0 , LSCSH2 = 0 , LS_AFm = 0
Xm4(T) = 0
Ym4(T) = - (logk_ch(T) - logk_csh2(T)  + Xm4(T)) + logk_h2so4(T)
Zm4(T) = - (4*logk_ch(T) - logk_afm(T)   + logk_ah3(T) + 4*Xm4(T) + (Ym4(T) - logk_h2so4(T)))











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
set size square 1.,1

T = 293


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
y0 = -8 + logk_h2so4(T)
y1 = +8 + logk_h2so4(T)
#y1 = +32 + logk_h2so4(T)
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
set arrow from X1(T),Y1(T),Z1(T) to X2(T),Y2(T),Z2(T) nohead lw 2 front
set arrow from X2(T),Y2(T),Z2(T) to X3(T),Y3(T),Z3(T) nohead lw 2 front
set arrow from X3(T),Y3(T),Z3(T) to X4(T),Y4(T),Z4(T) nohead lw 2 front
set arrow from X1(T),Y1(T),Z1(T) to X4(T),Y4(T),Z4(T) nohead lw 2 front
set arrow from X1(T),Y1(T),Z1(T) to X1(T),y0,Z1(T) nohead lw 2 front
set arrow from X2(T),Y2(T),Z2(T) to X2(T),y0,Z2(T) nohead lw 2 front
set arrow from X3(T),Y3(T),Z3(T) to x0,Y3(T)-x0+X3(T),Z3(T) nohead lw 2 front
#set arrow from X3(T),Y3(T),Z3(T) to X5(T),Y5(T),Z5(T) nohead lw 2 front
set arrow from X4(T),Y4(T),Z4(T) to X4(T),Y4(T),z0 nohead lw 2 front
#set arrow from X5(T),Y5(T),Z5(T) to x0,Y5(T),Z5(T) nohead lw 2 front
#set arrow from X5(T),Y5(T),Z5(T) to X5(T) + z0/3,Y5(T) -z0/3,z0 nohead lw 2 front
set arrow from Xm1(T),Ym1(T),Zm1(T) to Xm2(T),Ym2(T),Zm2(T) nohead lt 0 lw 4 front
set arrow from Xm2(T),Ym2(T),Zm2(T) to Xm3(T),Ym3(T),Zm3(T) nohead lt 0 lw 4 front
set arrow from Xm3(T),Ym3(T),Zm3(T) to Xm4(T),Ym4(T),Zm4(T) nohead lt 0 lw 4 front
set arrow from Xm1(T),Ym1(T),Zm1(T) to Xm4(T),Ym4(T),Zm4(T) nohead lt 0 lw 4 front

set arrow from X1(T),Y1(T),Z1(T) to Xm1(T),Ym1(T),Zm1(T) nohead lt 0 lw 4 front
set arrow from X2(T),Y2(T),Z2(T) to Xm2(T),Ym2(T),Zm2(T) nohead lt 0 lw 4 front
set arrow from X3(T),Y3(T),Z3(T) to Xm3(T),Ym3(T),Zm3(T) nohead lt 0 lw 4 front
set arrow from X4(T),Y4(T),Z4(T) to Xm4(T),Ym4(T),Zm4(T) nohead lt 0 lw 4 front



set label "Portlandite" at  0,y0,-12 left   front textcolor rgbcolor "black" font ",24"
set label "Ettringite"  at -1,-30,-4  center front textcolor rgbcolor "black" font ",24"
#set label "Monosulfo"  at -1,-31,-4  center front textcolor rgbcolor "black" font ",24"
set label "Gipsum"      at -2,-26,-8  center front textcolor rgbcolor "black" font ",24"
set label "Hydrogarnet" at  0,y0,-2  left   front textcolor rgbcolor "black" font ",24"
set label "Gibbsite"    at -4,y0+2,0   center front textcolor rgbcolor "black" font ",24"


set output 'PhaseDiagramOfCaO-Al2O3-SO3-H2O.eps'
set isosamples 100,150
#set isosamples 80,100

splot \
      (LS_AH3_AFt_(x,y,T))    with lines lt 2 title 'AFt' \
    , (LS_AH3_AH3_(x,y,T))    with lines lt 5 title 'AH3' \
    , (LS_AH3_AS3H6_(x,y,T))  with lines lt 1 title 'AS3H6' \
    , (LS_AH3_C3AH6_(x,y,T))  with lines lt 4 title 'C_3AH_6' \
    , (LS_AH3_AFm_(x,y,T))    with lines lt 3 title 'AFm' \
    , '++' using (0):($2):((($2 < logk_csh2(T) - logk_ch(T) + logk_h2so4(T)) && ((z0*$1)/x0 <= Surface(0,$2,T))) ? (z0*$1)/x0 : 1/0) with lines lt 6 title 'CH' \
    , '++' using (LS_CH_CSH2($2,T)):($2):((($2 > logk_csh2(T) - logk_ch(T) + logk_h2so4(T)) && ((z0*$1)/x0 <= Surface(LS_CH_CSH2($2,T),$2,T))) ? (z0*$1)/x0 : 1/0) with lines lt 8 title 'CSH_2' \



reset



# Plot 2
#-------

# Size
set size square 1,1


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
y0 = -8 + logk_h2so4(T)
y1 = +8 + logk_h2so4(T)
#y1 = +32 + logk_h2so4(T)
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


set hidden3d

# Legends
#set key lmargin  horizontal Left reverse samplen 2 spacing 1
#set key lmargin  horizontal Right samplen 2 spacing 1
#set key tmargin left  horizontal Right samplen 2 spacing 1 maxrows 2
set nokey


set term gif animate delay 30
#set term postscript eps enhanced color 20
set output 'PhaseDiagramOfCaO-Al2O3-SO3-H2O.gif'
set isosamples 100,150
#set isosamples 80,100

N = 20
do for [i = 0:N] {
T = 293 + (70. * i)/N
splot \
      (LS_AH3_AFt_(x,y,T))    with lines lt 2 title 'AFt' \
    , (LS_AH3_AH3_(x,y,T))    with lines lt 5 title 'AH3' \
    , (LS_AH3_AS3H6_(x,y,T))  with lines lt 1 title 'AS3H6' \
    , (LS_AH3_C3AH6_(x,y,T))  with lines lt 4 title 'C_3AH_6' \
    , (LS_AH3_AFm_(x,y,T))    with lines lt 3 title 'AFm' \
    , '++' using (0):($2):((($2 < logk_csh2(T) - logk_ch(T) + logk_h2so4(T)) && ((z0*$1)/x0 <= Surface(0,$2,T))) ? (z0*$1)/x0 : 1/0) with lines lt 6 title 'CH' \
    , '++' using (LS_CH_CSH2($2,T)):($2):((($2 > logk_csh2(T) - logk_ch(T) + logk_h2so4(T)) && ((z0*$1)/x0 <= Surface(LS_CH_CSH2($2,T),$2,T))) ? (z0*$1)/x0 : 1/0) with lines lt 8 title 'CSH_2' \

}


