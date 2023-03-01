# Time(1) Displacements(2) Stresses(5) Perturbated-displacements(14) Plastic-strains(17) Hardening variable(26) Yield function(27)


# Linetypes
# ---------
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


# To produce eps file
set term postscript eps enhanced color 20




# Units:
kPa = 1
MPa = 1.e3
Pa  = 1.e-3



# Yield function:
m = 0.9
m2 = m*m
pc = 0.1*MPa
Y(p,pc) = (p <= pc) ? m*sqrt(abs(p*(p - pc))) : 1/0




# Input Data file
file = 'toto'

# Data files
file1 = file.'.p1'



# 1st Plot
# --------

set output 'P-Q.eps'

# Legends
set key top right reverse Left samplen 2 spacing 2

# Size
set size square 1.,1.

# x-axis
set xlabel 'p (MPa)' font ",24"
x0 = 0
x1 = 2*pc
set autoscale x
set xrange [x0:x1] noreverse nowriteback
set xtics 0.02*MPa
#set format x "%1.0e"

# y-axis
set ylabel 'q (MPa)' font ",24"
y0 = 0
y1 = m*pc*3/4
set autoscale y
set yrange[y0:y1] noreverse nowriteback
#set ytics norotate 0,0.2,1.2

set grid

plot \
     Y(x,pc) lt 1 title 'Yield function' \
    ,Y(x,1.5*pc)  lt 3 notitle \
    ,file1 us (-($5+$9+$13)/3*Pa):(abs($5-$9)*Pa) lt 2 title 'p-q'

reset


# 2nd Plot
# --------
set output 'E-Q.eps'

# Legends
set key top right reverse Left samplen 2 spacing 2

# Size
set size square 1.,1.

# x-axis
set xlabel '{/Symbol e}@_q^p (left) and {/Symbol e}@_v^p (right) (10^{-3})' font ",24"
x0 = 0
x1 = 2*pc
set autoscale x
#set xrange [x0:x1] noreverse nowriteback
#set xtics 0.02*MPa
#set format x "%1.0e"

# y-axis
set ylabel 'q (MPa)' font ",24"
y0 = 0
y1 = m*pc*3/4
set autoscale y
set yrange[y0:y1] noreverse nowriteback
#set ytics norotate 0,0.2,1.2

set grid

plot \
     file1 us (-($17+$21+$25)*1.e3):(abs($5-$9)*Pa) ls 1 title '{/Symbol e}@_v^p-q' \
    ,file1 us (-abs($17-$21))*1.e3:(abs($5-$9)*Pa) ls 2 title '{/Symbol e}@_q^p-q'


