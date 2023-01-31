reset

set term pdfcairo linewidth 2 font ",14"
set output 'elnp.pdf'

#set size square 0.75,1.
set xlabel 'p (kPa)'
set ylabel '{/Symbol D} e (-)'
set grid
set key samplen 1 spacing 4

set logscale x
plot 'test_m15.p1' using (-$21/1000):($20 + $19) notitle \
  with lines lt 1 lw 2

set output # required for pdfcairo
