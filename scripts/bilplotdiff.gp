set term qt
#set term x11 
#set term wxt

tictac = 1

NbOfFiles = ARGC

#do for [i = 1:NbOfFiles] {
#filei = "ARG" . i
#}
file1 = ARG1
file2 = ARG2

if(ARG3 == "") {
  tictac = 1 ;
} else {
  tictac = ARG3 ;
}

# Styles
set pointsize 2
set style line 10  lw 2 lt 1 pi -10 ps 2 pt 1 lc rgb "blue"
set style line 20  lw 2 lt 1 pi -10 ps 2 pt 19 lc rgb "red"



# Linetypes
set linetype 1 lw 2 pi -10 ps 3  pt 1  linecolor rgb "red"
set linetype 2 lw 2 pi -10 ps 3  pt 6  linecolor rgb "light-magenta"
set linetype 3 lw 2 pi -10 ps 2  pt 2  linecolor rgb "purple"
set linetype 4 lw 2 pi -10 ps 2  pt 3  linecolor rgb "steelblue"
set linetype 5 lw 2 pi -10 ps 2  pt 4  linecolor rgb "aquamarine"
set linetype 6 lw 2 pi -10 ps 2  pt 5  linecolor rgb "bisque"
set linetype 7 lw 2 pi -10 ps 2  pt 6  linecolor rgb "bisque"
set linetype 8 lw 2 pi -10 ps 2  pt 7  linecolor rgb "light-goldenrod"
set linetype 9 lw 2 pi -10 ps 2  pt 8  linecolor rgb "light-goldenrod"


do for [i = 2:99] {
plot   file1 us 1:i w lp lt 1 title sprintf("%s i=%d",file1,i) \
      ,file2 us 1:i w lp lt 2 title sprintf("%s i=%d",file2,i)
      
pause tictac
}

