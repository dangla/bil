set term wxt
tictac = 1

file1 = "$0"
file2 = "$1"

do for [i = 2:99] {
plot   file1 us 1:i title sprintf("%s i=%d",file1,i) \
      ,file2 us 1:i title sprintf("%s i=%d",file2,i)
pause tictac
}

