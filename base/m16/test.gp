set term wxt
tictac = 1

file1 = "m16.t1"
file2 = "toto.t1"

do for [i = 4:29] {
plot   file1 us 1:i title sprintf("%s i=%d",file1,i) \
      ,file2 us 1:i title sprintf("%s i=%d",file2,i)
pause tictac
}
