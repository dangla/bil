#!/bin/sh

ar -x $1.a
rm -f $1.so
gcc -shared -gdwarf-2 -g3 -Wl,-soname=$1.so -o $1.so *.o $2
rm -f *.o

