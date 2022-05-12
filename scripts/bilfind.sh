#!/bin/sh
# Use [[:blank:]] as a regular expression for a blank.

grep --exclude-dir='backup*' --include='*.[c,h,f,y,l]' --include='*.sh' --include='*.gp' --include='*.py' --include='Makefile' --include='CMakeLists.txt' --include='*.cmake' --exclude='*.o' -r -T -E --color -e $1 ${BIL}

