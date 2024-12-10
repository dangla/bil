#!/bin/sh
# Use [[:blank:]] as a regular expression for a blank.
# $1 stands for a STRING to look for.
# Work with $1 = "STRING -A N".

grep --exclude-dir='backup*' --exclude-dir='Libraries' --include='*.[chfyl]' --include='*.[ch]pp' --include='*.sh' --include='*.gp' --include='*.py' --include='Makefile' --include='CMakeLists.txt' --include='*.cmake' --exclude='*.o' -r -T -E --color -e $1 ${BIL}

