#!/bin/sh

grep --exclude-dir='backup*' --include='*.[c,h,f]' --include='*.sh' --include='*.gp' --include='*.py' --include='Makefile' --exclude='*.o' -r -T -E --color -e $1 ${BIL}

