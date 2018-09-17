#!/bin/sh

grep --exclude-dir=backup --include="*.[c,h,cpp,hpp,f]" --include="Makefile" --exclude="*.o" -r -T -E --color -e "$1" ${BIL}/src

