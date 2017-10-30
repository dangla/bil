#!/bin/sh

grep --exclude-dir=backup --include="*.[c,h,cpp,hpp]" --include="Makefile" --exclude="*.o" -r -T -E --color -e "$1" ${BIL}/src

