#!/bin/sh

find ${BIL}/src -maxdepth 3 -name "*.[c,h,cpp,hpp,f]" | xargs wc -l
