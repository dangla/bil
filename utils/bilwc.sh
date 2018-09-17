#!/bin/sh

find ${BIL}/src -name "*.[c,h,cpp,hpp,f]" | xargs wc -l
