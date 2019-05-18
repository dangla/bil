#!/bin/sh

find ${BIL} -name "*.[c,h,cpp,hpp,f]" | xargs wc -l
