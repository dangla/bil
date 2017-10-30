#!/bin/sh

find ${BIL}/src -name "*.[c,h,cpp,hpp]" | xargs wc -l
