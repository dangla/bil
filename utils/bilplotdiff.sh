#!/bin/sh

gnuplot -e "call '${BIL}/utils/bilplotdiff.gp' \"$1\" \"$2\" \"$3\" "

