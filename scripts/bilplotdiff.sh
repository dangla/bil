#!/bin/sh

gnuplot -e "call '${BIL}/scripts/bilplotdiff.gp' \"$1\" \"$2\" \"$3\" "

