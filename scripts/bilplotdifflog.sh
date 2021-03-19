#!/bin/sh

gnuplot -e "call '${BIL}/scripts/bilplotdifflog.gp' \"$1\" \"$2\" \"$3\" "

