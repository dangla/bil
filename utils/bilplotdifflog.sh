#!/bin/sh

gnuplot -e "call '${BIL}/utils/bilplotdifflog.gp' \"$1\" \"$2\" \"$3\" "

