#!/bin/sh

find ${BIL} -name $1 | xargs less -N

