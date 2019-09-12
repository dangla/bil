#!/bin/sh

find ${BIL} -readable -type f -name $1 -printf "\"%p\"\n" | xargs less -N

