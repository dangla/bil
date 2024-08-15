#!/bin/sh

find ${BIL}/src -maxdepth 3 -name "*.*[chpf]" | xargs wc -l
