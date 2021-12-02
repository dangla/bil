#!/bin/sh
# Use [[:blank:]] as a regular expression for a blank.

valgrind --tool=memcheck --leak-check=full --read-var-info=yes bil $1 $2 $3 $4
