#!/bin/sh

bil $@
bil -pos GmshParsed $1
gmsh $1.pos*&

