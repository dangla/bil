#
SHELL = /bin/sh

showtitles:
	grep --include="*.c" --include="*.cpp" -r -T --color -e "#define TITLE" .


showauthors:
	grep --include="*.c" --include="*.cpp" -r -T --color -e "#define AUTHORS" .


awktitles:
	awk '$$2 == "TITLE" {gsub("#define TITLE","",$$0); gsub("\"","",$$0); printf("%s %s\n",FILENAME,$$0)}' */*.c
