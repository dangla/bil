#
#

dvi:
	latex bil

ps: dvi
	dvips bil -o

pdf: ps
	ps2pdf bil.ps bil.pdf 

all: ps pdf

clean:
	rm -f bil.dvi bil.log bil.aux bil.ps bil.pdf

