#
#
SHELL = /bin/sh



# Bil path
# --------
BIL_PATH := ${shell realpath ../..}
export BIL_PATH

include ../../make.inc


# Rules
# -----


#=======================================================================
# Default target executed when no arguments are given to make.
entry: lex yacc




#=======================================================================
# Target rules for cleaning

clean:
	rm -f Parser.yy.c Parser.tab.c Parser.tab.h ${wildcard *.o}



#=======================================================================
# Target rules for lexical analysis (lexing)

lex: Parser.yy.c

Parser.yy.c Parser.yy.cpp: Parser.l
	@echo "\nLexing: "$<
	${LEX} -P Parser_yy -o Parser.yy.c Parser.l



#=======================================================================
# Target rules for parsing

yacc: Parser.tab.c Parser.tab.h

Parser.tab.c Parser.tab.cpp Parser.tab.h: Parser.y
	@echo "\nParsing: "$<
	${YACC} -p Parser_yy -o Parser.tab.c -d Parser.y



#=======================================================================
# Target rules for exec

parse: Parser.tab.o Parser.yy.o
	${CC} -o parse Parser.tab.o Parser.yy.o -lm


#=======================================================================
# Target rules for testing

test:
	@echo '----------------------------------------------'
	@echo 'Directories'
	@echo '----------------------------------------------'
	@echo '1. All directories'
	@echo 'DIRSALL = '${DIRSALL}
	@echo '----------------------------------------------'
	@echo '2. Directories except those of Models and Modules'
	@echo 'DIRSOTHER = '${DIRSOTHER}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Sources'
	@echo '----------------------------------------------'
	@echo '1. All source files'
	@echo 'SRCALL = '${SRCALL}
	@echo '----------------------------------------------'
	@echo '2. Source Model Files'
	@echo 'SRCMODELS = '${SRCMODELS}
	@echo '----------------------------------------------'
	@echo '3. Source Module Files'
	@echo 'SRCMODULES = '${SRCMODULES}
	@echo '----------------------------------------------'
	@echo '4. Other source files .c'
	@echo 'SRCOTHERC = '${SRCOTHERC}
	@echo '----------------------------------------------'
	@echo '5. Other source files .cpp'
	@echo 'SRCOTHERCPP = '${SRCOTHERCPP}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Headers'
	@echo '----------------------------------------------'
	@echo '1. All Headers'
	@echo 'HEADERSALL = '${HEADERSALL}
	@echo '----------------------------------------------'
	@echo '2. Headers from sources .h'
	@echo 'HEADERSC = '${HEADERSC}
	@echo '----------------------------------------------'
	@echo '3. Headers from sources .hpp'
	@echo 'HEADERSCPP = '${HEADERSCPP}
	@echo '----------------------------------------------'
	@echo '4. Header files without associated source files'
	@echo 'HEADERSONLY = '${HEADERSONLY}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Objects'
	@echo '----------------------------------------------'
	@echo 'OBJALL = '${OBJALL}
	@echo '----------------------------------------------'

test1:
	@echo '----------------------------------------------'
	@echo 'Source files .cpp'
	@echo 'SRCALL = '${wildcard ${addsuffix *.c,${DIRSOTHER}}} ${wildcard ${addsuffix *.cpp,${DIRSOTHER}}}
	@echo '----------------------------------------------'
	@echo 'Object files'
	@echo 'OBJCPP from SRCOTHERC = '${SRCOTHERC:%.cpp=%.o}
	@echo 'OBJCPP from SRCOTHERC = '${filter %.o,${patsubst %.c,%.o,toto.cpp titi.c}}

test2:
	@echo 'Bil path = ' ${BIL_PATH}

