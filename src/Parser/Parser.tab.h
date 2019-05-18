/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_PARSER_YY_PARSER_TAB_H_INCLUDED
# define YY_PARSER_YY_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int Parser_yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    tDOUBLE = 258,
    tINT = 259,
    tSTRING = 260,
    tBIGSTR = 261,
    tEND = 262,
    tAFFECT = 263,
    tDOTS = 264,
    tPi = 265,
    tExp = 266,
    tLog = 267,
    tLog10 = 268,
    tSqrt = 269,
    tSin = 270,
    tAsin = 271,
    tCos = 272,
    tAcos = 273,
    tTan = 274,
    tRand = 275,
    tAtan = 276,
    tAtan2 = 277,
    tSinh = 278,
    tCosh = 279,
    tTanh = 280,
    tFabs = 281,
    tFloor = 282,
    tCeil = 283,
    tFmod = 284,
    tModulo = 285,
    tHypot = 286,
    tExit = 287,
    tGeometry = 288,
    tDimension = 289,
    tSymmetry = 290,
    tMesh = 291,
    tPeriodicity = 292,
    tMasterRegion = 293,
    tSlaveRegion = 294,
    tPeriodVector = 295,
    tModel = 296,
    tMaterial = 297,
    tParameters = 298,
    tCurves = 299,
    tCurves_log = 300,
    tField = 301,
    tType = 302,
    tValue = 303,
    tGradient = 304,
    tFile = 305,
    tFunction = 306,
    tInitialization = 307,
    tRegion = 308,
    tUnknown = 309,
    tBoundaryCondition = 310,
    tLoad = 311,
    tEquation = 312,
    tPoint = 313,
    tCoordinates = 314,
    tDate = 315,
    tObjectiveVariation = 316,
    tIterativeProcess = 317,
    tIterations = 318,
    tTolerance = 319,
    tRepetitions = 320,
    tTimeStep = 321,
    tDtini = 322,
    tDtmax = 323,
    tDtmin = 324,
    tUnits = 325,
    tLength = 326,
    tTime = 327,
    tMass = 328,
    tAFFECTPLUS = 329,
    tAFFECTMINUS = 330,
    tAFFECTTIMES = 331,
    tAFFECTDIVIDE = 332,
    tOR = 333,
    tAND = 334,
    tEQUAL = 335,
    tNOTEQUAL = 336,
    tLESSOREQUAL = 337,
    tGREATEROREQUAL = 338,
    tPLUSPLUS = 339,
    tMINUSMINUS = 340,
    UNARYPREC = 341
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 39 "Parser.y" /* yacc.c:1909  */

  char* c;
  int i;
  unsigned int u;
  double d;
  double v[3];

#line 149 "Parser.tab.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE Parser_yylval;

int Parser_yyparse (DataSet_t* dataset);

#endif /* !YY_PARSER_YY_PARSER_TAB_H_INCLUDED  */
