%{
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Common/Message.h"
#include "../Tools/Math.h"
#include "../DataSet/DataSet.h"
#include "Parser.tab.h"

void   parsestring(char endchar);
char*  strsave(char* ptr);
void   skipcomments(void);
void   skipline(void);

#define YY_INPUT(buf,result,max_size)					\
     {									\
       int c = '*', n;							\
       for ( n = 0; n < max_size &&					\
	       (c = getc( yyin )) != EOF && c != '\n'; ++n )		\
	 buf[n] = (char) c;						\
       if ( c == '\n' ){						\
	 buf[n++] = (char) c;						\
	 yylineno++;							\
       }								\
       if ( c == EOF && ferror( yyin ) )				\
	 Message_FatalError("Input in flex scanner failed");	        \
       result = n;							\
     }	


#define SAVEYYTEXT    \
        do { \
          Parser_yylval.c = strsave((char*)yytext) ;  \
          PRINTYYTEXT; \
          if(!Parser_yylval.c) { \
            Message_FatalError("impossible to save \"%s\"",yytext) ; \
          } \
        } while(0)


#define YYTEXTTOI \
        do { \
          Parser_yylval.u = atoi((char *)yytext); \
          PRINTYYTEXT; \
        } while(0)


#define YYTEXTTOD \
        do { \
          Parser_yylval.d = atof((char *)yytext); \
          PRINTYYTEXT; \
        } while(0)


#define PRINTYYTEXT \
        Message_Direct("%s\n",yytext)


#define PRINT(...) \
        Message_Direct(__VA_ARGS__)



%}

blank      [\ \t\n\r\f]
alpha	   [a-zA-Z\_]
special    [\.\-]
digit	   [0-9]
integer    [-+]?{digit}+
exp	   [Ee]{integer}
real       [-+]?({digit}+|{digit}+"."{digit}*|{digit}*"."{digit}+)({exp})?
string	   {alpha}({alpha}|{digit}|{special})*


%%


{blank}		          /* none */;
";"                       { SAVEYYTEXT; return tEND; }
"/*"			  skipcomments();
"//"|"#"                  skipline();
"\""			  { parsestring('\"'); return tBIGSTR; }
"\'"			  { parsestring('\''); return tBIGSTR; }
"="                       { SAVEYYTEXT; return tAFFECT; }
"+="                      { SAVEYYTEXT; return tAFFECTPLUS; }
"-="                      { SAVEYYTEXT; return tAFFECTMINUS; }
"*="                      { SAVEYYTEXT; return tAFFECTTIMES; }
"/="                      { SAVEYYTEXT; return tAFFECTDIVIDE; }
":"                       { SAVEYYTEXT; return tDOTS; }
"..."                     { SAVEYYTEXT; return tDOTS; }
"||"                      { SAVEYYTEXT; return tOR; }
"&&"                      { SAVEYYTEXT; return tAND; }
"++"                      { SAVEYYTEXT; return tPLUSPLUS; }
"--"                      { SAVEYYTEXT; return tMINUSMINUS; }
"=="                      { SAVEYYTEXT; return tEQUAL; }
"!="                      { SAVEYYTEXT; return tNOTEQUAL; }
"<="                      { SAVEYYTEXT; return tLESSOREQUAL; }
">="                      { SAVEYYTEXT; return tGREATEROREQUAL; }


Geometry                  { SAVEYYTEXT; return tGeometry; }
Dimension                 { SAVEYYTEXT; return tDimension; }
Symmetry                  { SAVEYYTEXT; return tSymmetry; }

Mesh                      { SAVEYYTEXT; return tMesh; }

Periodicity               { SAVEYYTEXT; return tPeriodicity; }
MasterRegion              { SAVEYYTEXT; return tMasterRegion; }
SlaveRegion               { SAVEYYTEXT; return tSlaveRegion; }
PeriodVector              { SAVEYYTEXT; return tPeriodVector; }

Model                     { SAVEYYTEXT; return tModel; }

Material                  { SAVEYYTEXT; return tMaterial; }
Parameters                { SAVEYYTEXT; return tParameters; }
Curves                    { SAVEYYTEXT; return tCurves; }

Field                     { SAVEYYTEXT; return tField; }
Function                  { SAVEYYTEXT; return tFunction; }
Initialization            { SAVEYYTEXT; return tInitialization; }
BoundaryCondition         { SAVEYYTEXT; return tBoundaryCondition; }
Load                      { SAVEYYTEXT; return tLoad; }
Point                     { SAVEYYTEXT; return tPoint; }
Date                      { SAVEYYTEXT; return tDate; }
ObjectiveVariation        { SAVEYYTEXT; return tObjectiveVariation; }
IterativeProcess          { SAVEYYTEXT; return tIterativeProcess; }
TimeStep                  { SAVEYYTEXT; return tTimeStep; }

Units                     { SAVEYYTEXT; return tUnits; }
Length                    { SAVEYYTEXT; return tLength; }
Time                      { SAVEYYTEXT; return tTime; }
Mass                      { SAVEYYTEXT; return tMass; }

Exit                      { SAVEYYTEXT; return tExit; }



Acos                    return tAcos;
ArcCos                  return tAcos;
Asin                    return tAsin;
ArcSin                  return tAsin;
Atan                    return tAtan;
ArcTan                  return tAtan;
Atan2                   return tAtan2;
ArcTan2                 return tAtan2;
Cosh                    return tCosh;
Cos                     return tCos;
Exp                     return tExp;
Fabs                    return tFabs;
Log                     return tLog;
Log10                   return tLog10;
Pi                      return tPi;
Sqrt                    return tSqrt;
Sin                     return tSin;
Sinh                    return tSinh;
Tanh                    return tTanh;
Tan                     return tTan;


{integer}                 { YYTEXTTOI; return tINT; }

{real}                    { YYTEXTTOD; return tDOUBLE; }

{string}		  { SAVEYYTEXT; return tSTRING; }

.                         { SAVEYYTEXT; return yytext[0]; }

%%

#undef Parser_yywrap

int Parser_yywrap() {return 1;}

int yyinput(void) {return input();}



void skipcomments(void)
{
  int c;

  while (1) {
    while ((c = yyinput()) != '*'){
      if(feof(yyin)){
	Message_FatalError("End of file in commented region");
        return;
      }
    }
    if ((c = yyinput()) == '/')
      return;
    unput(c);
  }
}



void parsestring(char endchar)
{
  int c;
  char tmp[1024];
  
  // Note that we keep special characters (end-of-line \n, tabs \t,
  // etc.) "as is" in the output string: see yyinput() above
  int i = 0;
  while ((c = yyinput()) != endchar) {
    if(feof(yyin)){
      Message_FatalError("End of file in string");
      break;
    }
    else if(i >= (int)sizeof(tmp)-1){
      Message_FatalError("String too long");
      break;
    }
    else{
      tmp[i++] = (char)c;
    }
  }
  tmp[i] = '\0';
  Parser_yylval.c = strsave(tmp);
}



char* strsave(char* ptr)
{
  char* d = (char*) malloc(strlen(ptr)+1) ;

  if(d == NULL) return NULL ;

  strcpy(d,ptr);

  return(d);
}



void skipline(void)
{
  int c;
  while ((c = yyinput()) != '\n'){
    if(feof(yyin)) return;
  }
}



void skip_until(const char* skip, const char* until)
{
  int i, nb_skip;
  int l, l_skip, l_until;
  char chars[256];

  nb_skip = 0;

  if(skip)
    l_skip = strlen(skip);
  else
    l_skip = 0;

  l_until = strlen(until);

  while(1){
    while (1){
      chars[0] = yyinput();
      if(feof(yyin)){
	Message_FatalError("Unexpected end of file");
	return;
      }
      if(chars[0] == until[0]) break;
      if(skip && chars[0] == skip[0]) break;
    }

    l = MAX(l_skip,l_until);
    if(l >= (int)sizeof(chars)){
      Message_FatalError("Search pattern too long in skip_until");
      return;
    }
    for(i=1; i<l; i++){
      chars[i] = yyinput();
      if(feof(yyin)){
	l = i;
	break;
      }
    }

    if(!strncmp(chars,until,l_until)){
      if(!nb_skip){
	return;
      }
      else{
	nb_skip--;
      }
    }
    else if(skip && !strncmp(chars,skip,l_skip)){
      nb_skip++;
    }
    else{
      for(i=1;i<l-1;i++){
	unput(chars[l-i]);
      }
    }

  }
}



void Parser_yyflush() { YY_FLUSH_BUFFER; }
