#ifndef PARSER_H
#define PARSER_H


/* vacuous declarations and typedef names */
struct Parser_s  ; typedef struct Parser_s Parser_t ;



//#include <map>
//#include <string>
//#include <vector>
#include "DataSet.h"
#include "Options.h"

extern Parser_t*  (Parser_Create)(void) ;
extern void       (Parser_Delete)(void*) ;
extern int        (Parser_ParseFile)(DataSet_t*) ;
extern DataSet_t* (Parser_CreateDataSet)(char*,Options_t*) ;



union Parser_yysymbols_u {
  char* c;
  int   i;
  double d;
  double v[3];
} ;

typedef union Parser_yysymbols_u Parser_yysymbols_t ;




#define Parser_Getyyin(PA)          ((PA)->yyin)
#define Parser_Getyylineno(PA)      ((PA)->yylineno)
#define Parser_Getyytext(PA)        ((PA)->yytext)
#define Parser_Getyyviewindex(PA)   ((PA)->yyviewindex)
#define Parser_Getyyname(PA)        ((PA)->yyname)
#define Parser_Getyyerrorstate(PA)  ((PA)->yyerrorstate)


int  Parser_yyparse(DataSet_t*);
int  Parser_yylex();
void Parser_yyflush();

// global parser variables that need to be exported
extern FILE* Parser_yyin;
extern int   Parser_yylineno;
extern char* Parser_yytext;
extern int   Parser_yyviewindex;
extern char* Parser_yyname;
extern int   Parser_yyerrorstate;
//extern std::map<std::string, std::vector<double> > Parser_yysymbols;



struct Parser_s {      /* Parser */
  FILE* yyin;
  int   yylineno;
  char* yytext;
  int   yyviewindex;
  char* yyname;
  int   yyerrorstate;
} ;


#endif
