%{
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "../Common/Message.h"
#include "../Common/Mry.h"
#include "../DataSet/Geometry.h"
#include "../DataSet/DataFile.h"
#include "../DataSet/DataSet.h"
#include "../Common/InternationalSystemOfUnits.h"
#include "Parser.h"
#include "../Libraries/map.h"

typedef map_t(Parser_yysymbols_t) map_yysymbols_t ;

#define Free  Mry_Free

// Global parser variables

char* Parser_yyname;
int Parser_yyerrorstate = 0;
int Parser_yyviewindex = 0;
map_double_t Parser_yysymbols;

// Static parser variables (accessible only in this file)

void yyerror(DataSet_t*,const char*);
void yymsg(int level, const char* fmt, ...);
void skip_until(const char* skip, const char* until);
//int PrintListOfDouble(char *format, List_T *list, char *buffer);
//fullMatrix<double> ListOfListOfDouble2Matrix(List_T *list);


%}

%union {
  char* c;
  int i;
  unsigned int u;
  double d;
  double v[3];
}

%parse-param {DataSet_t* dataset}

%token <d> tDOUBLE
%token <i> tINT
%token <c> tSTRING tBIGSTR

%token <c> tEND tAFFECT tDOTS

%token <d> tPi

%token tExp tLog tLog10 tSqrt tSin tAsin tCos tAcos tTan tRand
%token tAtan tAtan2 tSinh tCosh tTanh tFabs tFloor tCeil
%token tFmod tModulo tHypot 
%token tExit


%token <c> tGeometry tDimension tSymmetry
%token <c> tMesh
%token <c> tPeriodicity tMasterRegion tSlaveRegion tPeriodVector
%token <c> tModel
%token <c> tMaterial tParameters tCurves tCurves_log
%token <c> tField tType tValue tGradient tFile
%token <c> tFunction
%token <c> tInitialization tRegion tUnknown
%token <c> tBoundaryCondition
%token <c> tLoad tEquation
%token <c> tPoint tCoordinates
%token <c> tDate
%token <c> tObjectiveVariation
%token <c> tIterativeProcess tIterations tTolerance tRepetitions
%token <c> tTimeStep tDtini tDtmax tDtmin
%token <c> tUnits tLength tTime tMass


%type <c> BilFormatItem

%type <d> FExpr FExpr_Single 
%type <v> VExpr VExpr_Single
%type <i> NumericAffectation NumericIncrement


// Operators (with ascending priority): cf. C language
//
// Notes: - associativity (%left, %right)
//        - UNARYPREC is a dummy terminal to resolve ambiguous cases
//          for + and - (which exist in both unary and binary form)

%right   tAFFECT tAFFECTPLUS tAFFECTMINUS tAFFECTTIMES tAFFECTDIVIDE
%right   '?' tDOTS
%left    tOR
%left    tAND
%left    tEQUAL tNOTEQUAL
%left    '<' tLESSOREQUAL  '>' tGREATEROREQUAL
%left    '+' '-'
%left    '*' '/' '%'
%right   '!' tPLUSPLUS tMINUSMINUS UNARYPREC
%right   '^'
%left    '(' ')' '[' ']' '.' '#'


%start All

%%

All : 
    BilFormatItems
  | error tEND { yyerrok; return 1; }
;



//  B I L   F I L E   F O R M A T

BilFormatItems : 
    /* empty */
  | BilFormatItems BilFormatItem
    {
      Message_Direct("%s has been read\n",$2) ;
    }
;



BilFormatItem :
    Geometry
  | Mesh
  | Periodicity
  | tModel
  | tMaterial
  | tField
  | tFunction
  | tInitialization
  | tBoundaryCondition
  | tPoint
  | tDate
  | tObjectiveVariation
  | tIterativeProcess
  | tTimeStep
  | Units
  | Command
;



// U N I T S

Units :
    tUnits tAFFECT '{' UnitsItems '}' tEND
  | UnitsItems
;



UnitsItems :
    UnitsItem
  | UnitsItems UnitsItem
;



UnitsItem :
    tLength tAFFECT tSTRING tEND 
    {
      InternationalSystemOfUnits_UseAsLength($3) ;
      printf("%s = { %s }\n",$1,$3);
    }
  | tTime  tAFFECT tSTRING tEND
    {
      InternationalSystemOfUnits_UseAsTime($3) ;
      printf("%s = { %s }\n",$1,$3);
    }
  | tMass  tAFFECT tSTRING tEND
    {
      InternationalSystemOfUnits_UseAsMass($3) ;
      printf("%s = { %s }\n",$1,$3);
    }
;



// G E O M E T R Y

Geometry :
    tGeometry tAFFECT '{' GeometryItems '}' tEND
  | GeometryItems
;



GeometryItems :
    GeometryItem
  | GeometryItems GeometryItem
;



GeometryItem :
    tDimension tAFFECT tINT tEND 
    {
      Geometry_t* geom = DataSet_GetGeometry(dataset) ;
      
      if(!geom) {
        geom = Geometry_New() ;
        DataSet_GetGeometry(dataset) = geom ;
      }

      Geometry_GetDimension(geom) = $3 ;

      printf("%s = { %d } %s\n",$1,$3,$4);
    }
  | tSymmetry  tAFFECT tSTRING tEND
    {
      Geometry_t* geom = DataSet_GetGeometry(dataset) ;
      
      if(!geom) {
        geom = Geometry_New() ;
        DataSet_GetGeometry(dataset) = geom ;
      }

      {
        if(!strncasecmp($3,"plan",4))  {
          Geometry_SetPlaneSymmetry(geom) ;
        } else if(!strncasecmp($3,"axis",4))  {
          Geometry_SetCylindricalSymmetry(geom) ;
        } else if(!strncasecmp($3,"sphe",4))  {
          Geometry_SetSphericalSymmetry(geom) ;
        } else {
          Message_FatalError("Geometry_Create: geometry not available\n\
          Available geometries are: PLAN, AXIS, SPHE") ;
        }
      }

      printf("%s = { %s } %s\n",$1,$3,$4);
    }
;



// M E S H

Mesh :
    tMesh tAFFECT tBIGSTR tEND
    {
      Mesh_t* mesh = (Mesh_t*) Mry_New(Mesh_t) ;

      DataSet_GetMesh(dataset) = mesh ;
      Mesh_GetGeometry(mesh) = DataSet_GetGeometry(dataset) ;
      Mesh_Scan(mesh,$3) ;
      Mesh_CreateMore(mesh) ;
      {
        Materials_t* materials = DataSet_GetMaterials(dataset) ;
        if(materials) {
          Elements_LinkUp(Mesh_GetElements(mesh),materials) ;
          Elements_CreateMore(Mesh_GetElements(mesh)) ;
          Nodes_CreateMore(Mesh_GetNodes(mesh)) ;
          Mesh_SetEquationContinuity(mesh) ;
        }
      }
    }
;



// P E R I O D I C I T Y

Periodicity :
    tPeriodicity '[' tINT ']' tAFFECT '{' PeriodicityItems '}' tEND
;



PeriodicityItems :
    PeriodicityItem
  | PeriodicityItems PeriodicityItem
;



PeriodicityItem :
    tMasterRegion tAFFECT FExpr tEND 
    {
      printf("%s = { %d }\n",$1,$3);
    }
  | tSlaveRegion  tAFFECT FExpr tEND 
  | tPeriodVector tAFFECT VExpr tEND 
;



//  C O M M A N D  

Command :
     tExit tEND
    {
      YYABORT;
    }
;




//  A F F E C T A T I O N

NumericAffectation :
    tAFFECT        { $$ = 0; }
  | tAFFECTPLUS    { $$ = 1; }
  | tAFFECTMINUS   { $$ = 2; }
  | tAFFECTTIMES   { $$ = 3; }
  | tAFFECTDIVIDE  { $$ = 4; }
;

NumericIncrement :
    tPLUSPLUS      { $$ = 1; }
  | tMINUSMINUS    { $$ = -1; }
;


Affectation :

  // Variables

    tSTRING NumericAffectation FExpr tEND
    {
      double* v = map_get(&Parser_yysymbols,$1);

      if(!v){
        if(map_set(&Parser_yysymbols,$1,0.)) {
	  yymsg(0, "impossible to set '%s'", $1);
        }
        v = map_get(&Parser_yysymbols,$1);
      }

      {
	switch($2){
	case 0 : v[0]  = $3; break;
	case 1 : v[0] += $3; break;
	case 2 : v[0] -= $3; break;
	case 3 : v[0] *= $3; break;
	case 4 : 
	  if($3) {
            v[0] /= $3; 
	  } else {
            yymsg(0, "Division by zero in '%s /= %g'", $1, $3);
          }
	  break;
	}
      }

      Free($1);
    }
;









//  G E N E R A L

FExpr :
    FExpr_Single                     { $$ = $1;           }
  | '(' FExpr ')'                    { $$ = $2;           }
  | '-' FExpr %prec UNARYPREC        { $$ = -$2;          }
  | '+' FExpr %prec UNARYPREC        { $$ = $2;           }
  | '!' FExpr                        { $$ = !$2;          }
  | FExpr '-' FExpr                  { $$ = $1 - $3;      }
  | FExpr '+' FExpr                  { $$ = $1 + $3;      }
  | FExpr '*' FExpr                  { $$ = $1 * $3;      }
  | FExpr '/' FExpr
    { 
      if(!$3) {
	yymsg(0, "Division by zero in '%g / %g'", $1, $3);
      } else {
	$$ = $1 / $3;
      }
    }
  | FExpr '%' FExpr                  { $$ = (int)$1 % (int)$3;  }
  | FExpr '^' FExpr                  { $$ = pow($1, $3);  }
  | FExpr '<' FExpr                  { $$ = $1 < $3;      }
  | FExpr '>' FExpr                  { $$ = $1 > $3;      }
  | FExpr tLESSOREQUAL FExpr         { $$ = $1 <= $3;     }
  | FExpr tGREATEROREQUAL FExpr      { $$ = $1 >= $3;     }
  | FExpr tEQUAL FExpr               { $$ = $1 == $3;     }
  | FExpr tNOTEQUAL FExpr            { $$ = $1 != $3;     }
  | FExpr tAND FExpr                 { $$ = $1 && $3;     }
  | FExpr tOR FExpr                  { $$ = $1 || $3;     }
  | FExpr '?' FExpr tDOTS FExpr      { $$ = $1 ? $3 : $5; }
  | tExp    '(' FExpr ')'            { $$ = exp($3);      }
  | tLog    '(' FExpr ')'            { $$ = log($3);      }
  | tLog10  '(' FExpr ')'            { $$ = log10($3);    }
  | tSqrt   '(' FExpr ')'            { $$ = sqrt($3);     }
  | tSin    '(' FExpr ')'            { $$ = sin($3);      }
  | tAsin   '(' FExpr ')'            { $$ = asin($3);     }
  | tCos    '(' FExpr ')'            { $$ = cos($3);      }
  | tAcos   '(' FExpr ')'            { $$ = acos($3);     }
  | tTan    '(' FExpr ')'            { $$ = tan($3);      }
  | tAtan   '(' FExpr ')'            { $$ = atan($3);     }
  | tAtan2  '(' FExpr ',' FExpr ')'  { $$ = atan2($3, $5);}
  | tSinh   '(' FExpr ')'            { $$ = sinh($3);     }
  | tCosh   '(' FExpr ')'            { $$ = cosh($3);     }
  | tTanh   '(' FExpr ')'            { $$ = tanh($3);     }
  | tFabs   '(' FExpr ')'            { $$ = fabs($3);     }
  | tFloor  '(' FExpr ')'            { $$ = floor($3);    }
  | tCeil   '(' FExpr ')'            { $$ = ceil($3);     }
  | tFmod   '(' FExpr ',' FExpr ')'  { $$ = fmod($3, $5); }
  | tModulo '(' FExpr ',' FExpr ')'  { $$ = fmod($3, $5); }
  | tHypot  '(' FExpr ',' FExpr ')'  { $$ = sqrt($3 * $3 + $5 * $5); }
  | tRand   '(' FExpr ')'            { $$ = $3 * (double)rand() / (double)RAND_MAX; }
;



FExpr_Single :

  // Constants

    tINT      { $$ = $1; }
  | tDOUBLE   { $$ = $1; }
  | tPi       { $$ = 3.141592653589793; }

  // Variables

  | tSTRING
    {
      double* v = map_get(&Parser_yysymbols,$1);

      if(!v){
	yymsg(0, "Unknown variable '%s'", $1);
	$$ = 0.;
      } else {
        $$ = v[0];
      }

      Free($1);
    }
  | tSTRING '[' FExpr ']'
    {
      int index = (int) $3;
      double* v = map_get(&Parser_yysymbols,$1);

      if(!v){
	yymsg(0, "Unknown variable '%s'", $1);
	$$ = 0.;
      } else {
	$$ = v[index];
      }

      Free($1);
    }
;



VExpr :
    VExpr_Single
    {
      memcpy($$, $1, 3*sizeof(double));
    }
  | '-' VExpr %prec UNARYPREC
    {
      for(int i = 0; i < 3; i++) $$[i] = -$2[i];
    }
  | '+' VExpr %prec UNARYPREC
    { 
      for(int i = 0; i < 3; i++) $$[i] = $2[i];
    }
  | VExpr '-' VExpr
    { 
      for(int i = 0; i < 3; i++) $$[i] = $1[i] - $3[i];
    }
  | VExpr '+' VExpr
    {
      for(int i = 0; i < 3; i++) $$[i] = $1[i] + $3[i];
    }
;



VExpr_Single :
    '{' FExpr ',' FExpr ',' FExpr '}'
    { 
      $$[0] = $2;  $$[1] = $4;  $$[2] = $6;
    }
  | '(' FExpr ',' FExpr ',' FExpr ')'
    {
      $$[0] = $2;  $$[1] = $4;  $$[2] = $6;
    }
;





%%

/*
int PrintListOfDouble(char *format, List_T *list, char *buffer)
{
  int j, k;
  char tmp1[256], tmp2[256];

  j = 0;
  buffer[j] = '\0';

  while(j < (int)strlen(format) && format[j] != '%') j++;
  strncpy(buffer, format, j); 
  buffer[j]='\0'; 
  for(int i = 0; i < List_Nbr(list); i++){
    k = j;
    j++;
    if(j < (int)strlen(format)){
      if(format[j] == '%'){
	strcat(buffer, "%");
	j++;
      }
      while(j < (int)strlen(format) && format[j] != '%') j++;
      if(k != j){
	strncpy(tmp1, &(format[k]), j-k);
	tmp1[j-k] = '\0';
	sprintf(tmp2, tmp1, *(double*)List_Pointer(list, i)); 
	strcat(buffer, tmp2);
      }
    }
    else
      return List_Nbr(list)-i;
  }
  if(j != (int)strlen(format))
    return -1;
  return 0;
}



fullMatrix<double> ListOfListOfDouble2Matrix(List_T *list)
{
  int M = List_Nbr(list);
  int N = 0;
  for(int i = 0; i < M; i++){
    List_T *line = *(List_T**)List_Pointer_Fast(list, i);
    N = std::max(N, List_Nbr(line));
  }
  fullMatrix<double> mat(M, N);
  for(int i = 0; i < M; i++){
    List_T *line = *(List_T**)List_Pointer_Fast(list, i);
    for(int j = 0; j < List_Nbr(line); j++){
      double val;
      List_Read(line, j, &val);
      mat(i, j) = val;
    }
  }
  for(int i = 0; i < List_Nbr(list); i++)
    List_Delete(*(List_T**)List_Pointer(list, i));
  List_Delete(list);
  return mat;
}
*/


void yyerror(DataSet_t* dataset,const char* s)
{
  Message_FatalError("'%s', line %d : %s (%s)", Parser_yyname, Parser_yylineno - 1,s, Parser_yytext);
  Parser_yyerrorstate++;
}



void yymsg(int level, const char* fmt, ...)
{
  va_list args;
  char tmp[1024];

  va_start(args, fmt);
  vsprintf(tmp, fmt, args);
  va_end(args);

  if(level == 0){
    Message_FatalError("'%s', line %d : %s", Parser_yyname, Parser_yylineno - 1, tmp);
    Parser_yyerrorstate++;
  }
  else
    Message_Warning("'%s', line %d : %s", Parser_yyname, Parser_yylineno - 1, tmp);
}



Parser_t* (Parser_Create)(void)
{
  Parser_t* parser = (Parser_t*) malloc(sizeof(Parser_t)) ;
  
  assert(parser) ;
  
  return(parser) ;
}



void (Parser_Delete)(void* self)
{
  Parser_t** pparser = (Parser_t**) self ;
  Parser_t*  parser  = *pparser ;
  
  free(parser) ;
  *pparser = NULL ;
}



int Parser_ParseFile(DataSet_t* dataset)
{
  DataFile_t* datafile = DataSet_GetDataFile(dataset) ;
  // add 'b' for pure Windows programs: opening in text mode messes up
  // fsetpos/fgetpos (used e.g. for user-defined functions)

  Parser_yyname       = DataFile_GetFileName(datafile);
  Parser_yyin         = DataFile_OpenFile(datafile,"rb");
  Parser_yyerrorstate = 0;
  Parser_yylineno     = 1;
  Parser_yyviewindex  = 0;

  map_init(&Parser_yysymbols) ;

  while(!feof(Parser_yyin)) {
    if(Parser_yyparse(dataset)) break ;
    //Parser_yylex();
    
    if(Parser_yyerrorstate > 20){
      Message_FatalError("Too many errors: aborting...");
      Parser_yyflush();
      break;
    }
  }

  map_deinit(&Parser_yysymbols);

  fclose(Parser_yyin);

  return 1;
}





DataSet_t* Parser_CreateDataSet(char* filename,Options_t* options)
{
  DataSet_t*  dataset  = DataSet_New() ;
  
  DataSet_GetOptions(dataset) = options ;

  {
    DataFile_t* datafile = DataFile_Create(filename) ;
  
    DataSet_GetDataFile(dataset) = datafile ;
  
    if(DataFile_DoesNotExist(datafile)) {
      Message_FatalError("The file \"%s\" doesn't exist\n",filename) ;
      exit(EXIT_SUCCESS) ;
    }
  }

  Parser_ParseFile(dataset) ;

  return(dataset) ;
}

