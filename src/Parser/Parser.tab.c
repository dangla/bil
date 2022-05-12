/* A Bison parser, made by GNU Bison 3.5.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
   Inc.

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.5.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         Parser_yyparse
#define yylex           Parser_yylex
#define yyerror         Parser_yyerror
#define yydebug         Parser_yydebug
#define yynerrs         Parser_yynerrs
#define yylval          Parser_yylval
#define yychar          Parser_yychar

/* First part of user prologue.  */
#line 1 "Parser.y"

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



#line 115 "Parser.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
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
union YYSTYPE
{
#line 39 "Parser.y"

  char* c;
  int i;
  unsigned int u;
  double d;
  double v[3];

#line 262 "Parser.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE Parser_yylval;

int Parser_yyparse (DataSet_t* dataset);

#endif /* !YY_PARSER_YY_PARSER_TAB_H_INCLUDED  */



#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))

/* Stored state numbers (used for stacks). */
typedef yytype_uint8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  5
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1128

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  106
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  95
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  254

#define YYUNDEFTOK  2
#define YYMAXUTOK   341


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    92,     2,   102,     2,    91,     2,     2,
      97,    98,    89,    87,   105,    88,   101,    90,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      83,     2,    85,    78,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    99,     2,   100,    96,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   103,     2,   104,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    79,    80,    81,    82,    84,    86,    93,
      94,    95
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   112,   112,   113,   120,   122,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   154,   155,   161,   162,   168,   173,   178,   190,
     191,   197,   198,   204,   217,   248,   272,   278,   279,   285,
     289,   290,   298,   369,   370,   371,   372,   373,   374,   375,
     376,   377,   385,   386,   387,   388,   389,   390,   391,   392,
     393,   394,   395,   396,   397,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,   412,
     413,   414,   415,   416,   425,   426,   427,   431,   444,   463,
     467,   471,   475,   479,   488,   492
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "tDOUBLE", "tINT", "tSTRING", "tBIGSTR",
  "tEND", "tAFFECT", "tDOTS", "tPi", "tExp", "tLog", "tLog10", "tSqrt",
  "tSin", "tAsin", "tCos", "tAcos", "tTan", "tRand", "tAtan", "tAtan2",
  "tSinh", "tCosh", "tTanh", "tFabs", "tFloor", "tCeil", "tFmod",
  "tModulo", "tHypot", "tExit", "tGeometry", "tDimension", "tSymmetry",
  "tMesh", "tPeriodicity", "tMasterRegion", "tSlaveRegion",
  "tPeriodVector", "tModel", "tMaterial", "tParameters", "tCurves",
  "tCurves_log", "tField", "tType", "tValue", "tGradient", "tFile",
  "tFunction", "tInitialization", "tRegion", "tUnknown",
  "tBoundaryCondition", "tLoad", "tEquation", "tPoint", "tCoordinates",
  "tDate", "tObjectiveVariation", "tIterativeProcess", "tIterations",
  "tTolerance", "tRepetitions", "tTimeStep", "tDtini", "tDtmax", "tDtmin",
  "tUnits", "tLength", "tTime", "tMass", "tAFFECTPLUS", "tAFFECTMINUS",
  "tAFFECTTIMES", "tAFFECTDIVIDE", "'?'", "tOR", "tAND", "tEQUAL",
  "tNOTEQUAL", "'<'", "tLESSOREQUAL", "'>'", "tGREATEROREQUAL", "'+'",
  "'-'", "'*'", "'/'", "'%'", "'!'", "tPLUSPLUS", "tMINUSMINUS",
  "UNARYPREC", "'^'", "'('", "')'", "'['", "']'", "'.'", "'#'", "'{'",
  "'}'", "','", "$accept", "All", "BilFormatItems", "BilFormatItem",
  "Units", "UnitsItems", "UnitsItem", "Geometry", "GeometryItems",
  "GeometryItem", "Mesh", "Periodicity", "PeriodicityItems",
  "PeriodicityItem", "Command", "FExpr", "FExpr_Single", "VExpr",
  "VExpr_Single", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,    63,   333,
     334,   335,   336,    60,   337,    62,   338,    43,    45,    42,
      47,    37,    33,   339,   340,   341,    94,    40,    41,    91,
      93,    46,    35,   123,   125,    44
};
# endif

#define YYPACT_NINF (-106)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-5)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     218,     2,    19,   260,  -106,  -106,    20,    26,    33,    34,
      61,    -8,  -106,  -106,  -106,  -106,  -106,  -106,  -106,  -106,
    -106,  -106,  -106,   122,   133,   136,   137,  -106,  -106,   -47,
    -106,  -106,   -30,  -106,  -106,  -106,  -106,  -106,    12,   142,
     143,   144,   145,    48,   147,   148,   154,  -106,  -106,   -30,
     153,   155,   156,    64,   -47,   160,   161,   164,   -14,  -106,
    -106,  -106,   165,   -36,  -106,  -106,  -106,   168,    58,   169,
    -106,   100,  -106,   170,   172,   173,   -37,  -106,    82,    82,
     -75,   175,  -106,  -106,  -106,    73,  -106,    80,    86,    87,
     102,   103,   104,   105,   123,   124,   127,   128,   129,   146,
     150,   151,   152,   159,   166,   171,   174,   177,    82,    82,
      82,    82,   107,  -106,   126,   -75,   -75,    82,    82,     1,
    -106,  -106,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,   162,   162,   162,   502,  -106,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,  -106,  -106,  -106,    38,   256,  -106,
     -75,   -75,   479,   523,   544,   565,   586,   607,   628,   649,
     670,   691,   712,   733,   284,   754,   775,   796,   817,   838,
     859,   312,   340,   368,  -106,   149,  1002,  1018,  1032,  1032,
     -73,   -73,   -73,   -73,    46,    46,   162,   162,   162,   162,
      82,    82,  -106,  -106,  -106,  -106,  -106,  -106,  -106,  -106,
    -106,  -106,  -106,  -106,  -106,  -106,    82,  -106,  -106,  -106,
    -106,  -106,  -106,    82,    82,    82,    82,   396,   424,   880,
     901,   922,   943,   985,    82,    82,  -106,  -106,  -106,  -106,
     964,   452,  -106,  -106
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     2,     3,     1,     0,     0,     0,     0,
       0,     0,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,     0,     0,     0,     0,     5,    20,    23,
      24,     6,    30,    31,     7,     8,    21,    42,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    25,    32,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    33,
      34,    35,     0,     0,    26,    27,    28,     0,     0,     0,
      29,     0,    22,     0,     0,     0,     0,    37,     0,     0,
       0,     0,    38,    85,    84,    87,    86,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    43,     0,     0,     0,     0,     0,     0,
      89,    36,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    46,    45,    47,     0,    39,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    40,    91,    90,     0,     0,    41,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    44,     0,    61,    60,    58,    59,
      54,    56,    55,    57,    49,    48,    50,    51,    52,    53,
       0,     0,    93,    92,    88,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    83,    72,     0,    74,    75,    76,
      77,    78,    79,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    62,     0,     0,    73,    80,    81,    82,
       0,     0,    95,    94
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -106,  -106,  -106,  -106,  -106,   187,   -23,  -106,   193,   -25,
    -106,  -106,  -106,   181,  -106,   -79,  -106,  -105,  -106
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     2,     3,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    76,    77,    36,   112,   113,   119,   120
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     114,    73,    74,    75,     8,     9,    47,    48,   169,     4,
     165,   166,   115,   116,   158,   159,   160,   161,   162,     5,
       8,     9,   117,   163,    24,    25,    26,    37,   118,   144,
     145,   146,   147,    48,    38,    24,    25,    26,   167,   168,
      47,    39,    40,   172,   173,   174,   175,   176,   177,   178,
     179,   180,   181,   182,   183,   184,   185,   186,   187,   188,
     189,   190,   191,   192,   193,   212,   213,    81,    69,    41,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,    83,    84,    85,   170,   171,
      67,    42,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   148,    49,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
      43,   237,   238,   164,   163,   160,   161,   162,    73,    74,
      75,    44,   163,   210,    45,    46,    50,   239,    51,    53,
      52,    54,    55,    56,   240,   241,   242,   243,   236,    57,
      59,    71,    60,    61,    62,   250,   251,    64,    65,   108,
     109,    66,   122,    68,   110,    70,    72,   123,    78,   111,
      79,    80,   121,   124,   125,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,   162,   126,
     127,   128,   129,   163,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,    -4,     1,
     130,   131,   163,     0,   132,   133,   134,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,    63,    58,   135,     0,   163,     0,   136,   137,   138,
      -4,    -4,    -4,    -4,    -4,    -4,   139,    82,   163,    -4,
      -4,     0,     0,   140,    -4,     0,     0,     0,   141,    -4,
      -4,   142,     0,    -4,   143,     0,    -4,     0,    -4,    -4,
      -4,     0,     0,     0,    -4,     0,     0,     0,    -4,    -4,
      -4,    -4,     6,     7,     8,     9,    10,    11,     0,     0,
       0,    12,    13,     0,     0,     0,    14,     0,     0,     0,
       0,    15,    16,     0,     0,    17,     0,     0,    18,     0,
      19,    20,    21,     0,     0,     0,    22,     0,     0,     0,
      23,    24,    25,    26,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,     0,     0,
       0,     0,   163,     0,     0,     0,     0,     0,     0,     0,
       0,   211,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,     0,     0,     0,     0,
     163,     0,     0,     0,     0,     0,     0,     0,     0,   226,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,     0,
       0,     0,     0,     0,     0,     0,     0,   233,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
     161,   162,     0,     0,     0,     0,   163,     0,     0,     0,
       0,     0,     0,     0,     0,   234,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
       0,     0,     0,     0,   163,     0,     0,     0,     0,     0,
       0,     0,     0,   235,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,     0,     0,
       0,     0,   163,     0,     0,     0,     0,     0,     0,     0,
       0,   244,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,     0,     0,     0,     0,
     163,     0,     0,     0,     0,     0,     0,     0,     0,   245,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,     0,
       0,     0,     0,     0,     0,     0,   253,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,     0,     0,     0,     0,   163,     0,     0,     0,   214,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,     0,
     194,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,     0,     0,     0,     0,   163,
       0,   215,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,     0,     0,     0,     0,
     163,     0,   216,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   161,   162,     0,     0,     0,
       0,   163,     0,   217,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,     0,     0,
       0,     0,   163,     0,   218,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,   162,     0,
       0,     0,     0,   163,     0,   219,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
       0,     0,     0,     0,   163,     0,   220,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,     0,     0,     0,     0,   163,     0,   221,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
     161,   162,     0,     0,     0,     0,   163,     0,   222,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,     0,     0,     0,     0,   163,     0,   223,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,     0,
     224,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,     0,     0,     0,     0,   163,
       0,   225,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,     0,     0,     0,     0,
     163,     0,   227,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   161,   162,     0,     0,     0,
       0,   163,     0,   228,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,     0,     0,
       0,     0,   163,     0,   229,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,   162,     0,
       0,     0,     0,   163,     0,   230,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
       0,     0,     0,     0,   163,     0,   231,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,     0,     0,     0,     0,   163,     0,   232,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   159,   160,
     161,   162,     0,     0,     0,     0,   163,     0,   246,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,     0,     0,     0,     0,   163,     0,   247,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,     0,
     248,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,     0,     0,     0,     0,   163,
       0,   249,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,     0,     0,     0,     0,
     163,     0,   252,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   161,   162,     0,     0,     0,
       0,   163,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
       0,     0,     0,     0,   163,   154,   155,   156,   157,   158,
     159,   160,   161,   162,     0,     0,     0,     0,   163
};

static const yytype_int16 yycheck[] =
{
      79,    38,    39,    40,    34,    35,    29,    32,     7,     7,
     115,   116,    87,    88,    87,    88,    89,    90,    91,     0,
      34,    35,    97,    96,    71,    72,    73,     7,   103,   108,
     109,   110,   111,    58,     8,    71,    72,    73,   117,   118,
      63,     8,     8,   122,   123,   124,   125,   126,   127,   128,
     129,   130,   131,   132,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   142,   143,   170,   171,   104,   104,     8,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   163,     3,     4,     5,    87,    88,
     104,    99,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,     7,   103,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
       8,   210,   211,     7,    96,    89,    90,    91,    38,    39,
      40,     8,    96,   105,     8,     8,     4,   226,     5,     4,
       6,   103,     5,     5,   233,   234,   235,   236,     9,     5,
       7,   103,     7,     7,   100,   244,   245,     7,     7,    87,
      88,     7,    99,     8,    92,     7,     7,    97,     8,    97,
       8,     8,     7,    97,    97,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    97,
      97,    97,    97,    96,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,     0,     1,
      97,    97,    96,    -1,    97,    97,    97,    78,    79,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    54,    49,    97,    -1,    96,    -1,    97,    97,    97,
      32,    33,    34,    35,    36,    37,    97,    76,    96,    41,
      42,    -1,    -1,    97,    46,    -1,    -1,    -1,    97,    51,
      52,    97,    -1,    55,    97,    -1,    58,    -1,    60,    61,
      62,    -1,    -1,    -1,    66,    -1,    -1,    -1,    70,    71,
      72,    73,    32,    33,    34,    35,    36,    37,    -1,    -1,
      -1,    41,    42,    -1,    -1,    -1,    46,    -1,    -1,    -1,
      -1,    51,    52,    -1,    -1,    55,    -1,    -1,    58,    -1,
      60,    61,    62,    -1,    -1,    -1,    66,    -1,    -1,    -1,
      70,    71,    72,    73,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    -1,    -1,
      -1,    -1,    96,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   105,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,
      96,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   105,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   105,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    -1,    -1,    -1,    -1,    96,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   105,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      -1,    -1,    -1,    -1,    96,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   105,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    -1,    -1,
      -1,    -1,    96,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   105,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,
      96,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   105,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   104,    78,    79,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    -1,    -1,    -1,    -1,    96,    -1,    -1,    -1,   100,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,
      98,    78,    79,    80,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,
      -1,    98,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,
      96,    -1,    98,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    -1,    -1,    -1,
      -1,    96,    -1,    98,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    -1,    -1,
      -1,    -1,    96,    -1,    98,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    -1,
      -1,    -1,    -1,    96,    -1,    98,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      -1,    -1,    -1,    -1,    96,    -1,    98,    78,    79,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    -1,    -1,    -1,    -1,    96,    -1,    98,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    -1,    -1,    -1,    -1,    96,    -1,    98,    78,
      79,    80,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,    98,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,
      98,    78,    79,    80,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,
      -1,    98,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,
      96,    -1,    98,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    -1,    -1,    -1,
      -1,    96,    -1,    98,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    -1,    -1,
      -1,    -1,    96,    -1,    98,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    -1,
      -1,    -1,    -1,    96,    -1,    98,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      -1,    -1,    -1,    -1,    96,    -1,    98,    78,    79,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    -1,    -1,    -1,    -1,    96,    -1,    98,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    -1,    -1,    -1,    -1,    96,    -1,    98,    78,
      79,    80,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,    98,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    -1,
      98,    78,    79,    80,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,
      -1,    98,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    -1,    -1,    -1,    -1,
      96,    -1,    98,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    -1,    -1,    -1,
      -1,    96,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      -1,    -1,    -1,    -1,    96,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    -1,    -1,    -1,    -1,    96
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,   107,   108,     7,     0,    32,    33,    34,    35,
      36,    37,    41,    42,    46,    51,    52,    55,    58,    60,
      61,    62,    66,    70,    71,    72,    73,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   120,     7,     8,     8,
       8,     8,    99,     8,     8,     8,     8,   112,   115,   103,
       4,     5,     6,     4,   103,     5,     5,     5,   114,     7,
       7,     7,   100,   111,     7,     7,     7,   104,     8,   104,
       7,   103,     7,    38,    39,    40,   118,   119,     8,     8,
       8,   104,   119,     3,     4,     5,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    87,    88,
      92,    97,   121,   122,   121,    87,    88,    97,   103,   123,
     124,     7,    99,    97,    97,    97,    97,    97,    97,    97,
      97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
      97,    97,    97,    97,   121,   121,   121,   121,     7,    78,
      79,    80,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    96,     7,   123,   123,   121,   121,     7,
      87,    88,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,    98,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     105,   105,   123,   123,   100,    98,    98,    98,    98,    98,
      98,    98,    98,    98,    98,    98,   105,    98,    98,    98,
      98,    98,    98,   105,   105,   105,     9,   121,   121,   121,
     121,   121,   121,   121,   105,   105,    98,    98,    98,    98,
     121,   121,    98,   104
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int8 yyr1[] =
{
       0,   106,   107,   107,   108,   108,   109,   109,   109,   109,
     109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
     109,   109,   110,   110,   111,   111,   112,   112,   112,   113,
     113,   114,   114,   115,   115,   116,   117,   118,   118,   119,
     119,   119,   120,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   122,   122,   122,   122,   122,   123,
     123,   123,   123,   123,   124,   124
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     2,     0,     2,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     6,     1,     1,     2,     4,     4,     4,     6,
       1,     1,     2,     4,     4,     4,     9,     1,     2,     4,
       4,     4,     2,     1,     3,     2,     2,     2,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     5,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     6,     4,     4,     4,     4,     4,     4,
       6,     6,     6,     4,     1,     1,     1,     1,     4,     1,
       2,     2,     3,     3,     7,     7
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (dataset, YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, dataset); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, DataSet_t* dataset)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  YYUSE (dataset);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, DataSet_t* dataset)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep, dataset);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule, DataSet_t* dataset)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[+yyssp[yyi + 1 - yynrhs]],
                       &yyvsp[(yyi + 1) - (yynrhs)]
                                              , dataset);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, dataset); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen(S) (YY_CAST (YYPTRDIFF_T, strlen (S)))
#  else
/* Return the length of YYSTR.  */
static YYPTRDIFF_T
yystrlen (const char *yystr)
{
  YYPTRDIFF_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYPTRDIFF_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYPTRDIFF_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            else
              goto append;

          append:
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (yyres)
    return yystpcpy (yyres, yystr) - yyres;
  else
    return yystrlen (yystr);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYPTRDIFF_T *yymsg_alloc, char **yymsg,
                yy_state_t *yyssp, int yytoken)
{
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat: reported tokens (one for the "unexpected",
     one per "expected"). */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Actual size of YYARG. */
  int yycount = 0;
  /* Cumulated lengths of YYARG.  */
  YYPTRDIFF_T yysize = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[+*yyssp];
      YYPTRDIFF_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
      yysize = yysize0;
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYPTRDIFF_T yysize1
                    = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
                    yysize = yysize1;
                  else
                    return 2;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    /* Don't count the "%s"s in the final size, but reserve room for
       the terminator.  */
    YYPTRDIFF_T yysize1 = yysize + (yystrlen (yyformat) - 2 * yycount) + 1;
    if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
      yysize = yysize1;
    else
      return 2;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          ++yyp;
          ++yyformat;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, DataSet_t* dataset)
{
  YYUSE (yyvaluep);
  YYUSE (dataset);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (DataSet_t* dataset)
{
    yy_state_fast_t yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss;
    yy_state_t *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYPTRDIFF_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYPTRDIFF_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
# undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 3:
#line 113 "Parser.y"
               { yyerrok; return 1; }
#line 1793 "Parser.tab.c"
    break;

  case 5:
#line 123 "Parser.y"
    {
      Message_Direct("%s has been read\n",(yyvsp[0].c)) ;
    }
#line 1801 "Parser.tab.c"
    break;

  case 26:
#line 169 "Parser.y"
    {
      InternationalSystemOfUnits_UseAsLength((yyvsp[-1].c)) ;
      printf("%s = { %s }\n",(yyvsp[-3].c),(yyvsp[-1].c));
    }
#line 1810 "Parser.tab.c"
    break;

  case 27:
#line 174 "Parser.y"
    {
      InternationalSystemOfUnits_UseAsTime((yyvsp[-1].c)) ;
      printf("%s = { %s }\n",(yyvsp[-3].c),(yyvsp[-1].c));
    }
#line 1819 "Parser.tab.c"
    break;

  case 28:
#line 179 "Parser.y"
    {
      InternationalSystemOfUnits_UseAsMass((yyvsp[-1].c)) ;
      printf("%s = { %s }\n",(yyvsp[-3].c),(yyvsp[-1].c));
    }
#line 1828 "Parser.tab.c"
    break;

  case 33:
#line 205 "Parser.y"
    {
      Geometry_t* geom = DataSet_GetGeometry(dataset) ;
      
      if(!geom) {
        geom = Geometry_New() ;
        DataSet_GetGeometry(dataset) = geom ;
      }

      Geometry_GetDimension(geom) = (yyvsp[-1].i) ;

      printf("%s = { %d } %s\n",(yyvsp[-3].c),(yyvsp[-1].i),(yyvsp[0].c));
    }
#line 1845 "Parser.tab.c"
    break;

  case 34:
#line 218 "Parser.y"
    {
      Geometry_t* geom = DataSet_GetGeometry(dataset) ;
      
      if(!geom) {
        geom = Geometry_New() ;
        DataSet_GetGeometry(dataset) = geom ;
      }

      {
        if(!strncasecmp((yyvsp[-1].c),"plan",4))  {
          Geometry_SetPlaneSymmetry(geom) ;
        } else if(!strncasecmp((yyvsp[-1].c),"axis",4))  {
          Geometry_SetCylindricalSymmetry(geom) ;
        } else if(!strncasecmp((yyvsp[-1].c),"sphe",4))  {
          Geometry_SetSphericalSymmetry(geom) ;
        } else {
          Message_FatalError("Geometry_Create: geometry not available\n\
          Available geometries are: PLAN, AXIS, SPHE") ;
        }
      }

      printf("%s = { %s } %s\n",(yyvsp[-3].c),(yyvsp[-1].c),(yyvsp[0].c));
    }
#line 1873 "Parser.tab.c"
    break;

  case 35:
#line 249 "Parser.y"
    {
      Mesh_t* mesh = (Mesh_t*) Mry_New(Mesh_t) ;

      DataSet_GetMesh(dataset) = mesh ;
      Mesh_GetGeometry(mesh) = DataSet_GetGeometry(dataset) ;
      Mesh_Scan(mesh,(yyvsp[-1].c)) ;
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
#line 1894 "Parser.tab.c"
    break;

  case 39:
#line 286 "Parser.y"
    {
      printf("%s = { %e }\n",(yyvsp[-3].c),(yyvsp[-1].d));
    }
#line 1902 "Parser.tab.c"
    break;

  case 42:
#line 299 "Parser.y"
    {
      YYABORT;
    }
#line 1910 "Parser.tab.c"
    break;

  case 43:
#line 369 "Parser.y"
                                     { (yyval.d) = (yyvsp[0].d);           }
#line 1916 "Parser.tab.c"
    break;

  case 44:
#line 370 "Parser.y"
                                     { (yyval.d) = (yyvsp[-1].d);           }
#line 1922 "Parser.tab.c"
    break;

  case 45:
#line 371 "Parser.y"
                                     { (yyval.d) = -(yyvsp[0].d);          }
#line 1928 "Parser.tab.c"
    break;

  case 46:
#line 372 "Parser.y"
                                     { (yyval.d) = (yyvsp[0].d);           }
#line 1934 "Parser.tab.c"
    break;

  case 47:
#line 373 "Parser.y"
                                     { (yyval.d) = !(yyvsp[0].d);          }
#line 1940 "Parser.tab.c"
    break;

  case 48:
#line 374 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) - (yyvsp[0].d);      }
#line 1946 "Parser.tab.c"
    break;

  case 49:
#line 375 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) + (yyvsp[0].d);      }
#line 1952 "Parser.tab.c"
    break;

  case 50:
#line 376 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) * (yyvsp[0].d);      }
#line 1958 "Parser.tab.c"
    break;

  case 51:
#line 378 "Parser.y"
    { 
      if(!(yyvsp[0].d)) {
  yymsg(0, "Division by zero in '%g / %g'", (yyvsp[-2].d), (yyvsp[0].d));
      } else {
  (yyval.d) = (yyvsp[-2].d) / (yyvsp[0].d);
      }
    }
#line 1970 "Parser.tab.c"
    break;

  case 52:
#line 385 "Parser.y"
                                     { (yyval.d) = (int)(yyvsp[-2].d) % (int)(yyvsp[0].d);  }
#line 1976 "Parser.tab.c"
    break;

  case 53:
#line 386 "Parser.y"
                                     { (yyval.d) = pow((yyvsp[-2].d), (yyvsp[0].d));  }
#line 1982 "Parser.tab.c"
    break;

  case 54:
#line 387 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) < (yyvsp[0].d);      }
#line 1988 "Parser.tab.c"
    break;

  case 55:
#line 388 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) > (yyvsp[0].d);      }
#line 1994 "Parser.tab.c"
    break;

  case 56:
#line 389 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) <= (yyvsp[0].d);     }
#line 2000 "Parser.tab.c"
    break;

  case 57:
#line 390 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) >= (yyvsp[0].d);     }
#line 2006 "Parser.tab.c"
    break;

  case 58:
#line 391 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) == (yyvsp[0].d);     }
#line 2012 "Parser.tab.c"
    break;

  case 59:
#line 392 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) != (yyvsp[0].d);     }
#line 2018 "Parser.tab.c"
    break;

  case 60:
#line 393 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) && (yyvsp[0].d);     }
#line 2024 "Parser.tab.c"
    break;

  case 61:
#line 394 "Parser.y"
                                     { (yyval.d) = (yyvsp[-2].d) || (yyvsp[0].d);     }
#line 2030 "Parser.tab.c"
    break;

  case 62:
#line 395 "Parser.y"
                                     { (yyval.d) = (yyvsp[-4].d) ? (yyvsp[-2].d) : (yyvsp[0].d); }
#line 2036 "Parser.tab.c"
    break;

  case 63:
#line 396 "Parser.y"
                                     { (yyval.d) = exp((yyvsp[-1].d));      }
#line 2042 "Parser.tab.c"
    break;

  case 64:
#line 397 "Parser.y"
                                     { (yyval.d) = log((yyvsp[-1].d));      }
#line 2048 "Parser.tab.c"
    break;

  case 65:
#line 398 "Parser.y"
                                     { (yyval.d) = log10((yyvsp[-1].d));    }
#line 2054 "Parser.tab.c"
    break;

  case 66:
#line 399 "Parser.y"
                                     { (yyval.d) = sqrt((yyvsp[-1].d));     }
#line 2060 "Parser.tab.c"
    break;

  case 67:
#line 400 "Parser.y"
                                     { (yyval.d) = sin((yyvsp[-1].d));      }
#line 2066 "Parser.tab.c"
    break;

  case 68:
#line 401 "Parser.y"
                                     { (yyval.d) = asin((yyvsp[-1].d));     }
#line 2072 "Parser.tab.c"
    break;

  case 69:
#line 402 "Parser.y"
                                     { (yyval.d) = cos((yyvsp[-1].d));      }
#line 2078 "Parser.tab.c"
    break;

  case 70:
#line 403 "Parser.y"
                                     { (yyval.d) = acos((yyvsp[-1].d));     }
#line 2084 "Parser.tab.c"
    break;

  case 71:
#line 404 "Parser.y"
                                     { (yyval.d) = tan((yyvsp[-1].d));      }
#line 2090 "Parser.tab.c"
    break;

  case 72:
#line 405 "Parser.y"
                                     { (yyval.d) = atan((yyvsp[-1].d));     }
#line 2096 "Parser.tab.c"
    break;

  case 73:
#line 406 "Parser.y"
                                     { (yyval.d) = atan2((yyvsp[-3].d), (yyvsp[-1].d));}
#line 2102 "Parser.tab.c"
    break;

  case 74:
#line 407 "Parser.y"
                                     { (yyval.d) = sinh((yyvsp[-1].d));     }
#line 2108 "Parser.tab.c"
    break;

  case 75:
#line 408 "Parser.y"
                                     { (yyval.d) = cosh((yyvsp[-1].d));     }
#line 2114 "Parser.tab.c"
    break;

  case 76:
#line 409 "Parser.y"
                                     { (yyval.d) = tanh((yyvsp[-1].d));     }
#line 2120 "Parser.tab.c"
    break;

  case 77:
#line 410 "Parser.y"
                                     { (yyval.d) = fabs((yyvsp[-1].d));     }
#line 2126 "Parser.tab.c"
    break;

  case 78:
#line 411 "Parser.y"
                                     { (yyval.d) = floor((yyvsp[-1].d));    }
#line 2132 "Parser.tab.c"
    break;

  case 79:
#line 412 "Parser.y"
                                     { (yyval.d) = ceil((yyvsp[-1].d));     }
#line 2138 "Parser.tab.c"
    break;

  case 80:
#line 413 "Parser.y"
                                     { (yyval.d) = fmod((yyvsp[-3].d), (yyvsp[-1].d)); }
#line 2144 "Parser.tab.c"
    break;

  case 81:
#line 414 "Parser.y"
                                     { (yyval.d) = fmod((yyvsp[-3].d), (yyvsp[-1].d)); }
#line 2150 "Parser.tab.c"
    break;

  case 82:
#line 415 "Parser.y"
                                     { (yyval.d) = sqrt((yyvsp[-3].d) * (yyvsp[-3].d) + (yyvsp[-1].d) * (yyvsp[-1].d)); }
#line 2156 "Parser.tab.c"
    break;

  case 83:
#line 416 "Parser.y"
                                     { (yyval.d) = (yyvsp[-1].d) * (double)rand() / (double)RAND_MAX; }
#line 2162 "Parser.tab.c"
    break;

  case 84:
#line 425 "Parser.y"
              { (yyval.d) = (yyvsp[0].i); }
#line 2168 "Parser.tab.c"
    break;

  case 85:
#line 426 "Parser.y"
              { (yyval.d) = (yyvsp[0].d); }
#line 2174 "Parser.tab.c"
    break;

  case 86:
#line 427 "Parser.y"
              { (yyval.d) = 3.141592653589793; }
#line 2180 "Parser.tab.c"
    break;

  case 87:
#line 432 "Parser.y"
    {
      double* v = map_get(&Parser_yysymbols,(yyvsp[0].c));

      if(!v){
  yymsg(0, "Unknown variable '%s'", (yyvsp[0].c));
  (yyval.d) = 0.;
      } else {
        (yyval.d) = v[0];
      }

      Free((yyvsp[0].c));
    }
#line 2197 "Parser.tab.c"
    break;

  case 88:
#line 445 "Parser.y"
    {
      int index = (int) (yyvsp[-1].d);
      double* v = map_get(&Parser_yysymbols,(yyvsp[-3].c));

      if(!v){
  yymsg(0, "Unknown variable '%s'", (yyvsp[-3].c));
  (yyval.d) = 0.;
      } else {
  (yyval.d) = v[index];
      }

      Free((yyvsp[-3].c));
    }
#line 2215 "Parser.tab.c"
    break;

  case 89:
#line 464 "Parser.y"
    {
      memcpy((yyval.v), (yyvsp[0].v), 3*sizeof(double));
    }
#line 2223 "Parser.tab.c"
    break;

  case 90:
#line 468 "Parser.y"
    {
      for(int i = 0; i < 3; i++) (yyval.v)[i] = -(yyvsp[0].v)[i];
    }
#line 2231 "Parser.tab.c"
    break;

  case 91:
#line 472 "Parser.y"
    { 
      for(int i = 0; i < 3; i++) (yyval.v)[i] = (yyvsp[0].v)[i];
    }
#line 2239 "Parser.tab.c"
    break;

  case 92:
#line 476 "Parser.y"
    { 
      for(int i = 0; i < 3; i++) (yyval.v)[i] = (yyvsp[-2].v)[i] - (yyvsp[0].v)[i];
    }
#line 2247 "Parser.tab.c"
    break;

  case 93:
#line 480 "Parser.y"
    {
      for(int i = 0; i < 3; i++) (yyval.v)[i] = (yyvsp[-2].v)[i] + (yyvsp[0].v)[i];
    }
#line 2255 "Parser.tab.c"
    break;

  case 94:
#line 489 "Parser.y"
    { 
      (yyval.v)[0] = (yyvsp[-5].d);  (yyval.v)[1] = (yyvsp[-3].d);  (yyval.v)[2] = (yyvsp[-1].d);
    }
#line 2263 "Parser.tab.c"
    break;

  case 95:
#line 493 "Parser.y"
    {
      (yyval.v)[0] = (yyvsp[-5].d);  (yyval.v)[1] = (yyvsp[-3].d);  (yyval.v)[2] = (yyvsp[-1].d);
    }
#line 2271 "Parser.tab.c"
    break;


#line 2275 "Parser.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (dataset, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = YY_CAST (char *, YYSTACK_ALLOC (YY_CAST (YYSIZE_T, yymsg_alloc)));
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (dataset, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, dataset);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, dataset);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;


#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (dataset, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif


/*-----------------------------------------------------.
| yyreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, dataset);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[+*yyssp], yyvsp, dataset);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 502 "Parser.y"


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
  Parser_t* parser = (Parser_t*) self ;
  
  //free(parser) ;
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

