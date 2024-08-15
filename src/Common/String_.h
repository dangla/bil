#ifndef STRING_H
#define STRING_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <stdio.h>



/* vacuous declarations and typedef names */

/* class-like structures */
struct String_s       ; typedef struct String_s       String_t ;



static int   String_bytes ;
static char* String_pchar ;


//extern char*       (String_Create)         (const char*) ;
//extern void        (String_Delete)         (void*) ;
extern char*       (String_FindToken3)           (const char*,const char*,const char*) ;
extern char*       (String_FindAndSkipToken3)    (const char*,const char*,const char*) ;
extern char*       (String_FindNthToken3)        (const char*,const char*,const int) ;
extern char*       (String_FindNthToken4)        (const char*,const char*,const char*,const int) ;
extern int         (String_CountTokens2)         (const char*,const char*) ;
extern int         (String_CountTokens3)         (const char*,const char*,const char*) ;
extern int         (String_CountTokensAloneInOneLine2)(const char*,const char*) ;
extern int         (String_CountTokensAloneInOneLine3)(const char*,const char*,const char*) ;
extern char**      (String_BreakIntoTokens)      (const char*,const char*) ;
extern int         (String_NbOfTokens)           (char**) ;
extern char*       (String_CopyLine)             (const char*) ;
extern const char* (String_SkipRemainingComments)(const char*) ;
extern int         (String_NbOfUncommentedLines) (const char*,const char*) ;
extern int         (String_FindPositionIndex)    (const char*,const char* const*,const int) ;
extern char*       (String_RemoveComments)       (char*,char*) ;



#define String_MaxLengthOfKeyWord     (30)

#define String_MaxNbOfKeyWords        (10)
#define String_MaxLengthOfKeyWords    (String_MaxNbOfKeyWords*String_MaxLengthOfKeyWord)
#define String_MaxLengthOfLine        (500)




#include "Arg.h"
#include "Tuple.h"
#include "Algos.h"
#include "Logic.h"
#include "Utils.h"


/* Scan string
 * ----------- */
/** Scan a string with a given format and return the nb of characters read. */
#define String_Scan(STR,...) \
        Logic_IF(Logic_GE(Arg_NARG(__VA_ARGS__),2))\
        (String_ScanN,String_Scan2)(STR,__VA_ARGS__)


#define String_ScanStringUntil(STR,KEY,END) \
        String_Scan(STR,"%*[ ]%[^" END "]",KEY)


#define String_ScanAffectedKeyphrase(STR,KEY) \
        String_Scan(STR,"%*[ ]%[^=]",KEY)


#define String_FindAndScanExp(STR,TOK,DEL, ...) \
        ((String_pchar  = String_FindAndSkipToken(STR,TOK,DEL)) ? \
        ((String_pchar  = String_FindToken(String_pchar," "))   ? \
        ((String_Scan(String_pchar,__VA_ARGS__))) : 0) : 0)


/** Reads N data from the string STR with the format "FMT" and store into V . */
#define String_ScanArray(STR,N,FMT,V) \
        do { \
          char* String_c = STR ; \
          int String_i ; \
          for(String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,(V) + String_i) ; \
          } \
        } while(0)
        
/** Reads N data for each entries from the string STR with the format "FMT". */
#define String_ScanArrays(STR,N,FMT,...) \
        do { \
          char* String_c = STR ; \
          int String_i ; \
          for(String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,String_IncrementAll(__VA_ARGS__)) ; \
          } \
        } while(0)
        
#define String_ScanArrays0(STR,N,FMT,...) \
        do { \
          char* String_c = STR ; \
          int String_i ; \
          for(String_i = 0 ; String_i < (N) ; String_i++) { \
            String_c += String_Scan(String_c,FMT,__VA_ARGS__) ; \
            Algos_SEPWITH(String_IncrementAll0(__VA_ARGS__),;) ; \
          } \
        } while(0)


/** Gets the advanced position in the string. */
#define String_GetAdvancedPosition \
        (String_pchar)


#define String_GetPosition \
        (String_pchar)


/* Implementation */
#define String_Scan2(STR,FMT) \
        (sscanf(STR,String_Fmt(FMT),&String_bytes) , \
         String_pchar = STR + String_bytes , \
         String_bytes)
        
#define String_ScanN(STR,FMT, ...) \
        (sscanf(STR,String_Fmt(FMT),__VA_ARGS__,&String_bytes) , \
         String_pchar = STR + String_bytes , \
         String_bytes)

#define String_Fmt(FMT)  FMT"%n"

#define String_Increment(a) \
        ((a) + String_i)
        
#define String_IncrementAll(...) \
        Tuple_SEQ(Algos_MAP(Tuple_TUPLE(__VA_ARGS__),String_Increment))

#define String_Increment0(a) \
        (a++)
        
#define String_IncrementAll0(...) \
        Algos_MAP(Tuple_TUPLE(__VA_ARGS__),String_Increment0)

//        strcat(strcpy(String_Save,FMT),"%n")




/* Find characters
 * --------------- */
#define String_FindChar(STR,C) \
        ((STR) ? (char*) strchr(STR,C) : NULL)
        

#define String_FindAnyChar(STR,Cs) \
        ((STR) ? strpbrk(STR,Cs) : NULL)
        

#define String_FindEndOfLine(STR) \
        String_FindChar(STR,'\n')
        

#define String_FindEndOfString(STR) \
        String_FindChar(STR,'\0')



/* Skip tokens/characters
 * ---------------------- */
#define String_SpaceChars \
        " \f\n\r\t\v"

#define String_BlankChars \
        " \t\r"


#define String_SkipAnyChars(STR,Cs) \
        ((STR) ? (STR) + strspn(STR,Cs) : NULL)
        

#define String_SkipAnyOtherChars(STR,Cs) \
        ((STR) ? (STR) + strcspn(STR,Cs) : NULL)
        

#define String_SkipBlankChars(STR) \
        String_SkipAnyChars(STR,String_BlankChars)
        

#define String_SkipNonBlankChars(STR) \
        String_SkipAnyOtherChars(STR,String_BlankChars)
        

#define String_SkipSpaceChars(STR) \
        String_SkipAnyChars(STR,String_SpaceChars)
        

#define String_SkipNonSpaceChars(STR) \
        String_SkipAnyOtherChars(STR,String_SpaceChars)
        

#define String_SkipLine(STR) \
        ((String_pchar = String_FindEndOfLine(STR)),String_SkipSpaceChars(String_pchar))


#define String_SkipNextToken(STR) \
        String_SkipNonBlankChars(String_SkipBlankChars(STR))
        



/* Compare with characters
 * ----------------------- */
#define String_Is(...) \
        Utils_CAT_NARG(String_Is,__VA_ARGS__)(__VA_ARGS__)
        
#define String_IsNot(...) \
        !String_Is(__VA_ARGS__)

#define String_CaseIgnoredIs(...) \
        Utils_CAT_NARG(String_CaseIgnoredIs,__VA_ARGS__)(__VA_ARGS__)
        
#define String_BeginsWithAnyChar(STR,Cs) \
        ((STR) ? strspn(STR,Cs) : 0)

#define String_BeginsWithSingleLineComment(STR) \
        (String_Is(STR,"#",1) || String_Is(STR,"//",2))

#define String_BeginsWithMultiLineComment(STR) \
        String_Is(STR,"/*",2)

#define String_SkipMultiLineComment(STR) \
        String_FindAndSkipToken(STR,"*/")


/* Implementation */
#define String_Is2(STR,...) \
        ((STR) ? (!strcmp(STR,__VA_ARGS__)) : 0)
        
#define String_Is3(STR,...) \
        ((STR) ? (!strncmp(STR,__VA_ARGS__)) : 0)
        
#define String_CaseIgnoredIs2(STR,...) \
        ((STR) ? (!strcasecmp(STR,__VA_ARGS__)) : 0)
        
#define String_CaseIgnoredIs3(STR,...) \
        ((STR) ? (!strncasecmp(STR,__VA_ARGS__)) : 0)




/* Find tokens
 * ----------- */
#define String_FindToken(...) \
        Utils_CAT_NARG(String_FindToken,__VA_ARGS__)(__VA_ARGS__)
        
#define String_FindTokenEndingLines(...) \
        Utils_CAT_NARG(String_FindTokenEndingLines,__VA_ARGS__)(__VA_ARGS__)
        
#define String_FindNthToken(...) \
        Utils_CAT_NARG(String_FindNthToken,__VA_ARGS__)(__VA_ARGS__)

        
/* Implementation */
#define String_FindToken2(STR,TOK) \
        ((STR) ? (char*) strstr(STR,TOK) : NULL)

#define String_FindTokenEndingLines2(STR,TOK) \
        (String_pchar = String_FindToken2(STR,TOK), \
        (String_IsTheLastTokenOfTheLine(String_pchar) ? String_pchar : NULL))



/* Find and skip tokens
 * -------------------- */
        
#define String_FindAndSkipToken(...) \
        Utils_CAT_NARG(String_FindAndSkipToken,__VA_ARGS__)(__VA_ARGS__)

        
/* Implementation */
#define String_FindAndSkipToken2(STR,TOK) \
        ((String_pchar = String_FindToken(STR,TOK)) ? String_pchar + strlen(TOK) : NULL)



/* Count tokens
 * ------------ */
#define String_CountTokens(...) \
        Utils_CAT_NARG(String_CountTokens,__VA_ARGS__)(__VA_ARGS__)

#define String_CountTokensAloneInOneLine(...) \
        Utils_CAT_NARG(String_CountTokensAloneInOneLine,__VA_ARGS__)(__VA_ARGS__)



/* Test
 * ---- */
#define String_HasAtMostOneTokenBeforeEndOfLine(STR) \
        (String_SkipBlankChars(String_SkipNextToken(STR))[0] == '\n')




#if 0
#define String_GetStringLength(STR)                 ((STR)->length)
#define String_GetStringContent(STR)                ((STR)->head)
#define String_GetCurrentPositionInString(STR)      ((STR)->current)




struct String_s {
  char*     head ;         /* String content */
  int       length ;        /* String length */
  char*     current ;       /* Current position in the string */
} ;
#endif



#ifdef __CPLUSPLUS
}
#endif
#endif
