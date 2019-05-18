#ifndef STRING_H
#define STRING_H



#include <stdio.h>
#include <string.h>
#include <stdarg.h>

static int   String_bytes ;
//static char* String_str ;
//static char String_Save[100] ;


#define String_ScanAdv(STR, ...) \
        (STR += String_Scan(STR,__VA_ARGS__))


#define String_Scan(STR,FMT, ...) \
        (sscanf(STR,String_Fmt(FMT),__VA_ARGS__,&String_bytes) , String_bytes)
        

#define String_FindChar(STR,C) \
        (strchr(STR,C))
        

#define String_FindToken(STR,TOK) \
        (strstr(STR,TOK))
        

#define String_FindAnyChar(STR,Cs) \
        (strpbrk(STR,Cs))
        

#define String_SkipAnyChar(STR,Cs) \
        (STR + strspn(STR,Cs))
        

#define String_SkipEndOfLine(SRT) \
        (String_SkipAnyChar(STR," \n"))
        

#define String_FindEndOfLine(SRT) \
        (String_FindChar(SRT,'\n'))
        

#define String_NextLine(SRT) \
        String_SkipAnyChar(String_FindEndOfLine(SRT),"\n")


#define String_ScanAffectedKeyphrase(STR,KEY) \
        sscanf(STR,"%*[ ]%[^=]",KEY)


/* Implementation */
        
#define String_Fmt(FMT)  FMT"%n"

//        strcat(strcpy(String_Save,FMT),"%n")


#endif
