#ifndef TEXTFILE_H
#define TEXTFILE_H


/* vacuous declarations and typedef names */

/* class-like structures */
struct TextFile_s       ; typedef struct TextFile_s       TextFile_t ;


#include <stdio.h>

extern TextFile_t*     (TextFile_Create)(const char*) ;
extern void            (TextFile_Delete)(void*) ;
extern FILE*           (TextFile_OpenFile)(TextFile_t*,const char*) ;
extern void            (TextFile_CloseFile)(TextFile_t*) ;
extern void            (TextFile_StoreFilePosition)(TextFile_t*) ;
extern void            (TextFile_MoveToStoredFilePosition)(TextFile_t*) ;
extern char*           (TextFile_ReadLineFromCurrentFilePosition)(TextFile_t*,char*,int) ;
extern char*           (TextFile_ReadLineFromCurrentFilePositionInString)(TextFile_t*,char*,int) ;
extern long int        (TextFile_CountNbOfCharacters)(TextFile_t*) ;
extern long int        (TextFile_CountNbOfEatenCharacters)(TextFile_t*) ;
extern int             (TextFile_CountTheMaxNbOfCharactersPerLine)(TextFile_t*) ;
extern int             (TextFile_Exists)(TextFile_t*) ;
extern void            (TextFile_CleanTheStream)(TextFile_t*) ;
extern FILE*           (TextFile_FileStreamCopy)(TextFile_t*) ;
extern char*           (TextFile_FileCopy)(TextFile_t*) ;
extern char*           (TextFile_StoreFileContent)(TextFile_t*) ;


/*
 *  Function-like macros
 */

#define TextFile_AllocateInBuffer(TF,sz) \
        Buffer_Allocate(TextFile_GetBuffer(TF),(sz))
        
#define TextFile_FreeBuffer(TF) \
        Buffer_Free(TextFile_GetBuffer(TF))

#define TextFile_Rewind(TF) \
        do { \
          rewind(TextFile_GetFileStream(TF)) ; \
          TextFile_GetCurrentPositionInString(TF) = 0 ; \
        } while(0) 

#define TextFile_DoesNotExist(TF) \
        (!TextFile_Exists(TF))


#include "String.h"

#define TextFile_Scan(TF, ...) \
        String_Scan(TextFile_GetCurrentPositionInFileContent(TF),__VA_ARGS__)
        
#define TextFile_SkipLine(TF) \
        String_SkipLine(TextFile_GetCurrentPositionInFileContent(TF))
        

/** Reads N data of size sizeof(V) with the format "FMT" from 
 *  the current position in the string and advance accordingly. */
#define TextFile_ReadArrayFromCurrentFilePosition(TF,N,FMT,V) \
        do { \
          char* TextFile_c = TextFile_GetCurrentPositionInFileContent(TF) ; \
          String_ReadArray(TextFile_c,N,FMT,V) ; \
        } while(0)


/* Remove comments */
#define TextFile_RemoveComments(TF) \
        String_RemoveComments(TextFile_GetFileContent(TF),TextFile_GetFileContent(TF))



#define TextFile_MaxLengthOfTextLine      (500)
#define TextFile_SizeOfBuffer             (500*sizeof(char))
#define TextFile_MaxLengthOfFileName      (200)



#define TextFile_GetFileName(TF)          ((TF)->filename)
#define TextFile_GetFileStream(TF)        ((TF)->stream)
#define TextFile_GetFilePosition(TF)      ((TF)->pos)
#define TextFile_GetFileContent(TF)       ((TF)->filecontent)
#define TextFile_GetPreviousPositionInString(TF)   ((TF)->prestrpos)
#define TextFile_GetCurrentPositionInString(TF)    ((TF)->curstrpos)


/* Not used
#define TextFile_GetTextLine(TF)          ((TF)->line)
#define TextFile_GetBuffer(TF)            ((TF)->buffer)
#define TextFile_GetNbOfCharacters(TF)    ((TF)->ccount)
#define TextFile_GetNbOfWords(TF)         ((TF)->wcount)
#define TextFile_GetNbOfLines(TF)         ((TF)->lcount)
#define TextFile_GetMaxNbOfCharactersPerLine(TF)    ((TF)->linelength)
*/

#define TextFile_GetCurrentPositionInFileContent(TF) \
        (TextFile_GetFileContent(TF) + TextFile_GetCurrentPositionInString(TF))

#define TextFile_SetCurrentPositionInFileContent(TF,C) \
        do { \
          TextFile_GetCurrentPositionInString(TF) = C - TextFile_GetFileContent(TF) ; \
        } while(0)





//#include "Buffer.h"

struct TextFile_s {           /* File */
  char*     filename ;        /* Name of the file */
  char*     filecontent ;
  FILE*     stream ;          /* Current file stream if any */
  fpos_t*   pos ;             /* Previous stored file position of the stream */
  int       prestrpos ;       /* Previous position in the string file content */
  int       curstrpos ;       /* Current position in the string file content */
  /* char*     line ;            *//* memory space for a line */
  /* Buffer_t* buffer ;          *//* Buffer */
  /* long int ccount ;           *//* Nb of characters in file */
  /* long int wcount ;           *//* Nb of words in file */
  /* long int lcount ;           *//* Nb of lines in file */
  /* int linelength ;            *//* Length of the longest line */
} ;

#endif
