#ifndef TEXTFILE_H
#define TEXTFILE_H


/* vacuous declarations and typedef names */

/* class-like structures */
struct TextFile_s       ; typedef struct TextFile_s       TextFile_t ;


#include <stdio.h>

extern TextFile_t*     (TextFile_Create)(char*) ;
extern void            (TextFile_Delete)(TextFile_t**) ;
extern FILE*           (TextFile_OpenFile)(TextFile_t*,const char*) ;
extern void            (TextFile_CloseFile)(TextFile_t*) ;
extern void            (TextFile_StoreFilePosition)(TextFile_t*) ;
extern void            (TextFile_MoveToStoredFilePosition)(TextFile_t*) ;
extern char*           (TextFile_ReadLineFromCurrentFilePosition)(TextFile_t*,char*,int) ;
extern long int        (TextFile_CountNbOfCharacters)(TextFile_t*) ;
extern int             (TextFile_CountTheMaxNbOfCharactersPerLine)(TextFile_t*) ;
extern int             (TextFile_Exists)(TextFile_t*) ;
extern void            (TextFile_CleanTheStream)(TextFile_t*) ;
extern FILE*           (TextFile_FileStreamCopy)(TextFile_t*) ;
extern char*           (TextFile_FileCopy)(TextFile_t*) ;
extern char*           (TextFile_StoreFileContent)(TextFile_t*) ;


#define TextFile_MaxLengthOfTextLine      (500)
#define TextFile_SizeOfBuffer             (500*sizeof(char))
#define TextFile_MaxLengthOfFileName      (200)



#define TextFile_GetFileName(TF)          ((TF)->filename)
#define TextFile_GetFileStream(TF)        ((TF)->stream)
#define TextFile_GetFilePosition(TF)      ((TF)->pos)
#define TextFile_GetFileContent(TF)       ((TF)->filecontent)
#define TextFile_GetTextLine(TF)          ((TF)->line)
#define TextFile_GetBuffer(TF)            ((TF)->buffer)
#define TextFile_GetNbOfCharacters(TF)    ((TF)->ccount)
#define TextFile_GetNbOfWords(TF)         ((TF)->wcount)
#define TextFile_GetNbOfLines(TF)         ((TF)->lcount)
#define TextFile_GetMaxNbOfCharactersPerLine(TF)    ((TF)->linelength)


#define TextFile_AllocateInBuffer(TF,sz) \
        Buffer_Allocate(TextFile_GetBuffer(TF),(sz))
        
#define TextFile_FreeBuffer(TF) \
        Buffer_Free(TextFile_GetBuffer(TF))

/* Function-like macros */
#define TextFile_Rewind(TF) \
        rewind(TextFile_GetFileStream(TF))

#define TextFile_DoesNotExist(TF) \
        (!TextFile_Exists(TF))



#include "Buffer.h"

struct TextFile_s {           /* File */
  char*     filename ;        /* Name of the file */
  FILE*     stream ;          /* Current file stream if any */
  fpos_t*   pos ;             /* Latest stored file position of the stream */
  char*     filecontent ;
  /* char*     line ;            *//* memory space for a line */
  /* Buffer_t* buffer ;          *//* Buffer */
  /* long int ccount ;           *//* Nb of characters in file */
  /* long int wcount ;           *//* Nb of words in file */
  /* long int lcount ;           *//* Nb of lines in file */
  /* int linelength ;            *//* Length of the longest line */
} ;

#endif
