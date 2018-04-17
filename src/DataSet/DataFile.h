#ifndef DATAFILE_H
#define DATAFILE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct DataFile_s     ; typedef struct DataFile_s     DataFile_t ;


#include <stdio.h>

extern DataFile_t*  (DataFile_Create)(char*) ;
extern void         (DataFile_Delete)(DataFile_t**) ;
extern int          (DataFile_CountNbOfKeyWords)(DataFile_t*,const char*,const char*) ;
extern void         (DataFile_SetFilePositionAfterKey)(DataFile_t*,const char*,const char*,short int) ;
extern char*        (DataFile_ReadLineFromCurrentFilePosition)(DataFile_t*) ;
extern double*      (DataFile_ReadDoublesFromCurrentFilePosition)(DataFile_t*,double*,int) ;
extern void*        (DataFile_ReadDataFromCurrentFilePosition)(DataFile_t*,void*,int,size_t,const char*) ;



#include "TextFile.h"

#define DataFile_MaxLengthOfFileName    (TextFile_MaxLengthOfFileName)
#define DataFile_MaxLengthOfTextLine    (TextFile_MaxLengthOfTextLine)
#define DataFile_MaxLengthOfKeyWord     (30)

#define DataFile_MaxNbOfKeyWords        (10)
#define DataFile_MaxLengthOfKeyWords    (DataFile_MaxNbOfKeyWords*DataFile_MaxLengthOfKeyWord)


#define DataFile_GetTextFile(DF)              ((DF)->textfile)
#define DataFile_GetTextLine(DF)              ((DF)->line)
#define DataFile_GetInitialization(DF)        ((DF)->initialization)
#define DataFile_GetFileContent(DF)          ((DF)->filecontent)



#define DataFile_GetFileName(DF) \
        TextFile_GetFileName(DataFile_GetTextFile(DF))

#define DataFile_GetFileStream(DF) \
        TextFile_GetFileStream(DataFile_GetTextFile(DF))

#define DataFile_GetFilePosition(DF) \
        TextFile_GetFilePosition(DataFile_GetTextFile(DF))


/*
 *  Function-like macros
 */
#define DataFile_OpenFile(DF,MODE) \
        TextFile_OpenFile(DataFile_GetTextFile(DF),MODE)

#define DataFile_CloseFile(DF) \
        TextFile_CloseFile(DataFile_GetTextFile(DF))

#define DataFile_Exists(DF) \
        TextFile_Exists(DataFile_GetTextFile(DF))

#define DataFile_DoesNotExist(DF) \
        (!DataFile_Exists(DF))

#define DataFile_StoreFilePosition(DF) \
        TextFile_StoreFilePosition(DataFile_GetTextFile(DF))

#define DataFile_MoveToStoredFilePosition(DF) \
        TextFile_MoveToStoredFilePosition(DataFile_GetTextFile(DF))

#define DataFile_Rewind(DF) \
        rewind(DataFile_GetFileStream(DF))
        

/* Test initialization */
#define DataFile_ContextIsDefaultInitialization(DF) \
        (DataFile_GetInitialization(DF) == 0)

#define DataFile_ContextIsFullInitialization \
        DataFile_ContextIsDefaultInitialization
        
#define DataFile_ContextIsPartialInitialization(DF) \
        (DataFile_GetInitialization(DF) == 1)

#define DataFile_ContextIsNoInitialization(DF) \
        (DataFile_GetInitialization(DF) == 2)

#define DataFile_ContextIsInitialization(DF) \
        (DataFile_GetInitialization(DF) < 2)
        
        
/* Set initialization */
#define DataFile_ContextSetToFullInitialization(DF) \
        do {DataFile_GetInitialization(DF) = 0 ;} while(0)
        
#define DataFile_ContextSetToPartialInitialization(DF) \
        do {DataFile_GetInitialization(DF) = 1 ;} while(0)
        
#define DataFile_ContextSetToNoInitialization(DF) \
        do {DataFile_GetInitialization(DF) = 2 ;} while(0)


struct DataFile_s {
  TextFile_t* textfile ;      /* Text file */
  char* line ;                /* memory space for a line */
  int   initialization ;
  char* filecontent ;

} ;

#endif
