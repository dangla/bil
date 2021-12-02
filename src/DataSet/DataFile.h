#ifndef DATAFILE_H
#define DATAFILE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct DataFile_s     ; typedef struct DataFile_s     DataFile_t ;


#include <stdio.h>

extern DataFile_t*  (DataFile_Create)(char*) ;
extern void         (DataFile_Delete)(void*) ;
extern char*        (DataFile_SetFilePositionAfterKey)(DataFile_t*,const char*,const char*,short int) ;
extern char*        (DataFile_ReadLineFromCurrentFilePosition)(DataFile_t*) ;
extern char*        (DataFile_ReadLineFromCurrentFilePositionInString)(DataFile_t*) ;
//extern void*        (DataFile_ReadArray)(DataFile_t*,const char*,void*,int,size_t) ;
extern int*         (DataFile_ReadInversePermutationOfNodes)(DataFile_t*,int) ;



/*
 *  Function-like macros
 */
#include "TextFile.h"

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
        TextFile_Rewind(DataFile_GetTextFile(DF))
        
#define DataFile_ReadArrayFromCurrentFilePosition(DF, ...) \
        TextFile_ReadArrayFromCurrentFilePosition(DataFile_GetTextFile(DF),__VA_ARGS__)
        
#define DataFile_ReadDoublesFromCurrentFilePosition(DF, ...) \
        DataFile_ReadArrayFromCurrentFilePosition(DF,"%le",__VA_ARGS__)

#define DataFile_RemoveComments(DF) \
        TextFile_RemoveComments(DataFile_GetTextFile(DF))



#include "String.h"

/* Tokens in file content */
#define DataFile_FindToken(DF, ...) \
        String_FindToken(DataFile_GetFileContent(DF),__VA_ARGS__)
        
#define DataFile_FindNthToken(DF, ...) \
        String_FindNthToken(DataFile_GetFileContent(DF),__VA_ARGS__)
        
#define DataFile_CountTokens(DF, ...) \
        String_CountTokens(DataFile_GetFileContent(DF),__VA_ARGS__)
        
#define DataFile_CountNbOfKeyWords(DF, ...) \
        String_CountTokens(DataFile_GetFileContent(DF),__VA_ARGS__)
        
        

#define DataFile_MaxLengthOfFileName    (TextFile_MaxLengthOfFileName)
#define DataFile_MaxLengthOfTextLine    (TextFile_MaxLengthOfTextLine)
#define DataFile_MaxLengthOfKeyWord     (30)

#define DataFile_MaxNbOfKeyWords        (10)
#define DataFile_MaxLengthOfKeyWords    (DataFile_MaxNbOfKeyWords*DataFile_MaxLengthOfKeyWord)



#define DataFile_GetTextFile(DF)              ((DF)->textfile)
#define DataFile_GetTextLine(DF)              ((DF)->line)
#define DataFile_GetInitialization(DF)        ((DF)->initialization)
#define DataFile_GetMaxLengthOfTextLine(DF)   ((DF)->linelength)
#define DataFile_GetParent(DF)                ((DF)->parent)



#define DataFile_GetFileName(DF) \
        TextFile_GetFileName(DataFile_GetTextFile(DF))
        
#define DataFile_GetFileContent(DF) \
        TextFile_GetFileContent(DataFile_GetTextFile(DF))

#define DataFile_GetFileStream(DF) \
        TextFile_GetFileStream(DataFile_GetTextFile(DF))

#define DataFile_GetFilePosition(DF) \
        TextFile_GetFilePosition(DataFile_GetTextFile(DF))

#define DataFile_GetCurrentPositionInString(DF) \
        TextFile_GetCurrentPositionInString(DataFile_GetTextFile(DF))

#define DataFile_GetCurrentPositionInFileContent(DF) \
        TextFile_GetCurrentPositionInFileContent(DataFile_GetTextFile(DF))

#define DataFile_SetCurrentPositionInFileContent(DF,C) \
        TextFile_SetCurrentPositionInFileContent(DataFile_GetTextFile(DF),C) 
        

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
        


#include "DataSet.h"

/* The dataset */
#define DataFile_GetDataSet(DF) \
        ((DataSet_t*) DataFile_GetParent(DF))


/* The sequential index */
#define DataFile_GetSequentialIndex(DF) \
        DataSet_GetSequentialIndex(DataFile_GetDataSet(DF))


#define DataFile_GetNbOfSequences(DF) \
        DataSet_GetNbOfSequences(DataFile_GetDataSet(DF))


struct DataFile_s {
  TextFile_t* textfile ;      /* Text file */
  char* line ;                /* memory space for a line */
  int   initialization ;
  int   linelength ;          /* Length of the longest line */
  void* parent ;
} ;

#endif
