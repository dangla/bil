#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>

#include "Message.h"
#include "Buffer.h"
#include "TextFile.h"


TextFile_t*   (TextFile_Create)(char* filename)
{
  TextFile_t* textfile   = (TextFile_t*) malloc(sizeof(TextFile_t)) ;
  
  if(!textfile) arret("TextFile_Create(0)") ;
  
  
  /* Initialization */
  {
    TextFile_GetFileStream(textfile) = NULL ;
  }
  

  /* Memory space for the file name */
  {
    char* name = (char*) malloc(TextFile_MaxLengthOfFileName*sizeof(char)) ;
    
    if(!name) {
      arret("TextFile_Create(1)") ;
    }
    
    if(filename) {
      if(strlen(filename) > TextFile_MaxLengthOfFileName) {
        arret("TextFile_Create(3)") ;
      }
      strcpy(name,filename) ;
    }
    
    TextFile_GetFileName(textfile) = name ;
  }
  
  
  /* Memory space for the file positions */
  {
    fpos_t* pos = (fpos_t*) malloc(sizeof(fpos_t)) ;
    
    if(!pos) arret("TextFile_Create(2)") ;
    
    TextFile_GetFilePosition(textfile) = pos ;
  }
  
    
  /* Space allocation for buffer */
  /*
  TextFile_GetBuffer(textfile) = Buffer_Create(TextFile_SizeOfBuffer) ;
  */
  
  
  /* Memory space for the text line to be read */
  /*
  {
    int n = TextFile_MaxLengthOfTextLine ;
    char* line = (char*) malloc(n*sizeof(char)) ;
    
    if(!line) arret("TextFile_Create(5)") ;
    
    TextFile_GetTextLine(textfile) = line ;
    TextFile_GetMaxNbOfCharactersPerLine(textfile) = n ;
  }
  */
  
  return(textfile) ;
}


void (TextFile_Delete)(TextFile_t** textfile)
{
  free(TextFile_GetFileName(*textfile)) ;
  free(TextFile_GetFilePosition(*textfile)) ;
  
  free(*textfile) ;
}



int (TextFile_Exists)(TextFile_t* textfile)
{
  char* filename = TextFile_GetFileName(textfile) ;
  FILE* str ;
  
  /* After C, we can have multiple streams pointing to 
   * the same file open at the same time. So we open a
   * stream without worrying and close it afterwards. */
  str = fopen(filename,"r") ;
  
  if(!str) {
    return(0) ;
  }
  
  fclose(str) ;
  
  return(1) ;
}


FILE* (TextFile_OpenFile)(TextFile_t* textfile,const char* mode)
{
  char* filename = TextFile_GetFileName(textfile) ;
  FILE* str ;
  
  TextFile_CloseFile(textfile) ;
  
  str = fopen(filename,mode) ;
  
  TextFile_GetFileStream(textfile) = str ;
    
  if(!str) {
    Message_RuntimeError("TextFile_OpenFile: failed to open %s\n",filename) ;
  }
  
  return(str) ;
}


void (TextFile_CloseFile)(TextFile_t* textfile)
{
  FILE* str = TextFile_GetFileStream(textfile) ;
  
  if(!str) return ;
  
  fclose(str) ;
  
  TextFile_GetFileStream(textfile) = NULL ;
}


void (TextFile_CleanTheStream)(TextFile_t* textfile)
{
  FILE* str = TextFile_GetFileStream(textfile) ;
  
  if(!str) return ;
  
  fflush(str) ;
}


void (TextFile_StoreFilePosition)(TextFile_t* textfile)
/** Store the current file position of the stream for further reading. */
{
  FILE* str  = TextFile_GetFileStream(textfile) ;
  fpos_t* pos = TextFile_GetFilePosition(textfile) ;
  
  /* Store the current file position of the stream in pos */
  if(fgetpos(str,pos)) {
    arret("TextFile_StoreFilePosition") ;
  }
}



void (TextFile_MoveToStoredFilePosition)(TextFile_t* textfile)
/** Set the file position of the stream to stored value. */
{
  FILE* str  = TextFile_GetFileStream(textfile) ;
  fpos_t* pos = TextFile_GetFilePosition(textfile) ;
  
  /* Set the file position of the stream to the stored position */
  if(fsetpos(str,pos)) {
    arret("TextFile_MoveToStoredFilePosition") ;
  }
}



char* (TextFile_ReadLineFromCurrentFilePosition)(TextFile_t* textfile,char* line,int n)
/** Reads a line from the stream of textfile at the current position 
 *  and stores it into the string pointed to by line. It stops when 
 *  either (n-1) characters are read, the newline character is read,
 *  or the end-of-file is reached, whichever comes first. 
 *  The file of textfile is assumed open for reading. 
 *  Return a pointer to the string line on success or NULL on failure. */
{
  FILE* str  = TextFile_GetFileStream(textfile) ;
  char* c ;
  
  
  /* Reads a non empty line from the stream after the current position */
  do {
    if(feof(str)) return(NULL) ;
    
    c = fgets(line,n,str) ;
    
    if(!c) {
      return(NULL) ;
    }
    
    /* Eliminate the first blank characters */
    if(*c == ' ') c += strspn(c," ") ;
    
  /* } while((*c == '\n') || (*c == '#')) ; */
  } while((*c == '\n')) ;
  /* } while(0) ; */
  
  return(c) ;
}


long int TextFile_CountNbOfCharacters(TextFile_t* textfile)
{
  FILE* str = TextFile_OpenFile(textfile,"r") ;
  long int count = 0 ;
  char c ;
  
  if(!str) {
    arret("TextFile_CountNbOfCharacters") ;
  }
  
  while((c = fgetc(str)) != EOF) {
    
    count++ ;
    
  }
  
  TextFile_CloseFile(textfile) ;
  
  return(count) ;
}


int TextFile_CountTheMaxNbOfCharactersPerLine(TextFile_t* textfile)
{
  FILE* str = TextFile_OpenFile(textfile,"r") ;
  int linelength = 0 ;
  int ll = 0 ;
  char c ;
  
  if(!str) {
    arret("TextFile_CountTheMaxNbOfCharactersPerLine") ;
  }
  
  while((c = fgetc(str)) != EOF) {
    
    if(c != '\n') {
      ll++ ;
    } else {
      if(ll > linelength) linelength = ll ;
      ll = 0 ;
    }
  }
  
  TextFile_CloseFile(textfile) ;
  
  return(linelength) ;
}




char* (TextFile_FileCopy)(TextFile_t* textfile)
/** Copy the file in a temporary file
 *  Return its name. */
{
  char*   targetfile = tmpnam(NULL) ; /* temporary filename */
  FILE*   target  = fopen(targetfile,"w") ;
  
  if(!target) {
    arret("TextFile_FileCopy(2)") ;
  }
  
  {
    char c ;
    FILE*   source  = TextFile_OpenFile(textfile,"r") ;
  
    if(!source) {
      arret("TextFile_FileCopy(1)") ;
    }
    
    while( (c = fgetc(source)) != EOF ) fputc(c,target) ;
  }
    
  TextFile_CloseFile(textfile) ;
  fclose(target) ;
  
  return(targetfile) ;
}


FILE* (TextFile_FileStreamCopy)(TextFile_t* textfile)
/** Copy the file in a temporary binary file
 *  Return the stream of the temporary file */
{
  FILE*   target = tmpfile() ; /* temporary file */
  FILE*   source = TextFile_OpenFile(textfile,"r") ;
  
  if(!source) {
    arret("TextFile_FileStreamCopy(0)") ;
  }
  
  {
    char c ;
    
    while( (c = fgetc(source)) != EOF ) fputc(c,target) ;
  }
    
  TextFile_CloseFile(textfile) ;

  rewind(target) ;
  
  return(target) ;
}
