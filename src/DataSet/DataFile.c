#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "DataFile.h"
#include "Message.h"




DataFile_t*  (DataFile_Create)(char* filename)
{
  DataFile_t* datafile = (DataFile_t*) malloc(sizeof(DataFile_t)) ;
  
  if(!datafile) assert(datafile) ;
  
  
  /* Memory space for textfile */
  {
    TextFile_t* textfile = TextFile_Create(filename) ;
    
    DataFile_GetTextFile(datafile) = textfile ;
  }
  
  
  /* Memory space for line */
  {
    TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
    int n = TextFile_CountTheMaxNbOfCharactersPerLine(textfile) ;
    size_t sz = n*sizeof(char) ;
    char* line = (char*) malloc(sz) ;
    
    if(!line) {
      assert(line) ;
    }
    
    DataFile_GetTextLine(datafile) = line ;
    DataFile_GetMaxLengthOfTextLine(datafile) = n ;
  }
  
  
  /* The datafile content */
  {
    TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
    
    TextFile_StoreFileContent(textfile) ;
    DataFile_GetFileContent(datafile) = TextFile_GetFileContent(textfile) ;
  }
  
  return(datafile) ;
}


void (DataFile_Delete)(void* self)
{
  DataFile_t** pdatafile = (DataFile_t**) self ;
  DataFile_t*   datafile = *pdatafile ;
  
  TextFile_Delete(&DataFile_GetTextFile(datafile)) ;
  free(DataFile_GetTextLine(datafile)) ;
  free(datafile) ;
}


int (DataFile_CountNbOfKeyWords)(DataFile_t* datafile,const char* cle,const char* del)
/** Open the file for reading and return the number of times any token of the series
 *  of tokens that are delimited by any character of "del" in "cle,
 *  occurs in the filename "datafile". */
{
  int    n = 0 ;
  char*   tok[DataFile_MaxNbOfKeyWords] ;
  short int ntok = 1 ;
  
  /* Break "cle" into a series of tokens */
  {
    char  cle1[DataFile_MaxLengthOfKeyWords] ;
  
    strcpy(cle1,cle) ;
    
    tok[0] = strtok(cle1,del) ;
  
    if(tok[0] == NULL) {
      arret("DataFile_CountNbOfKeyWords(1)") ;
    }
  
    while((tok[ntok] = strtok(NULL,del))) ntok++ ;
  
    if(ntok > DataFile_MaxNbOfKeyWords - 1) {
      arret("DataFile_CountNbOfKeyWords(2)") ;
    }
  }
  
  /* Compute the nb of times we find the tokens in datafile */
  {
    char*   line ;
    
    DataFile_OpenFile(datafile,"r") ;
  
    while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && \
        (line[0] != '#')) {
      short int itok = 0 ;
    
      while(itok < ntok && strncmp(line,tok[itok],strlen(tok[itok]))) itok++ ;
      
      if(itok < ntok) n++ ;
      
    }
  
    DataFile_CloseFile(datafile) ;
  }
  
  return(n) ;
}



void (DataFile_SetFilePositionAfterKey)(DataFile_t* datafile,const char* cle,const char* del,short int n)
/** The file "datafile" is assumed open for reading. Then set the file position of
 *  its stream just after the n^th occurence of any token of the series of tokens 
 *  that are delimited by any character of "del" in "cle". */
{
  short int count = 0 ;
  char*   tok[DataFile_MaxNbOfKeyWords] ;
  short int ntok = 1 ;
  
  /* Break "cle" into a series of tokens */
  {
    char   cle1[DataFile_MaxLengthOfKeyWords] ;
  
    strcpy(cle1,cle) ;
  
    tok[0] = strtok(cle1,del) ;
    
    while((tok[ntok] = strtok(NULL,del))) ntok++ ;
  
    if(ntok > DataFile_MaxNbOfKeyWords) {
      arret("DataFile_SetFilePositionAfterKey") ;
    }
  }


  /* Compute the nb of times we find the tokens in datafile */
  {
    char*  line ;
    
    DataFile_Rewind(datafile) ;
  
    while((count < n) && \
          (line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && \
          (line[0] != '#')) {
      short int itok = 0 ;
    
      while((itok < ntok) && strncmp(line,tok[itok],strlen(tok[itok]))) itok++ ;
      
      if(itok < ntok) count++ ;
    
    }
  
  
    if(count < n) {
      arret("DataFile_SetFilePositionAfterKey(1): %s not found",cle) ;
    }
  
    DataFile_StoreFilePosition(datafile) ;
  }
}



char* (DataFile_ReadLineFromCurrentFilePosition)(DataFile_t* datafile)
/** Reads the first non-commented line from the stream at the current position.
 *  Return a pointer to the string line if succeeded or stop if failed. */
{
  char* line = DataFile_GetTextLine(datafile) ;
  TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
  int n = DataFile_GetMaxLengthOfTextLine(datafile) ;
  char* c ;
  
  do {
    
    c = TextFile_ReadLineFromCurrentFilePosition(textfile,line,n) ;
      
  } while((c) && (c[0] == '#')) ;

  return(c) ;
}


double* (DataFile_ReadDoublesFromCurrentFilePosition)(DataFile_t* datafile,double* v,int n)
/** Reads n doubles from the stream at the current position.
 *  Return the pointer to double. */
{
  FILE* str  = DataFile_GetFileStream(datafile) ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    fscanf(str,"%le",v + i) ;
  }
  
  return(v) ;
}


void* (DataFile_ReadDataFromCurrentFilePosition)(DataFile_t* datafile,void* v,int n,size_t sz,const char* fmt)
/** Reads n data of size "sz" with the format "fmt" from the stream 
 *  at the current position. Return the pointer to data. */
/** NOT YET CHECKED */
{
  FILE* str  = DataFile_GetFileStream(datafile) ;
  char* c = (char*) v ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    fscanf(str,fmt,c + i*sz) ;
  }
  
  return(v) ;
}



/* Intern Functions */




#ifdef NOTDEFINED
static long int     (DataFile_GetPosition)(char*,const char*) ;
static long int     (DataFile_GetNthPosition)(char*,const char*,short int) ;

long int (DataFile_GetPosition)(char* filename,const char* cle)
/** Open the file "filename" for reading and return the file position of
 *  its stream just before the first occurence of any token of a series  
 *  of tokens delimited by a space character in the string "cle" */
{
  char   line[DataFile_MaxLengthOfTextLine] ;
  FILE*   ficd ;
  long int pos ;
  char*   tok[10] ;
  char   cle1[10*DataFile_MaxLengthOfKeyWord] ;
  short int ntok = 1 ;
  short int itok ;
  
  strcpy(cle1,cle) ;
  
  tok[0] = strtok(cle1," ") ;
  while((tok[ntok] = strtok(NULL," "))) ntok++ ;
  
  if(ntok > 9) arret("DataFile_GetPosition") ;
  
  ficd = fopen(filename,"r") ;
  if(!ficd) arret("DataFile_GetPosition : can\'t open the file") ;

  do {
    /* position dans le fichier */
    pos = ftell(ficd) ;
    if(!fgets(line,sizeof(line),ficd)) arret("DataFile_GetPosition (1) : error or end of file") ;
    itok = 0 ;
    while(itok < ntok && strncmp(line,tok[itok],strlen(tok[itok]))) itok++ ;
  /* } while(strncmp(line,cle,strlen(cle))) ; */
  } while(itok == ntok) ;
  
  /* if(!strncmp(line,cle,strlen(cle))) { */
  if(itok < ntok) {
    /* fprintf(stdout,"Enter in %s\n",cle) ; */
    fprintf(stdout,"Enter in %s\n",tok[itok]) ;
  } else {
    fprintf(stdout,"%s not found\n",cle) ;
    arret("DataFile_GetPosition (2) : not found") ;
    pos = 0 ;
  }
  
  fclose(ficd) ;
  
  return(pos) ;
}



long int (DataFile_GetNthPosition)(char* filename,const char* cle,short int n)
/** Open the file "filename" for reading and return the file position of
 *  its stream just before the n^th occurence of any token among those 
 *  which are delimited by a space character in the string "cle" */
{
  char   line[DataFile_MaxLengthOfTextLine] ;
  FILE*   ficd ;
  short int count = 1 ;
  long int pos ;
  char*   tok[10] ;
  char   cle1[10*DataFile_MaxLengthOfKeyWord] ;
  short int ntok = 1 ;
  short int itok ;
  
  strcpy(cle1,cle) ;
  
  tok[0] = strtok(cle1," ") ;
  while((tok[ntok] = strtok(NULL," "))) ntok++ ;
  
  if(ntok > 9) arret("DataFile_GetPNthosition") ;
  
  ficd = fopen(filename,"r") ;
  if(!ficd) arret("DataFile_GetNthPosition (1) : can\'t open the file") ;

  do {
    /* position dans le fichier */
    pos = ftell(ficd) ;
    if(!fgets(line,sizeof(line),ficd)) arret("DataFile_GetNthPosition (2)") ;
    itok = 0 ;
    while(itok < ntok && strncmp(line,tok[itok],strlen(tok[itok]))) itok++ ;
  /* } while(strncmp(line,cle,strlen(cle)) || count++ < n) ; */
  } while(itok == ntok || count++ < n) ;
  
  /* if(!strncmp(line,cle,strlen(cle))) { */
  if(itok < ntok) {
    /* fprintf(stdout,"Enter in %s %d\n",cle,n) ; */
    fprintf(stdout,"Enter in %s %d\n",tok[itok],n) ;
  } else {
    fprintf(stdout,"%s not found\n",cle) ;
    arret("DataFile_GetNthPosition (3) : not found") ;
    pos = 0 ;
  }
  
  fclose(ficd) ;
  
  return(pos) ;
}
#endif
