#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "DataFile.h"
#include "Message.h"
#include "String_.h"
#include "Mry.h"




DataFile_t*  (DataFile_Create)(char* filename)
{
  DataFile_t* datafile = (DataFile_t*) Mry_New(DataFile_t) ;
  
  /* Memory space for textfile */
  {
    TextFile_t* textfile = TextFile_Create(filename) ;
    
    DataFile_GetTextFile(datafile) = textfile ;
  }
  
  
  /* Memory space for line */
  {
    int n = DataFile_MaxLengthOfTextLine ;
    
    if(filename) {
      TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
      
      n = TextFile_CountTheMaxNbOfCharactersPerLine(textfile) ;
    }
    
    DataFile_GetMaxLengthOfTextLine(datafile) = n ;
    
    {
      char* line = (char*) Mry_New(char[n+1]) ;
    
      DataFile_GetTextLine(datafile) = line ;
    }
  }
  
  
  /* The datafile content */
  {
    TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
    
    TextFile_StoreFileContent(textfile) ;
    //DataFile_GetFileContent(datafile) = TextFile_GetFileContent(textfile) ;
  }
  
  return(datafile) ;
}



void (DataFile_Delete)(void* self)
{
  DataFile_t* datafile = (DataFile_t*) self ;
  
  {
    TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
    
    if(textfile) {
      TextFile_Delete(textfile) ;
      free(textfile) ;
      DataFile_GetTextFile(datafile) = NULL ;
    }
  }
  
  {
    char* line = DataFile_GetTextLine(datafile) ;
    
    if(line) {
      free(line) ;
      DataFile_GetTextLine(datafile) = NULL ;
    }
  }
}



char* (DataFile_SetFilePositionAfterKey)(DataFile_t* datafile,const char* cle,const char* del,short int n)
/** The file "datafile" is assumed open for reading. Then set the file position of
 *  its stream just after the n^th occurence of any token of the series of tokens 
 *  that are delimited by any character of "del" in "cle". */
{
  char**    tok = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;


  /* Compute the nb of times we find the tokens in datafile */
  #if 1
  {
    short int count = 0 ;
    char*  line ;
    
    DataFile_Rewind(datafile) ;
  
    while((count < n) && \
          (line = DataFile_ReadLineFromCurrentFilePosition(datafile))) {
      short int itok = 0 ;
    
      while((itok < ntok) && strncmp(line,tok[itok],strlen(tok[itok]))) itok++ ;
      
      if(itok < ntok) count++ ;
    
    }
  
  
    if(count < n) {
      arret("DataFile_SetFilePositionAfterKey(1): %s not found",cle) ;
    }
  
    DataFile_StoreFilePosition(datafile) ;
  }
  #endif
  
  #if 1
  {
    char*  line = NULL ;
    short int itok = 0 ;
    
    while((itok < ntok) && !(line = DataFile_FindNthToken(datafile,tok[itok],n))) itok++ ;
    
    if(line) {
      char* c = String_FindAndSkipToken(line,tok[itok]) ;
        
      DataFile_SetCurrentPositionInFileContent(datafile,c) ;
  
      //DataFile_StoreFilePosition(datafile) ;
        
      return(c) ;
    } else {
      arret("DataFile_SetFilePositionAfterKey(1): %s not found",cle) ;
    }
  }
  #endif
  
  return(NULL) ;
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



char* (DataFile_ReadLineFromCurrentFilePositionInString)(DataFile_t* datafile)
/** Reads the first non-commented line from the datafile at the current position of its string.
 *  Return a pointer to the string line if succeeded or stop if failed. */
{
  char* line = DataFile_GetTextLine(datafile) ;
  TextFile_t* textfile = DataFile_GetTextFile(datafile) ;
  int n = DataFile_GetMaxLengthOfTextLine(datafile) ;
  char* c ;
  
  do {
    
    c = TextFile_ReadLineFromCurrentFilePositionInString(textfile,line,n) ;
      
  } while((c) && (c[0] == '#')) ;

  return(c) ;
}






int*  DataFile_ReadInversePermutationOfNodes(DataFile_t* datafile,int n_no)
/** Read the inverse permutation vector of nodes in a file if it exists 
 *  or initialize it with identity function. Return a pointer to n_no int. 
 **/
{
  int* perm = (int*) Mry_New(int[n_no]) ;

  {
    char   nom_iperm[DataFile_MaxLengthOfFileName] ;
    
    {
      char*  filename = DataFile_GetFileName(datafile) ;
    
      if(strlen(filename) + 12 > DataFile_MaxLengthOfFileName) {
        arret("DataFile_ReadInversePermutationOfNodes") ;
      }
    
      sprintf(nom_iperm,"%s.graph.iperm",filename) ;
    }
  
    {
      FILE* fic_iperm = fopen(nom_iperm,"r") ;
  
      if(!fic_iperm) {
        int  i ;
    
        for(i = 0 ; i < n_no ; i++) perm[i] = i ;
    
      } else {
        int  i ;
    
        for(i = 0 ; i < n_no ; i++) {
          int   j ;
      
          fscanf(fic_iperm,"%d",&j) ;
          perm[j] = i ;
        }
      }
    
      if(fic_iperm) fclose(fic_iperm) ;
    }
  }
  
  return(perm) ;
}



/* Not used from here */

#if 0
int (DataFile_CountNbOfKeyWords)(DataFile_t* datafile,const char* cle,const char* del)
/** Return the number of times any token of the series of tokens
 *  that are delimited by any character of "del" in "cle,
 *  occurs in the filename "datafile". */
{
  char**  tok  = String_BreakIntoTokens(cle,del) ;
  int     ntok = String_NbOfTokens(tok) ;
  int     n = 0 ;
  
  /* Compute the nb of times we find the tokens in datafile */
  {
    char* c  = DataFile_GetFileContent(datafile) ;
    int itok ;
    
    for(itok = 0 ; itok < ntok ; itok++) {
      n += String_CountTokens(c,tok[itok]) ;
    }
  }
  
  return(n) ;
}
#endif



#if 0
void* (DataFile_ReadArray)(DataFile_t* datafile,const char* fmt,void* v,int n,size_t sz)
/** Reads n data of size "sz" with the format "fmt" from the stream 
 *  at the current position. Return the pointer to data. */
/** NOT YET CHECKED */
{
  //FILE* str  = DataFile_GetFileStream(datafile) ;
  char* c = (char*) v ;
  char* cur = DataFile_GetCurrentPositionInFileContent(datafile) ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    cur += String_Scan(cur,fmt,c + i*sz) ;
    //fscanf(str,fmt,c + i*sz) ;
  }
  
  return(v) ;
}
#endif



/* Intern Functions */




#if 0
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
