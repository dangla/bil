#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Message.h"
#include "Mesh.h"
#include "BilVersion.h"
#include "OutputFile.h"
#include "TextFile.h"
#include "Models.h"
#include "Mry.h"


char   OutputFile_TypeOfCurrentFile ;



OutputFile_t*   (OutputFile_Create)(char* filename,int nfiles)
{
  OutputFile_t* outputfile = (OutputFile_t*) Mry_New(OutputFile_t[nfiles]) ;
  
  {
    char* c = filename + strlen(filename) - 1 ;
    int j = (int) (*c - '0') ;
    int i ;
    
    for(i = 0 ; i < nfiles ; i++) {
      sprintf(c,"%d",j+i) ;
      
      {
        TextFile_t* textfile = TextFile_Create(filename) ;
          
        OutputFile_GetTextFile(outputfile + i) = textfile ;
      }
    }
  }
    
  return(outputfile) ;
}



void   (OutputFile_Delete)(void* self,int nfiles)
{
  OutputFile_t** poutputfile = (OutputFile_t**) self ;
  OutputFile_t*   outputfile = *poutputfile ;
  int i ;
  
  for(i = 0 ; i < nfiles ; i++) {
    TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i) ;
    
    TextFile_Delete(&textfile) ;
  }
  
  free(outputfile) ;
  *poutputfile = NULL ;
}
