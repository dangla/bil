#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "Fields.h"

/*
static void   lit_grille(FieldGrid_t* ,int,char*) ;
*/
//static void           Field_ReadGrid(FieldGrid_t*,int,char*) ;



Fields_t* Fields_New(const int n_fields)
{
  Fields_t* fields = (Fields_t*) Mry_New(Fields_t) ;
  
  
  Fields_GetNbOfFields(fields) = n_fields ;
    
  {
    if(n_fields > 0) {
      Field_t* field = (Field_t*) Mry_New(Field_t[n_fields]) ;
      int i ;
      
      for(i = 0 ; i < n_fields ; i++) {
        Field_t* fld = Field_New() ;
        
        field[i] = fld[0] ;
      }
      
      Fields_GetField(fields) = field ;
    }
  }
  
  return(fields) ;
}



Fields_t* Fields_Create(DataFile_t* datafile)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"CHMP,FLDS,Fields",",") ;
  int n_fields = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Fields_t* fields  = Fields_New(n_fields) ;
  
  Message_Direct("Enter in %s","Fields") ;
  Message_Direct("\n") ;
    
    
  if(n_fields <= 0) {
    return(fields) ;
  }



  {
    int i ;
    
    c = String_SkipLine(c) ;
    
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    /* Read the fields */
    for(i = 0 ; i < n_fields ; i++) {
      Field_t* field = Fields_GetField(fields) + i ;
      
  
      Message_Direct("Enter in %s %d","Field",i+1) ;
      Message_Direct("\n") ;
      
      Field_Scan(field,datafile) ;
    }
  }
  
  return(fields) ;
}
