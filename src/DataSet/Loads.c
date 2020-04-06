#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Message.h"
#include "Geometry.h"
#include "DataFile.h"
#include "IntFcts.h"
#include "Functions.h"
#include "Fields.h"
#include "String.h"
#include "Mry.h"
#include "Loads.h"





Loads_t* Loads_New(const int n_loads)
{
  Loads_t* loads = (Loads_t*) Mry_New(Loads_t) ;
  
  
  Loads_GetNbOfLoads(loads) = n_loads ;


  if(n_loads > 0) {
    Load_t* load  = (Load_t*) Mry_New(Load_t[n_loads]) ;
    int i ;

    for(i = 0 ; i < n_loads ; i++) {
      Load_t* ld  = Load_New() ;
      
      load[i] = ld[0] ;
    }

    Loads_GetLoad(loads) = load ;
  }
  

  return(loads) ;
}




Loads_t* Loads_Create(DataFile_t* datafile,Fields_t* fields,Functions_t* functions)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"CHAR,LOAD,Loads",",") ;
  int n_loads = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Loads_t* loads = Loads_New(n_loads) ;
  
  
  Message_Direct("Enter in %s","Loads") ;
  Message_Direct("\n") ;
  
  
  if(n_loads <= 0) {
    return(loads) ;
  }
  


  /* Scan the datafile */
  {
    int i ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i = 0 ; i < n_loads ; i++) {
      Load_t *load = Loads_GetLoad(loads) + i ;
    
      Message_Direct("Enter in %s %d","Load",i+1) ;
      Message_Direct("\n") ;
      
      
      Load_GetFields(load) = fields ;
      Load_GetFunctions(load) = functions ;
    
      Load_Scan(load,datafile) ;
    }
  }
  
  return(loads) ;
}
