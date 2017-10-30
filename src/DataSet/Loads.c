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
#include "Loads.h"


Loads_t* Loads_Create(DataFile_t *datafile,Fields_t *fields, Functions_t *functions)
{
  Loads_t *loads       = (Loads_t*) malloc(sizeof(Loads_t)) ;
  int n_loads ;
  int    i ;
  
  if(!loads) arret("Loads_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"CHAR,LOAD,Loads",",",1) ;
  
  Message_Direct("Enter in %s","Loads") ;
  Message_Direct("\n") ;
  
  
  {
    char *line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    
    n_loads = atoi(line) ;
    
    Loads_GetNbOfLoads(loads) = n_loads ;
    
    if(n_loads <= 0) return(loads) ;
  }
  
  
  {
    Load_t *load = (Load_t *) malloc(n_loads*sizeof(Load_t)) ;
    
    if(!load) arret("Loads_Create (1)") ;
    
    Loads_GetLoad(loads) = load ;
  }
  
  
  {
    char *type =  (char *) malloc(n_loads*Load_MaxLengthOfKeyWord*sizeof(char)) ;
    
    if(!type) arret("Loads_Create (2)") ;
  
    for(i = 0 ; i < n_loads ; i++) {
      Load_t *load = Loads_GetLoad(loads) + i ;
      Load_GetType(load) = type + i*Load_MaxLengthOfKeyWord ;
    }
  }
  
  
  {
    char *name_eqn = (char *) malloc(n_loads*Load_MaxLengthOfKeyWord*sizeof(char)) ;
    
    if(!name_eqn) arret("Loads_Create (3)") ;
  
    for(i = 0 ; i < n_loads ; i++) {
      Load_t *load = Loads_GetLoad(loads) + i ;
      Load_GetNameOfEquation(load) = name_eqn + i*Load_MaxLengthOfKeyWord ;
    }
  }


  for(i = 0 ; i < n_loads ; i++) {
    Load_t *load = Loads_GetLoad(loads) + i ;
    char   *line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    char   *pline ;
    
    /* Region */
    if((pline = strstr(line,"Reg"))) {
      pline = strchr(pline,'=') + 1 ;
      Load_GetRegionIndex(load) = atoi(pline) ;
      
    } else {
      arret("Loads_Create (4) : no Region") ;
    }
    
    /* Equation */
    if((pline = strstr(line,"Equ"))) {
      pline = strchr(pline,'=') + 1 ;
      sscanf(pline," %s",Load_GetNameOfEquation(load)) ;
      
      if(strlen(Load_GetNameOfEquation(load)) > Load_MaxLengthOfKeyWord) {
        arret("Loads_Create (5) : mot trop long") ;
      }
      
      if(isdigit(Load_GetNameOfEquation(load)[0])) {
        if(atoi(Load_GetNameOfEquation(load)) < 1) {
          arret("Loads_Create (6) : numero non positif") ;
        }
      }
      
    } else {
      arret("Loads_Create (4) : no Equation") ;
    }
    
    /* Type */
    if((pline = strstr(line,"Type"))) {
      pline = strchr(pline,'=') + 1 ;
      sscanf(pline," %s",Load_GetType(load)) ;
      
      if(strlen(Load_GetType(load)) > Load_MaxLengthOfKeyWord)  {
        arret("Loads_Create (4) : mot trop long") ;
      }
      
    } else {
      arret("Loads_Create (4) : no Type") ;
    }
    
    /* Field */
    if((pline = strstr(line,"Field")) || (pline = strstr(line,"Champ"))) {
      int n_fields = Fields_GetNbOfFields(fields) ;
      int    ich ;
      pline = strchr(pline,'=') + 1 ;
      sscanf(pline," %d",&ich) ;
    
      if(ich > n_fields) {
        arret("Loads_Create (7) : champ non defini") ;
      } else if(ich > 0) {
        Field_t *field = Fields_GetField(fields) ;
        Load_GetField(load) = field + ich - 1 ;
      } else {
        Load_GetField(load) = NULL ;
      }
      
    } else {
        arret("Loads_Create (4) : no Field") ;
    }
    
    /* Function */
    if((pline = strstr(line,"Func")) || (pline = strstr(line,"Fonc"))) {
      int n_fcts = Functions_GetNbOfFunctions(functions) ;
      int    ifn ;
      pline = strchr(pline,'=') + 1 ;
      sscanf(pline,"%d",&ifn) ;
      
      if(ifn > n_fcts) {
        arret("Loads_Create (8) : fonction non definie") ;
      } else if(ifn > 0) {
        Function_t *function = Functions_GetFunction(functions) ;
        Load_GetFunction(load) = function + ifn - 1 ;
      } else {
        Load_GetFunction(load) = NULL ;
      }
    
    } else {
      arret("Loads_Create (4) : no Function") ;
    }
    
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(loads) ;
}
