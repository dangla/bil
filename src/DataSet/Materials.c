#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <assert.h>
#include "Message.h"
#include "DataFile.h"
#include "Materials.h"
#include "Models.h"
#include "Curves.h"


/* Extern functions */

Materials_t* (Materials_New)(const int n_mats)
{
  Materials_t* materials   = (Materials_t*) malloc(sizeof(Materials_t)) ;
  
  assert(materials) ;


  Materials_GetNbOfMaterials(materials) = n_mats ;

  /* Allocate the materials */
  {
    Materials_GetMaterial(materials) = Material_Create(n_mats) ;
  }
  
  
  /* Allocate the space for the models used by the materials */
  {
    /* We create the space for n_mats models max */
    int n_models = n_mats ;
    Models_t* usedmodels = Models_New(n_models) ;
    
    Materials_GetUsedModels(materials) = usedmodels ;
  }
  
  
  return(materials) ;
}



Materials_t* (Materials_Create)(DataFile_t* datafile,Geometry_t* geom)
{
  int n_mats = DataFile_CountNbOfKeyWords(datafile,"MATE,Material",",") ;
  Materials_t* materials = Materials_New(n_mats) ;
  
  
  DataFile_OpenFile(datafile,"r") ;

  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
    
      DataFile_SetFilePositionAfterKey(datafile,"MATE,Material",",",i + 1) ;
  
      Message_Direct("Enter in %s %d","Material",i+1) ;
      Message_Direct("\n") ;

      /* Which model ? */
      {
        char   codename[Material_MaxLengthOfKeyWord] ;
        char*  code = codename + 1 ;
        char*  line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
      
        if(!strncmp(line,"Model",5)) {
          sscanf(line,"%*s = %s",codename + 1) ;
        } else {
          sscanf(line,"%s",codename + 1) ;
        }
      
        if(isdigit(codename[1])) {
          codename[0] = 'm' ;
          code = codename ;
        }
      
        /* Code name of the model */
        strcpy(Material_GetCodeNameOfModel(mat),code) ;
      }
      
      
      /* Find or append a model and point to it */
      {
        char*  code = Material_GetCodeNameOfModel(mat) ;
        Models_t* usedmodels = Materials_GetUsedModels(materials) ;
        Model_t* matmodel = Models_FindOrAppendModel(usedmodels,code,geom,datafile) ;
        
        Material_GetModel(mat) = matmodel ;
        //Material_GetModel(mat) = Models_FindModel(models,code) ;
      }

      /* for compatibility with old version */
      if(Material_GetModel(mat)) {
        mat->eqn = Material_GetNameOfEquation(mat) ;
        mat->inc = Material_GetNameOfUnknown(mat) ;
      }

      /* Input material data */
      /* A model pointing to a null pointer serves to build curves only */
      Material_GetNbOfProperties(mat) = Material_ReadProperties(mat,datafile) ;
    
      /* for compatibility with old version */
      if(Material_GetModel(mat)) {
        if(Material_GetNbOfEquations(mat) == 0) {
          Material_GetNbOfEquations(mat) = mat->neq ;
        }
      }
      mat->nc = Material_GetNbOfCurves(mat) ;
    
      if(!Material_GetModel(mat)) {
        Message_Info("Materials_Create(1): Model not known") ;
        exit(EXIT_SUCCESS) ;
      }

    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(materials) ;
}
