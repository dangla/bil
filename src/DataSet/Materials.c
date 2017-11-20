#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "Message.h"
#include "DataFile.h"
#include "Materials.h"
#include "Models.h"
#include "Curves.h"


/* Extern functions */

Materials_t* (Materials_Create)(DataFile_t* datafile,Geometry_t* geom)
{
  Materials_t* materials   = (Materials_t*) malloc(sizeof(Materials_t)) ;
  
  if(!materials) arret("Materials_Create(0)") ;

  /* Nb of materials */
  {
    int n_mats = DataFile_CountNbOfKeyWords(datafile,"MATE,Material",",") ;
    
    Materials_GetNbOfMaterials(materials) = n_mats ;

    /* Allocate the materials */
    Materials_GetMaterial(materials) = Material_Create(n_mats) ;
  }
  
  
  /* Allocate the space for the models used by the materials */
  {
    Models_t* usedmodels = (Models_t*) malloc(sizeof(Models_t)) ;
    
    if(!usedmodels) arret("Materials_Create(1)") ;
    
    /* We create the space for n_mats models max */
    {
      int n_mats = Materials_GetNbOfMaterials(materials) ;
      int n_models = n_mats ;
      
      Models_GetModel(usedmodels) = Model_Create(n_models) ;
      Models_GetMaxNbOfModels(usedmodels) = n_models ;
      Models_GetNbOfModels(usedmodels) = 0 ;
    }
    
    Materials_GetUsedModels(materials) = usedmodels ;
  }
  
  
  DataFile_OpenFile(datafile,"r") ;

  /* Initialization */
  {
    int n_mats = Materials_GetNbOfMaterials(materials) ;
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
  
  
  /* Used models */
  #if 0
  {
    Models_t* usedmodels = (Models_t*) malloc(sizeof(Models_t)) ;
    
    if(!usedmodels) arret("Materials_Create (2)") ;
    
    Materials_GetUsedModels(materials) = usedmodels ;
    
    /* We allocate space memory for n_mats models */
    {
      Model_t* usedmodel  = (Model_t*) calloc(n_mats,sizeof(Model_t)) ;
    
      if(!usedmodel)  arret("Materials_Create (3)") ;
      
      Models_GetModel(usedmodels) = usedmodel ;
    
      /* We initialize */   
      Models_GetNbOfModels(usedmodels) = 0 ;
    }
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      Model_t* model_i = Material_GetModel(mat) ;
      char* codename_i = Material_GetCodeNameOfModel(mat) ;
      
      if(!Models_FindModel(usedmodels,codename_i)) {
        int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
        Model_t* usedmodel = Models_GetModel(usedmodels) ;
        
        usedmodel[n_usedmodels] = *model_i ;
        n_usedmodels += 1 ;
        Models_GetNbOfModels(usedmodels) = n_usedmodels ;
      }
    }
  }
  #endif
  
  return(materials) ;
}
