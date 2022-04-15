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
#include "Mry.h"


/* Extern functions */

Materials_t* (Materials_New)(const int n_mats,Models_t* models)
{
  Materials_t* materials   = (Materials_t*) Mry_New(Materials_t) ;


  Materials_GetNbOfMaterials(materials) = n_mats ;

  /* Allocate the materials */
  {
    Material_t* material   = (Material_t*) Mry_New(Material_t[n_mats]) ;
    int    i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat   = Material_New() ;
      
      material[i] = mat[0] ;
      free(mat) ;
    }
    
    Materials_GetMaterial(materials) = material ;
  }
  
  
  /* Allocate the space for the models used by the materials */
  {
    /* We create the space for n_mats models max */
    int n_models = n_mats ;
    Models_t* usedmodels = (models) ? models : Models_New(n_models) ;
    
    Materials_GetUsedModels(materials) = usedmodels ;
  }
  
  
  /* All materials share the same pointer to usedmodels */
  {
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int    i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat   = Materials_GetMaterial(materials) + i ;
      
      Material_GetUsedModels(mat) = usedmodels ;
    }
  }
  
  
  return(materials) ;
}




#if 1
Materials_t* (Materials_Create)(DataFile_t* datafile,Geometry_t* geom,Fields_t* fields,Functions_t* functions,Models_t* models)
{
  int n_mats = DataFile_CountTokens(datafile,"MATE,Material",",") ;
  Materials_t* materials = Materials_New(n_mats,models) ;
  
  
  Message_Direct("Enter in %s","Materials") ;
  Message_Direct("\n") ;
  

  
  /* Open for Material_ScanProperties1/2 */
  DataFile_OpenFile(datafile,"r") ;

  /* Scan the datafile */
  {
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      char* c = DataFile_FindNthToken(datafile,"MATE,Material",",",i + 1) ;
      
      c = String_SkipLine(c) ;
      
      /* This for Material_ScanProperties1/2 */
      {
        DataFile_SetFilePositionAfterKey(datafile,"MATE,Material",",",i + 1) ;
        {
          char* c1 = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
        }
      }
      
      DataFile_SetCurrentPositionInFileContent(datafile,c) ;
  
      Message_Direct("Enter in %s %d","Material",i+1) ;
      Message_Direct("\n") ;
      
      Material_GetFields(mat) = fields ;
      Material_GetFunctions(mat) = functions ;
      
      Material_Scan(mat,datafile,geom) ;
    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(materials) ;
}
#endif



void (Materials_Delete)(void* self)
{
  Materials_t* materials = (Materials_t*) self ;
  
  {
    int n_mats = Materials_GetNbOfMaterials(materials) ;
    Material_t* material = Materials_GetMaterial(materials) ;
    
    if(material) {
      int i ;
      
      for(i = 0 ; i < n_mats ; i++) {
        Material_Delete(material+i) ;
      }
      //Mry_Delete(material,n_mats,Material_Delete) ;
      free(material) ;
      Materials_GetMaterial(materials) = NULL ;
    }
  }
  
  {
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    
    if(usedmodels) {
      Models_Delete(usedmodels) ;
      free(usedmodels) ;
      Materials_GetUsedModels(materials) = NULL ;
    }
  }
}



#if 0
Materials_t* (Materials_Create)(DataFile_t* datafile,Geometry_t* geom,Fields_t* fields,Functions_t* functions,Models_t* models)
{
  int n_mats = DataFile_CountNbOfKeyWords(datafile,"MATE,Material",",") ;
  Materials_t* materials = Materials_New(n_mats,models) ;
  
  
  /* Fields and functions */
  {
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      
      Material_GetFields(mat) = fields ;
      Material_GetFunctions(mat) = functions ;
    }
  }
  
  
  DataFile_OpenFile(datafile,"r") ;
  
  Message_Direct("Enter in %s","Materials") ;
  Message_Direct("\n") ;

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
        Message_FataError("Materials_Create(1): Model not known") ;
      }

    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(materials) ;
}
#endif
