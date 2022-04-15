#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Models.h"
#include "Model.h"

#include "Message.h"
#include "Geometry.h"
#include "DataFile.h"
#include "ObVals.h"
#include "Views.h"
#include "Mry.h"
#include "Math_.h"


static Models_t* (Models_CreateAll)(Geometry_t*) ;


Models_t* (Models_New)(const int n_models)
{
  Models_t* models = (Models_t*) Mry_New(Models_t) ;
  
  Models_GetMaxNbOfModels(models) = n_models ;
  Models_GetNbOfModels(models) = 0 ;
  Models_GetModel(models) = NULL ;
    
  if(n_models) {
    Model_t* model = (Model_t*) Mry_New(Model_t[n_models]) ;
      
    {
      int i ;
  
      for(i = 0 ; i < n_models ; i++) {
        Model_t* mod = Model_New() ;
    
        model[i] = mod[0] ;
        free(mod) ;
      }
    }
      
    Models_GetModel(models) = model ;
  }

  return(models) ;
}



void  (Models_Delete)(void* self)
{
  Models_t* models = (Models_t*) self ;
  
  {
    int n_models = Models_GetMaxNbOfModels(models) ;
    Model_t* model = Models_GetModel(models) ;
    
    {
      int i ;
      
      for(i = 0 ; i < n_models ; i++) {
        Model_t* modeli = model + i ;
        
        Model_Delete(modeli) ;
      }
    }
    
    free(model) ;
    Models_GetModel(models) = NULL ;
    Models_GetMaxNbOfModels(models) = 0 ;
    Models_GetNbOfModels(models) = 0 ;
  }
}



Models_t* (Models_Create)(DataFile_t* datafile,Geometry_t* geom)
{
  int n_mats = DataFile_CountTokens(datafile,"MATE,Material",",") ;
  int n_models = DataFile_CountMainTokens(datafile,"Model",",") ;
  int nmax = Math_Max(n_mats,n_models) ;
  Models_t* models  = Models_New(nmax) ;
    
    
  if(n_models <= 0) {
    return(models) ;
  }
  
  
  Message_Direct("Enter in %s","Models") ;
  Message_Direct("\n") ;


  {
    int i ;
    
    /* Read the models */
    for(i = 0 ; i < n_models ; i++) {
      Model_t* model = Models_GetModel(models) + i ;
      char* c = DataFile_FindNthToken(datafile,"Model",",",i + 1) ;
      
      c = String_SkipLine(c) ;
      
      DataFile_SetCurrentPositionInFileContent(datafile,c) ;
      
      Message_Direct("Enter in %s %d","Model",i+1) ;
      Message_Direct("\n") ;
      
      Model_Scan(model,datafile,geom) ;
      
      Models_GetNbOfModels(models) += 1 ;
    }
  }
  
  return(models) ;
}



Models_t* (Models_CreateAll)(Geometry_t* geom)
/** Create the models found in "ListOfModels.h"  */
{
  Models_t* models = Models_New(Models_NbOfModels) ;
  
  {
    int n = Models_NbOfModels ;
    const char* modelnames[] = {Models_ListOfNames} ;
    int   i ;
    
    Models_GetNbOfModels(models) = n ;
  
    for(i = 0 ; i < n ; i++) {
      Model_t* model_i = Models_GetModel(models) + i ;
      
      Model_Initialize(model_i,modelnames[i],geom,NULL) ;
    }
  }
  
  return(models) ;
}



void (Models_Print)(char* codename,FILE* ficd)
{
  Geometry_t geom = {3} ;
  Models_t* models = Models_CreateAll(&geom) ;
  int n_models = Models_GetNbOfModels(models) ;
  Model_t* model = Models_GetModel(models) ;

  if(!codename) { /* all */
    int i ;
    
    if(!ficd) {
      printf("  Model         |  Authors       |  Short Title\n") ;
      /*     "  14c...........|  14c...........|  ...........\n" */
      printf("----------------|----------------|-------------\n") ;
    }

    for(i = 0 ; i < n_models ; i++) {
      Model_t* model_i = model + i ;
      
      printf("  %-14.14s|",Model_GetCodeNameOfModel(model_i)) ;
      printf("  %-14.14s|",Model_GetNameOfAuthors(model_i)) ;
      printf("  %-s",Model_GetShortTitle(model_i)) ;
      printf("\n") ;
      
      if(ficd) {
        printf("Model = %s : ",Model_GetCodeNameOfModel(model_i)) ;
        Model_PrintModelProp(model_i,ficd) ;
        printf("=============================================") ;
        printf("\n\n") ;
      }
    }
    
  } else {
    Model_t* model_i = Models_FindModel(models,codename) ;
    
    if(model_i) {
      printf("Model = %s : ",codename) ;
      Model_PrintModelProp(model_i,ficd) ;
    }
  }
  
  Models_Delete(models) ;
  free(models) ;
}



int (Models_FindModelIndex)(Models_t* models,const char* codename)
{
  Model_t* model = Models_GetModel(models) ;
  int n_models = Models_GetNbOfModels(models) ;
  int j = 0 ;
  
  while(j < n_models && strcmp(Model_GetCodeNameOfModel(model + j),codename)) j++ ;
  
  if(j < n_models) {
    return(j) ;
  }

  return(-1) ;
}




Model_t* (Models_FindModel)(Models_t* models,const char* codename)
{
  int j = Models_FindModelIndex(models,codename) ;
  
  if(j >= 0) {
    Model_t* model = Models_GetModel(models) ;
    
    return(model + j) ;
  }

  return(NULL) ;
}



Model_t* (Models_FindOrAppendModel)(Models_t* models,const char* codename,Geometry_t* geom,DataFile_t* datafile)
{
  Model_t* model = Models_FindModel(models,codename) ;
      
  if(!model) {
    Model_t* model0 = Models_GetModel(models) ;
    int nmax = Models_GetMaxNbOfModels(models) ;
    int n = Models_GetNbOfModels(models) ;
    
    if(n < nmax) {
      model = model0 + n ;
      Model_Initialize(model,codename,geom,datafile) ;
      Models_GetNbOfModels(models) = n + 1 ;
    } else {
      arret("Models_FindOrAppendModel: cannot append") ;
    }
  }

  return(model) ;
}
