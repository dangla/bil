#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Model.h"
#include "Models.h"

#include "Message.h"
#include "Geometry.h"
#include "DataFile.h"
#include "ObVals.h"
#include "Views.h"
#include "Mry.h"


extern  Model_SetModelProperties_t  Models_ListOfSetModelProp ;




Model_t* (Model_New)(void)
{
  Model_t* model = (Model_t*) Mry_New(Model_t) ;
  
  
  {
    /* Allocation of space for name of equations */
    {
      char** names = (char**) Mry_New(char*,Model_MaxNbOfEquations) ;
      
      Model_GetNameOfEquation(model) = names ;
    }
      
    {
      char* name = (char*) Mry_New(char,Model_MaxLengthOfKeyWord*Model_MaxNbOfEquations) ;
    
      {
        int j ;
      
        for(j = 0 ; j < Model_MaxNbOfEquations ; j++) {
          Model_GetNameOfEquation(model)[j] = name + j*Model_MaxLengthOfKeyWord ;
        }
      }
    }
    
    
    /* Allocation of space for name of unknowns */
    {
      char** names = (char**) Mry_New(char*,Model_MaxNbOfEquations) ;
    
      Model_GetNameOfUnknown(model) = names ;
    }
      
    {
      char* name = (char*) Mry_New(char,Model_MaxLengthOfKeyWord*Model_MaxNbOfEquations) ;
    
      {
        int j ;
      
        for(j = 0 ; j < Model_MaxNbOfEquations ; j++) {
          Model_GetNameOfUnknown(model)[j]  = name + j*Model_MaxLengthOfKeyWord ;
        }
      }
    }
    
    
    /* Allocation of space for the sequential indexes of unknowns/equations */
    {
      int* ind = (int*) Mry_New(int,Model_MaxNbOfEquations) ;
    
      Model_GetSequentialIndexOfUnknown(model) = ind ;
    }
    
    
    /* Allocation of space for code name of the model */
    {
      char* name = (char*) Mry_New(char,Model_MaxLengthOfKeyWord) ;
      
      Model_GetCodeNameOfModel(model) = name ;
    }
  
    
    /* Allocation of space for short title of the model */
    {
      char* name = (char*) Mry_New(char,Model_MaxLengthOfShortTitle) ;
      
      Model_GetShortTitle(model) = name ;
    
      Model_CopyShortTitle(model,"\0") ;
    }
  
    
    /* Allocation of space for name of authors */
    {
      char* name = (char*) Mry_New(char,Model_MaxLengthOfAuthorNames) ;
    
      Model_GetNameOfAuthors(model) = name ;
      
      Model_CopyNameOfAuthors(model,"\0") ;
    }
    
    
    /* Allocation of space for objective values */
    {
      ObVal_t* obval = (ObVal_t*) Mry_New(ObVal_t,Model_MaxNbOfEquations) ;
    
      Model_GetObjectiveValue(model) = obval ;
    }
    
    
    /* Allocation of space for views */
    {
      Views_t* views = Views_Create(Model_MaxNbOfViews) ;
      
      Model_GetViews(model) = views ;
    }
  }
  
  return(model) ;
}



void  (Model_Delete)(void* self)
{
  Model_t* model = (Model_t*) self ;
  
  {
    {
      char* name = Model_GetNameOfEquation(model)[0] ;
    
      if(name) {
        free(name) ;
      }
    }
    
    {
      char** names = Model_GetNameOfEquation(model) ;
      
      if(names) {
        free(names) ;
      }
    }
    
    {
      char* name = Model_GetNameOfUnknown(model)[0] ;
    
      if(name) {
        free(name) ;
      }
    }
    
    {
      char** names = Model_GetNameOfUnknown(model) ;
      
      if(names) {
        free(names) ;
      }
    }
    
    {
      int* ind = Model_GetSequentialIndexOfUnknown(model) ;
    
      if(ind) {
        free(ind) ;
      }
    }
    
    {
      char* name = Model_GetCodeNameOfModel(model) ;
      
      if(name) {
        free(name) ;
      }
    }
    
    {
      char* name = Model_GetShortTitle(model) ;
      
      if(name) {
        free(name) ;
      }
    }
    
    {
      char* name = Model_GetNameOfAuthors(model) ;
      
      if(name) {
        free(name) ;
      }
    }
    
    {
      ObVal_t* obval = Model_GetObjectiveValue(model) ;
      
      if(obval) {
        free(obval) ;
      }
    }
    
    {
      Views_t* views = Model_GetViews(model) ;
      
      if(views) {
        Views_Delete(views) ;
        free(views) ;
      }
    }
  }
}



Model_t* (Model_Initialize)(Model_t* model,const char* codename,Geometry_t* geom,DataFile_t* datafile)
{
  int n_models = Models_NbOfModels ;
  const char* modelnames[] = {Models_ListOfNames} ;
  Model_SetModelProperties_t* xModel_SetModelProperties[] = {Models_ListOfSetModelProp} ;
  int i = 0 ;
  
  Model_GetGeometry(model) = geom ; /* Important! */
  Model_GetDataFile(model) = datafile ;
  
  while(i < n_models && strcmp(modelnames[i],codename)) i++ ;
    
  if(i < n_models) {
    Model_CopyCodeNameOfModel(model,modelnames[i]) ;
    Model_GetSetModelProperties(model) = xModel_SetModelProperties[i] ;
    Model_SetModelProperties(model) ; /* Call to SetModelProperties */
    
    return(model) ;
  } else {
    
    Message_Warning("No model named %s",codename) ;
  }
  
  //return(NULL) ;
  return(model) ;
}



void (Model_Scan)(Model_t* model,DataFile_t* datafile,Geometry_t* geom)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  
  /* Code name of the model */
  {
    char   codename[Model_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Name",","," = %s",codename) ;
    
    if(n) {
      Model_Initialize(model,codename,geom,datafile) ;
    } else {
      arret("Model_Scan") ;
    }
  }
      
  /* Name of equations */
  {
    int n = String_FindAndScanExp(line,"Equations",","," = ") ;
        
    if(n) {
      int neq = Model_GetNbOfEquations(model) ;
      char* pline = String_GetAdvancedPosition ;
      int i ;
          
      for(i = 0 ; i < neq ; i++) {
        char  name[Model_MaxLengthOfKeyWord] ;
        
        pline += String_Scan(pline,"%s",name) ;
        
        Model_CopyNameOfEquation(model,i,name) ;
      }
    }
  }
      
  /* Name of unknowns */
  {
    int n = String_FindAndScanExp(line,"Unknowns",","," = ") ;
        
    if(n) {
      int neq = Model_GetNbOfEquations(model) ;
      char* pline = String_GetAdvancedPosition ;
      int i ;
          
      for(i = 0 ; i < neq ; i++) {
        char  name[Model_MaxLengthOfKeyWord] ;
        
        pline += String_Scan(pline,"%s",name) ;
        
        Model_CopyNameOfUnknown(model,i,name) ;
      }
    }
  }
}
