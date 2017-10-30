#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "Context.h"
#include "Modules.h"


#include "ListOfModules.h"

extern  Module_SetModuleProp_t XMODULES(_SetModuleProp) ;



Modules_t* Modules_Create(void)
/** Create the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = (Modules_t*) malloc(sizeof(Modules_t)) ;
  
  if(!modules) arret("Modules_Create") ;
  
  Modules_GetNbOfModules(modules) = NB_MODULES ;
  Modules_GetModule(modules) = Module_Create(NB_MODULES) ;
  
  
  /* MODULENAMES and XMODULES are defined in "ListOfModules.h" */
  {
    const char* modulenames[NB_MODULES] = {MODULENAMES} ;
    Module_SetModuleProp_t* xModule_SetModuleProp[NB_MODULES] = {XMODULES(_SetModuleProp)} ;
    int   i ;
    
    for(i = 0 ; i < NB_MODULES ; i++) {
      Module_t* module_i = Modules_GetModule(modules) + i ;
    
      Module_CopyCodeNameOfModule(module_i,modulenames[i]) ;
      Module_GetSetModuleProp(module_i) = xModule_SetModuleProp[i] ;
      Module_SetModuleProp(module_i) ;
    }
  }
  
  return(modules) ;
}



void  Modules_Delete(Modules_t** modules)
{
  int n_modules = Modules_GetNbOfModules(*modules) ;
  int i ;
  
  for(i = 0 ; i < n_modules ; i++) {
    Module_t* module_i = Modules_GetModule(*modules) + i ;
    free(Module_GetCodeNameOfModule(module_i)) ;
    free(Module_GetShortTitle(module_i)) ;
  }
  
  free(Modules_GetModule(*modules)) ;
  free(*modules) ;
}



Module_t* Module_Create(int n_module)
/* Create modules */
{
  Module_t* module = (Module_t*) calloc(n_module,sizeof(Module_t)) ;
  int i ;
  
  if(!module) arret("Module_Create") ;
  
  
  for(i = 0 ; i < n_module ; i++) {
    Module_t* module_i = module + i ;
    
    
    /* Allocation of space for the code name */
    {
      size_t sz = Module_MaxLengthOfKeyWord*sizeof(char) ;
      char* code = (char*) malloc(sz) ;
      
      if(!code) arret("Module_Create (1)") ;
      
      Module_GetCodeNameOfModule(module_i) = code ;
    }
    
    
    /* Allocation of space for the short title */
    {
      size_t sz = Module_MaxLengthOfShortTitle*sizeof(char) ;
      char* title = (char*) malloc(sz) ;
      
      if(!title) arret("Module_Create (2)") ;
      
      Module_GetShortTitle(module_i) = title ;
    
      Module_CopyShortTitle(module_i,"\0") ;
    }
    
    
    /* Allocation of space for the author names */
    {
      size_t sz = Module_MaxLengthOfAuthorNames*sizeof(char) ;
      char* names = (char*) malloc(sz) ;
      
      if(!names) arret("Module_Create (3)") ;
      
      Module_GetNameOfAuthors(module_i) = names ;
      
      Module_CopyNameOfAuthors(module_i,"\0") ;
    }
  }
  
  return(module) ;
}



void Modules_Print(char* codename)
{
  Modules_t* modules = Modules_Create() ;
  int n_modules = Modules_GetNbOfModules(modules) ;
  Module_t* module = Modules_GetModule(modules) ;

  if(!codename) { /* all */
    int i ;
    
    printf("  Module     | Short Title\n") ;
    printf("-------------|------------\n") ;

    for(i = 0 ; i < n_modules ; i++) {
      Module_t* module_i = module + i ;
      
      printf("  %-10s | ",Module_GetCodeNameOfModule(module_i)) ;
      printf("%-s",Module_GetShortTitle(module_i)) ;
      printf("\n") ;
    }
    
  } else {
    Module_t* module_i = Modules_FindModule(modules,codename) ;
    
    printf("Module = %s : ",codename) ;
    printf("%-s",Module_GetShortTitle(module_i)) ;
    printf("\n") ;
  }
  
  Modules_Delete(&modules) ;
}




Module_t* Modules_FindModule(Modules_t* modules,const char* codename)
{
  Module_t*  module = Modules_GetModule(modules) ;
  int n_modules = Modules_GetNbOfModules(modules) ;
  int j = 0 ;
  
  while(j < n_modules && strcmp(Module_GetCodeNameOfModule(module + j),codename)) j++ ;
  
  if(j < n_modules) {
    Module_t* module_j = module + j ;
    return(module_j) ;
  }
  
  arret("Modules_FindModule: module %s not known",codename) ;

  return(NULL) ;
}
