#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Message.h"
#include "Context.h"
#include "Modules.h"


#include "ListOfModules.h"

extern  Module_SetModuleProp_t XMODULES(_SetModuleProp) ;


static void* Modules_Initialize(void*) ;



Modules_t* Modules_Create(void)
/** Create the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = (Modules_t*) malloc(sizeof(Modules_t)) ;
  
  assert(modules) ;
  
  Modules_Initialize(modules) ;
  
  return(modules) ;
}



void* Modules_Initialize(void* self)
/** Initialize the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = (Modules_t*) self ;
  
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



void  Modules_Delete(void* self)
{
  Modules_t** modules = (Modules_t**) self ;
  int n_modules = Modules_GetNbOfModules(*modules) ;
  int i ;
  
  for(i = 0 ; i < n_modules ; i++) {
    Module_t* module_i = Modules_GetModule(*modules) + i ;
    free(Module_GetCodeNameOfModule(module_i)) ;
    free(Module_GetShortTitle(module_i)) ;
  }
  
  free(Modules_GetModule(*modules)) ;
  free(*modules) ;
  *modules = NULL ;
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
