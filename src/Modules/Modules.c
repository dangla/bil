#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Message.h"
#include "Context.h"
#include "Mry.h"
#include "Modules.h"



static Modules_t* Modules_New(const int) ;



Modules_t* Modules_New(const int n_modules)
{
  Modules_t* modules = (Modules_t*) Mry_New(Modules_t) ;
  
  {
    Modules_GetNbOfModules(modules) = n_modules ;
    Modules_GetModule(modules) = Module_Create(n_modules) ;
  }
  
  return(modules) ;
}



Modules_t* Modules_Create(void)
/** Create the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = Modules_New(Modules_NbOfModules) ;
  
  {
    int n = Modules_NbOfModules ;
    const char* modulenames[] = {Modules_ListOfNames} ;
    int   i ;
  
    Modules_GetNbOfModules(modules) = n ;
    
    for(i = 0 ; i < n ; i++) {
      Module_t* module_i = Modules_GetModule(modules) + i ;
      
      Module_Initialize(module_i,modulenames[i]) ;
    }
  }
  
  return(modules) ;
}





void  Modules_Delete(void* self)
{
  Modules_t** pmodules = (Modules_t**) self ;
  Modules_t*   modules = *pmodules ;
  int n_modules = Modules_GetNbOfModules(modules) ;
  int i ;
  
  for(i = 0 ; i < n_modules ; i++) {
    Module_t* module_i = Modules_GetModule(modules) + i ;
    free(Module_GetCodeNameOfModule(module_i)) ;
    free(Module_GetShortTitle(module_i)) ;
  }
  
  free(Modules_GetModule(modules)) ;
  free(modules) ;
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


#if 0
Modules_t* Modules_Create(void)
/** Create the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = (Modules_t*) Mry_New(Modules_t) ;
  
  Modules_Initialize(modules) ;
  
  return(modules) ;
}
#endif
#if 0
static void* Modules_Initialize(void*) ;
void* Modules_Initialize(void* self)
/** Initialize the modules found in "ListOfModules.h"  */
{
  Modules_t* modules = (Modules_t*) self ;
  
  Modules_GetNbOfModules(modules) = Modules_NbOfModules ;
  Modules_GetModule(modules) = Module_Create(Modules_NbOfModules) ;
  
  {
    const char* modulenames[] = {Modules_ListOfNames} ;
    Module_SetModuleProp_t* xModule_SetModuleProp[] = {Modules_ListOfSetModuleProp} ;
    int   i ;
    
    for(i = 0 ; i < Modules_NbOfModules ; i++) {
      Module_t* module_i = Modules_GetModule(modules) + i ;
    
      Module_CopyCodeNameOfModule(module_i,modulenames[i]) ;
      Module_GetSetModuleProp(module_i) = xModule_SetModuleProp[i] ;
      Module_SetModuleProp(module_i) ;
    }
  }
  
  return(modules) ;
}
#endif
