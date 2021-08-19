#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Message.h"
#include "Context.h"
#include "Mry.h"
#include "Module.h"
#include "Modules.h"



extern  Module_SetModuleProp_t  Modules_ListOfSetModuleProp ;




static Module_t* Module_New(void) ;




Module_t* Module_New(void)
{
  Module_t* module = (Module_t*) Mry_New(Module_t) ;
  
  {
    /* Allocation of space for the code name */
    {
      char* code = (char*) Mry_New(char[Module_MaxLengthOfKeyWord]) ;
      
      Module_GetCodeNameOfModule(module) = code ;
    }
    
    
    /* Allocation of space for the short title */
    {
      char* title = (char*) Mry_New(char[Module_MaxLengthOfShortTitle]) ;
      
      Module_GetShortTitle(module) = title ;
    
      Module_CopyShortTitle(module,"\0") ;
    }
    
    
    /* Allocation of space for the author names */
    {
      char* names = (char*) Mry_New(char[Module_MaxLengthOfAuthorNames]) ;
      
      Module_GetNameOfAuthors(module) = names ;
      
      Module_CopyNameOfAuthors(module,"\0") ;
    }
    
    /* Initialize the nb of sequences and the sequential index */
    {
      Module_GetNbOfSequences(module) = 1 ;
      Module_GetSequentialIndex(module) = 0 ;
    }
  }
  
  return(module) ;
}



Module_t* Module_Create(int n_module)
/** Create modules */
{
  Module_t* module = (Module_t*) Mry_New(Module_t[n_module]) ;
  
  {
    int i ;
    
    for(i = 0 ; i < n_module ; i++) {
      Module_t* mod = Module_New() ;
    
      module[i] = mod[0] ;
    }
  }
  
  return(module) ;
}



Module_t* (Module_Initialize)(Module_t* module,const char* codename)
{
  int NbOfModules = Modules_NbOfModules ;
  const char* modulenames[] = {Modules_ListOfNames} ;
  Module_SetModuleProp_t* xModule_SetModuleProp[] = {Modules_ListOfSetModuleProp} ;
  int   i = 0 ;
  
  while(i < NbOfModules && strcmp(modulenames[i],codename)) i++ ;
    
  if(i < NbOfModules) {
    Module_CopyCodeNameOfModule(module,modulenames[i]) ;
    Module_GetSetModuleProp(module) = xModule_SetModuleProp[i] ;
    Module_SetModuleProp(module) ; /* Call to SetModuleProp */
    
    return(module) ;
  }
  
  return(module) ;
}




#if 0
//void* Module_Initialize(void*) ;
void* Module_Initialize(void* self)
{
  Module_t* module_i = (Module_t*) self ;
  
  {
    /* Allocation of space for the code name */
    {
      char* code = (char*) Mry_New(char[Module_MaxLengthOfKeyWord]) ;
      
      Module_GetCodeNameOfModule(module_i) = code ;
    }
    
    
    /* Allocation of space for the short title */
    {
      char* title = (char*) Mry_New(char[Module_MaxLengthOfShortTitle]) ;
      
      Module_GetShortTitle(module_i) = title ;
    
      Module_CopyShortTitle(module_i,"\0") ;
    }
    
    
    /* Allocation of space for the author names */
    {
      char* names = (char*) Mry_New(char[Module_MaxLengthOfAuthorNames]) ;
      
      Module_GetNameOfAuthors(module_i) = names ;
      
      Module_CopyNameOfAuthors(module_i,"\0") ;
    }
  }
  
  return(module_i) ;
}
#endif
