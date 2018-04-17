#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Message.h"
#include "Context.h"
#include "Module.h"



void* Module_Initialize(void*) ;



Module_t* Module_Create(int n_module)
/* Create modules */
{
  Module_t* module = (Module_t*) calloc(n_module,sizeof(Module_t)) ;
  
  assert(module) ;
  
  {
    int i ;
    
    for(i = 0 ; i < n_module ; i++) {
      Module_t* module_i = module + i ;
    
      Module_Initialize(module_i) ;
    }
  }
  
  return(module) ;
}



void* Module_Initialize(void* self)
{
  Module_t* module_i = (Module_t*) self ;
  
  {
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
  
  return(module_i) ;
}
