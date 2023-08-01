#include <stdio.h>
#include <stdlib.h>
#include "Mry.h"
#include "Message.h"
#include "Options.h"
#include "VectorStorageFormat.h"


VectorStorageFormat_t* (VectorStorageFormat_Create)(Options_t* options)
{
  VectorStorageFormat_t* vsf = (VectorStorageFormat_t*) Mry_New(VectorStorageFormat_t) ;
  
  VectorStorageFormat_GetOptions(vsf) = options ;
  
  /*  Vector storage format */
  {
    if(Options_ResolutionMethodIsCrout(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(Array) ;
    
    } else if(Options_ResolutionMethodIsSuperLU(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(Array) ;
    
    } else if(Options_ResolutionMethodIsSuperLUMT(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(Array) ;
    
    } else if(Options_ResolutionMethodIsSuperLUDist(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(Array) ;

    } else if(Options_ResolutionMethodIsMA38(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(Array) ;

    } else if(Options_ResolutionMethodIsPetscKSP(options)) {
      VectorStorageFormat_GetType(vsf) = VectorStorageFormat_Type(PetscVec) ;
    
    } else {
      arret("VectorStorageFormat_Create: method not available") ;
    }
  }
  
  return(vsf) ;
}



void (VectorStorageFormat_Delete)(void* self)
{
  VectorStorageFormat_t* vsf = (VectorStorageFormat_t*) self ;
}
