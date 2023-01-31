#include <stdio.h>
#include <stdlib.h>
#include "Mry.h"
#include "Message.h"
#include "ResolutionMethod.h"


ResolutionMethod_t* (ResolutionMethod_Create)(Options_t* options)
{
  ResolutionMethod_t* rm = (ResolutionMethod_t*) Mry_New(ResolutionMethod_t) ;
  
  
  ResolutionMethod_GetOptions(rm) = options ;
  
  
  /*  Method */
  {
    if(Options_ResolutionMethodIsCrout(options)) {
      ResolutionMethod_GetType(rm) = ResolutionMethod_Type(CROUT) ;
    
    } else if(Options_ResolutionMethodIsSuperLU(options)) {
      ResolutionMethod_GetType(rm) = ResolutionMethod_Type(SuperLU) ;
    
    } else if(Options_ResolutionMethodIsSuperLUMT(options)) {
      ResolutionMethod_GetType(rm) = ResolutionMethod_Type(SuperLUMT) ;
    
    } else if(Options_ResolutionMethodIsSuperLUDist(options)) {
      ResolutionMethod_GetType(rm) = ResolutionMethod_Type(SuperLUDist) ;

    } else if(Options_ResolutionMethodIsMA38(options)) {
      ResolutionMethod_GetType(rm) = ResolutionMethod_Type(MA38) ;
      
    } else {
      arret("ResolutionMethod_Create: method not available") ;
    }
  }
  
  
  return(rm) ;
}



void (ResolutionMethod_Delete)(void* self)
{
  ResolutionMethod_t* rm = (ResolutionMethod_t*) self ;
}
