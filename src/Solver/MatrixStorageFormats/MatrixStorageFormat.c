#include <stdio.h>
#include <stdlib.h>
#include "Mry.h"
#include "Message.h"
#include "MatrixStorageFormat.h"


MatrixStorageFormat_t* (MatrixStorageFormat_Create)(Options_t* options)
{
  MatrixStorageFormat_t* msf = (MatrixStorageFormat_t*) Mry_New(MatrixStorageFormat_t) ;
  
  MatrixStorageFormat_GetOptions(msf) = options ;
  
  /*  Matrix storage format */
  {
    if(Options_ResolutionMethodIsCrout(options)) {
      MatrixStorageFormat_GetType(msf) = MatrixStorageFormat_Type(LDUSKL) ;
    
    } else if(Options_ResolutionMethodIsSuperLU(options)) {
      MatrixStorageFormat_GetType(msf) = MatrixStorageFormat_Type(SuperLU) ;
    
    } else if(Options_ResolutionMethodIsSuperLUMT(options)) {
      MatrixStorageFormat_GetType(msf) = MatrixStorageFormat_Type(SuperLU) ;
    
    } else if(Options_ResolutionMethodIsSuperLUDist(options)) {
      MatrixStorageFormat_GetType(msf) = MatrixStorageFormat_Type(SuperLU) ;

    } else if(Options_ResolutionMethodIsMA38(options)) {
      MatrixStorageFormat_GetType(msf) = MatrixStorageFormat_Type(Coordinate) ;
    
    } else {
      arret("MatrixStorageFormat_Create: method not available") ;
    }
  }
  
  return(msf) ;
}



void (MatrixStorageFormat_Delete)(void* self)
{
  MatrixStorageFormat_t* msf = (MatrixStorageFormat_t*) self ;
}
