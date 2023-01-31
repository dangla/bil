#include <stdio.h>
#include <stdlib.h>
#include "Mry.h"
#include "Buffers.h"


Buffers_t* (Buffers_Create)(size_t size)
{
  Buffers_t* buffers = (Buffers_t*) Mry_New(Buffers_t) ;
  
  {
    const int n = Buffers_MaxNbOfBuffers ;
    
    Buffers_GetNbOfBuffers(buffers) = n ;
  
    if(n){
      Buffers_GetBuffer(buffers) = Mry_Create(Buffer_t,n,Buffer_Create(size)) ;
    }
  }

  return(buffers) ;
}



void (Buffers_Delete)(void* self)
{
  Buffers_t* buffers = (Buffers_t*) self ;
  
  {
    int n = Buffers_GetNbOfBuffers(buffers) ;
    Buffer_t* buffer = Buffers_GetBuffer(buffers) ;
    
    Mry_Delete(buffer,n,Buffer_Delete) ;
    
    if(buffer) free(buffer) ;
  }
}
