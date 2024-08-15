#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "Mry.h"
#include "Buffer.h"
#include "Message.h"


Buffer_t* (Buffer_Create)(size_t size)
{
  Buffer_t* buffer = (Buffer_t*) Mry_New(Buffer_t) ;
  
  Buffer_GetBeginOfBuffer(buffer) = (void*) Mry_Allocate(size) ;
  
  Buffer_GetEndOfBuffer(buffer) = (char*) Buffer_GetBeginOfBuffer(buffer) + size ;
  Buffer_GetSize(buffer) = size ;
  
  Buffer_Free(buffer) ;
  
  return(buffer) ;
}



void (Buffer_Free)(Buffer_t* buffer)
{
  Buffer_GetHeadOfBuffer(buffer)  = Buffer_GetBeginOfBuffer(buffer) ;
  Buffer_GetTailOfBuffer(buffer)  = Buffer_GetHeadOfBuffer(buffer) ;
  Buffer_GetAvailableSize(buffer) = Buffer_GetSize(buffer) ;
}



void (Buffer_Delete)(void* self)
{
  Buffer_t* buffer = (Buffer_t*) self ;
  
  {
    void* v = Buffer_GetBeginOfBuffer(buffer) ;
  
    if(v) {
      free(v) ;
      Buffer_GetBeginOfBuffer(buffer) = NULL ;
    }
  }
}



void* (Buffer_Allocate)(Buffer_t* buffer,size_t sz)
{
  char* head = (char*) Buffer_GetHeadOfBuffer(buffer) ;
  char* end  = (char*) Buffer_GetEndOfBuffer(buffer) ;
  
  if(head + sz > end) {
    //size_t dead_sz = end - head ;
    long int dead_sz = end - head ;
    Buffer_GetAvailableSize(buffer) -= dead_sz ;
    head = (char*) Buffer_GetBeginOfBuffer(buffer) ;
  }
  
  if(sz > Buffer_GetAvailableSize(buffer)) {
    arret("Buffer_Allocate: not enough memory") ;
  }
  
  Buffer_GetHeadOfBuffer(buffer)  = head + sz ;
  Buffer_GetAvailableSize(buffer) -= sz ;

  return((void*) head) ;
}



void (Buffer_FreeFrom)(Buffer_t* buffer,char* p)
{
  char* begin = (char*) Buffer_GetBeginOfBuffer(buffer) ;
  char* end   = (char*) Buffer_GetEndOfBuffer(buffer) ;
  char* head  = (char*) Buffer_GetHeadOfBuffer(buffer) ;
  char* tail  = (char*) Buffer_GetTailOfBuffer(buffer) ;
  
  /* p must lie between the begin and the end of buffer */
  if(p < begin || p >= end) {
    arret("Buffer_FreeFrom: pointer out of range") ;
  }
  
  /* Case 1: the free zone is after the head OR before the tail */
  if(head >= tail) {
    /* If p lies into the free zone then stop */
    if(p >= head || p < tail) {
      arret("Buffer_FreeFrom(1): pointer out of range") ;
      
    } else {
      size_t free_sz = head - p ;
      Buffer_GetHeadOfBuffer(buffer) = p ;
      Buffer_GetAvailableSize(buffer) += free_sz ; 
    }
    
  /* Case 2: the free zone is after the head AND before the tail */
  } else {
    /* If p lies into the free zone then stop */
    if(p >= head && p < tail) {
      arret("Buffer_FreeFrom(2): pointer out of range") ;
      
    } else if(p < head) {
      size_t free_sz = head - p ;
      Buffer_GetHeadOfBuffer(buffer) = p ;
      Buffer_GetAvailableSize(buffer) += free_sz ;
      
    } else if(p >= tail) {
      size_t busy_sz = p - tail ;
      Buffer_GetHeadOfBuffer(buffer) = p ;
      Buffer_GetAvailableSize(buffer) = Buffer_GetSize(buffer) - busy_sz ;
    }
  }
}

