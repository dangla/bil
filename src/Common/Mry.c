#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Mry.h"
#include "Message.h"



void* Mry_Malloc(size_t size)
{
  void* ptr ;

  if(!size) return (NULL);
  
  ptr = malloc(size) ;
  
  assert(ptr) ;
  
  return (ptr);
}



void* Mry_Calloc(size_t num, size_t size)
{
  void* ptr;

  if(!size) return (NULL);
  
  ptr = calloc(num, size);
  
  assert(ptr) ;
  
  return (ptr);
}



void* Mry_Realloc(void* ptr, size_t size)
{
  if(!size) return (NULL);
  
  ptr = realloc(ptr, size);
  
  assert(ptr) ;
  
  return (ptr);
}



void Mry_Free(void* ptr)
{
  if(ptr == NULL) return;
  free(ptr);
}
