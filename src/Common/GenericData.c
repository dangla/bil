#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Message.h"
#include "TypeId.h"
#include "GenericData.h"
#include "GenericObject.h"
#include "Tools/Math.h"




/* Global functions */

GenericData_t* (GenericData_New)(void)
{
  GenericData_t* gdat = (GenericData_t*) malloc(sizeof(GenericData_t)) ;
  
  assert(gdat) ;
  
  /* Allocation for the name */
  {
    size_t sz = (GenericData_MaxLengthOfKeyWord + 1)*sizeof(char) ;
    char* name = (char*) malloc(sz) ;
    
    assert(name) ;
    
    GenericData_GetName(gdat) = name ;
  }
  
  {
    GenericData_GetTypeId(gdat) = TypeId_Create(undefined) ;
    strcpy(GenericData_GetName(gdat),"\0") ;
    GenericData_GetNbOfData(gdat) = 0 ;
    GenericData_GetData(gdat) = NULL ;
    GenericData_GetNextGenericData(gdat) = NULL ;
    GenericData_GetPreviousGenericData(gdat) = NULL ;
  }
  
  GenericData_GetDelete(gdat) = GenericData_Delete ;
  
  return(gdat) ;
}




GenericData_t* (GenericData_Create_)(int n,void* data,TypeId_t typ,char const* name)
{
  GenericData_t* gdat = GenericData_New() ;
  
  GenericData_Initialize_(gdat,n,data,typ,name) ;
  
  return(gdat) ;
}




void (GenericData_Delete)(void* self)
{
  GenericData_t** pgdat = (GenericData_t**) self ;
  
  while(*pgdat) GenericData_Remove(pgdat) ;
}



void (GenericData_Remove)(GenericData_t** pgdat)
{
  GenericData_t* gdat = *pgdat ;
  GenericData_t* prev = GenericData_GetPreviousGenericData(gdat) ;
  GenericData_t* next = GenericData_GetNextGenericData(gdat) ;
  
  /* Connect prev and next */
  if(prev) GenericData_GetNextGenericData(prev) = next ;
  if(next) GenericData_GetPreviousGenericData(next) = prev ;
  
  {
    TypeId_t typ = GenericData_GetTypeId(gdat) ;
    void* data = GenericData_GetData(gdat) ;
    
    TypeId_Delete(typ,&data) ;
    //GenericObject_Delete(&data) ;
  }
  
  free(gdat) ;
  
  /* *pgdat set to NULL if the content is empty */
  if(prev) {
    *pgdat = prev ;
  } else {
    *pgdat = next ;
  }
}



void (GenericData_Initialize_)(GenericData_t* gdat,int n,void* data,TypeId_t typ,char const* name)
{
  GenericData_GetTypeId(gdat) = typ ;
  GenericData_GetNbOfData(gdat) = n ;
  GenericData_GetData(gdat) = data ;
    
  {
    int len = MAX(strlen(name),GenericData_MaxLengthOfKeyWord) ;
      
    strncpy(GenericData_GetName(gdat),name,len) ;
    GenericData_GetName(gdat)[len] = '\0' ;
  }
}



void (GenericData_InsertBefore)(GenericData_t* cur,GenericData_t* gdat)
/** Insert "gdat" in the linked list before "cur". */
{
  
  if(!cur) {
    exit ;
    return ;
  }
  
  /* Insert gdat between prev and cur: prev - gdat - cur */
  {
    GenericData_t* prev = GenericData_GetPreviousGenericData(cur) ;
    
    GenericData_GetPreviousGenericData(gdat) = prev ;
    GenericData_GetNextGenericData(gdat)     = cur ;
    if(prev) GenericData_GetNextGenericData(prev) = gdat ;
    if(cur) GenericData_GetPreviousGenericData(cur) = gdat ;
  }
}



void (GenericData_InsertAfter)(GenericData_t* cur,GenericData_t* gdat)
/** Insert "gdat" in the linked list after "cur". */
{
  
  if(!cur) {
    exit ;
    return ;
  }
    
  /* Insert gdat between cur and next: cur - gdat - next */
  {
    GenericData_t* next = GenericData_GetNextGenericData(cur) ;
    
    GenericData_GetPreviousGenericData(gdat) = cur ;
    GenericData_GetNextGenericData(gdat)     = next ;
    if(next) GenericData_GetPreviousGenericData(next) = gdat ;
    if(cur) GenericData_GetNextGenericData(cur) = gdat ;
  }
}



GenericData_t* (GenericData_Merge)(GenericData_t* a,GenericData_t* b)
/** Merge a and b into one generic data and return it. */
{
  
  if(!a) {
    return(b) ;
  }
    
  /* Insert the complete b between a and next(a): 
   * a - first(b) ... last(b) - next(a) */
  if(a) {
    GenericData_t* next  = GenericData_GetNextGenericData(a) ;
    GenericData_t* first = GenericData_First(b) ;
    GenericData_t* last  = GenericData_Last(b) ;
    
    GenericData_GetNextGenericData(a) = b ;
    GenericData_GetPreviousGenericData(first) = a ;
    
    GenericData_GetNextGenericData(last) = next ;
    if(next) GenericData_GetPreviousGenericData(next) = b ;
    
    return(a) ;
  }
}



GenericData_t* (GenericData_Last)(GenericData_t* gdat)
/** Return the last generic data or NULL pointer. */
{
  if(gdat) {
    GenericData_t* next ;
    
    while((next = GenericData_GetNextGenericData(gdat))) gdat = next ;
  }
  
  return(gdat) ;
}



GenericData_t* (GenericData_First)(GenericData_t* gdat)
/** Return the first generic data or NULL pointer. */
{
  if(gdat) {
    GenericData_t* prev ;
    
    while((prev = GenericData_GetPreviousGenericData(gdat))) gdat = prev ;
  }
  
  return(gdat) ;
}




GenericData_t* (GenericData_Find_)(GenericData_t* gdat,TypeId_t typ,char const* name)
/** Return the generic data named as "name" or NULL pointer. */
{
  {
    GenericData_t* next = gdat ;
    
    while(next) {
      if(GenericData_Is(next,typ,name)) return(next) ;
      next = GenericData_GetNextGenericData(next) ;
    }
  }
  
  {
    GenericData_t* prev = gdat ;
    
    while(prev) {
      if(GenericData_Is(prev,typ,name)) return(prev) ;
      prev = GenericData_GetPreviousGenericData(prev) ;
    }
  }
  
  return(NULL) ;
}
