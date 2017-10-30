#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "TypeId.h"
#include "GenericData.h"




/* Global functions */

GenericData_t* (GenericData_New)(void)
{
  GenericData_t* gdat = (GenericData_t*) malloc(sizeof(GenericData_t)) ;
  
  if(!gdat) {
    arret("GenericData_New(1)") ;
  }
  
  {
    GenericData_GetTypeId(gdat) = TypeId_Create(undefined) ;
    GenericData_GetNbOfData(gdat) = 0 ;
    GenericData_GetData(gdat) = NULL ;
    GenericData_GetNextGenericData(gdat) = NULL ;
    GenericData_GetPreviousGenericData(gdat) = NULL ;
  }
  
  return(gdat) ;
}



void (GenericData_Delete)(GenericData_t** pgdat)
{
  GenericData_t* gdat = *pgdat ;
  GenericData_t* prev = GenericData_GetPreviousGenericData(gdat) ;
  GenericData_t* next = GenericData_GetNextGenericData(gdat) ;
  
  /* Connect prev and next */
  if(prev) GenericData_GetNextGenericData(prev) = next ;
  if(next) GenericData_GetPreviousGenericData(next) = prev ;
  
  {
    void* data = GenericData_GetData(gdat) ;
      
    free(data) ;
  }
  
  free(gdat) ;
}



void (GenericData_Initialize_)(GenericData_t* gdat,int n,void* data,TypeId_t typ)
{
  {
    GenericData_GetTypeId(gdat) = typ ;
    GenericData_GetNbOfData(gdat) = n ;
    GenericData_GetData(gdat) = data ;
  }
}



GenericData_t* (GenericData_Create_)(int n,size_t sz,TypeId_t typ)
{
  GenericData_t* gdat = GenericData_New() ;
  void* data = malloc(n * sz) ;
  
  if(!data) {
    arret("GenericData_Create_(1)") ;
  }
  
  GenericData_Initialize_(gdat,n,data,typ) ;
  
  
  return(gdat) ;
}



void (GenericData_InsertBefore)(GenericData_t* cur,GenericData_t* gdat)
/** Insert "gdat" in the linked list before "cur". */
{
  
  if(!cur) {
    arret("GenericData_InsertBefore(1)") ;
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
    arret("GenericData_InsertAfter(1)") ;
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
