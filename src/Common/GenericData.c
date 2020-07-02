#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Message.h"
#include "TypeId.h"
#include "GenericData.h"
//#include "GenericObject.h"
#include "Tools/Math.h"
#include "Mry.h"


static void           (GenericData_Remove)(GenericData_t**) ;


/* Global functions */

GenericData_t* (GenericData_New)(void)
{
  GenericData_t* gdat = (GenericData_t*) Mry_New(GenericData_t) ;
  
  /* Allocation for the name */
  {
    char* name = (char*) Mry_New(char[GenericData_MaxLengthOfKeyWord + 1]) ;
    
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




GenericData_t* (GenericData_Create_)(int n,void* data,TypeId_t typ,const char* name)
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



void (GenericData_Initialize_)(GenericData_t* gdat,int n,void* data,TypeId_t typ,const char* name)
{
  GenericData_GetTypeId(gdat) = typ ;
  GenericData_GetNbOfData(gdat) = n ;
  GenericData_GetData(gdat) = data ;
    
  {
    int len = MIN(strlen(name),GenericData_MaxLengthOfKeyWord) ;
      
    strncpy(GenericData_GetName(gdat),name,len) ;
    GenericData_GetName(gdat)[len] = '\0' ;
  }
}



GenericData_t* (GenericData_Append)(GenericData_t* a,GenericData_t* b)
/** Append "b" to "a". 
 *  Return the appended generic data (ie "a" or "b") or NULL pointer. */
{
  
  if(!a) {
    return(b) ;
  }
    
  /* Add b at the end of a:
   * a ... last(a) - b */
  {
    GenericData_t* lasta  = GenericData_Last(a) ;
    
    /* The condition is useless since "lasta" must not be NULL.  */
    if(lasta) GenericData_GetNextGenericData(lasta) = b ;
    if(b) GenericData_GetPreviousGenericData(b) = lasta ;
    
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




GenericData_t* (GenericData_Find_)(GenericData_t* gdat,TypeId_t typ,const char* name)
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



#if 0
void (GenericData_InsertBefore)(GenericData_t* a,GenericData_t* b)
/** Insert "b" in the linked list before "a". */
{
  
  if(!a) {
    exit(0) ;
    return ;
  }
  
  /* Insert b between prev(a) and a: prev(a) - b - a */
  {
    GenericData_t* preva = GenericData_GetPreviousGenericData(a) ;
    
    GenericData_GetPreviousGenericData(b) = preva ;
    GenericData_GetNextGenericData(b)     = a ;
    if(preva) GenericData_GetNextGenericData(preva) = b ;
    if(a) GenericData_GetPreviousGenericData(a) = b ;
  }
}




void (GenericData_InsertAfter)(GenericData_t* a,GenericData_t* b)
/** Insert "b" in the linked list after "a". */
{
  
  if(!a) {
    exit(0) ;
    return ;
  }
    
  /* Insert b between a and next(a): a - b - next(a) */
  {
    GenericData_t* nexta = GenericData_GetNextGenericData(a) ;
    
    GenericData_GetPreviousGenericData(b) = a ;
    GenericData_GetNextGenericData(b)     = nexta ;
    if(nexta) GenericData_GetPreviousGenericData(nexta) = b ;
    if(a) GenericData_GetNextGenericData(a) = b ;
  }
}
#endif
