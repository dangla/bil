#include "Session.h"
#include "GenericData.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



static Session_t* cursession = NULL ;



/* Global functions */
Session_t*  (Session_GetCurrentInstance)(void)
{
  assert(cursession) ;
  
  return(cursession) ;
}



Session_t*   (Session_New)(void)
{
  Session_t* session = (Session_t*) calloc(1,sizeof(Session_t)) ;
  
  assert(session) ;
  
  return(session) ;
}



void    (Session_Delete)(void)
{
  while(cursession) Session_Close() ;
}



Session_t*   (Session_Open)(void)
{
  Session_t* prev  = cursession ;
  
  cursession = Session_New() ;
  
  Session_GetPreviousSession(cursession) = prev ;
    
  if(prev) {
    Session_GetIndex(cursession) = Session_GetIndex(prev) + 1 ;
  } else {
    Session_GetIndex(cursession) = 1 ;
  }
  
  return(cursession) ;
}



Session_t*   (Session_Close)(void)
{
  GenericData_t* gdat = Session_GetGenericData(cursession) ;

  if(gdat) {
    GenericData_Delete(&gdat) ;
  }

  {
    Session_t* garbage = cursession ;
    
    cursession = Session_GetPreviousSession(garbage) ;

    free(garbage) ;
  }

  return(cursession) ;
}
