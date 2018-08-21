#include <string.h>
#include "Message.h"
#include "Views.h"


Views_t* Views_Create(int n)
{
  Views_t* views = (Views_t*) malloc(sizeof(Views_t)) ;
  
  if(!views) {
    arret("Views_Create") ;
  }
  
  Views_GetNbOfViews(views) = n ;
  
  Views_GetView(views) = View_Create(n) ;
  
  return(views);
}


void Views_Delete(void* self)
{
  Views_t** pviews = (Views_t**) self ;
  Views_t*   views = *pviews ;
  View_t* view = Views_GetView(views) ;
  
  View_Delete(&view) ;
  free(views) ;
  *pviews = NULL ;
}
