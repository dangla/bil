#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "Views.h"


Views_t* (Views_Create)(int n)
{
  Views_t* views = (Views_t*) Mry_New(Views_t) ;
  
  Views_GetNbOfViews(views) = n ;
  
  {
    View_t* view = (View_t*) Mry_New(View_t[n]) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      View_t* vw = View_Create() ;
  
      view[i] = vw[0] ;
      free(vw) ;
    }
  
    Views_GetView(views) = view ;
  }
  
  return(views);
}



void (Views_Delete)(void* self)
{
  Views_t* views = (Views_t*) self ;
  
  {
    int n = Views_GetNbOfViews(views) ;
    View_t* view = Views_GetView(views) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      View_Delete(view + i) ;
    }
    
    free(view) ;
    Views_GetView(views) = NULL ;
  }
}
