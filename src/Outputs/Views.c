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


void Views_Delete(Views_t** views)
{
  View_t* view = Views_GetView(*views) ;
  
  View_Delete(&view) ;
  free(*views) ;
}


View_t* View_Create(int n)
{
  View_t* view = (View_t*) malloc(n*sizeof(View_t)) ;
  
  if(!view) {
    arret("View_Create(1)") ;
  }
  
  
  /* Allocate memory for the names */
  {
    char* text = (char*) malloc(n*View_MaxLengthOfViewName*sizeof(char)) ;
    int i ;
    
    if(!text) {
      arret("View_Create(3)") ;
    }
    
    for(i = 0 ; i < n ; i++) {
      View_GetNameOfView(view + i) = text + View_MaxLengthOfViewName*i ;
    }
  }
  
  return(view) ;
}


void View_Delete(View_t** view)
{
  free(View_GetNameOfView(*view)) ;
  free(*view) ;
}
