#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "View.h"



View_t* (View_Create)(void)
{
  View_t* view = (View_t*) Mry_New(View_t) ;

  /* Allocate memory for the name */
  {
    char* text = (char*) Mry_New(char[View_MaxLengthOfViewName]) ;
    
    View_GetNameOfView(view) = text ;
  }
  
  return(view) ;
}



void (View_Delete)(void* self)
{
  View_t* view = (View_t*) self ;
  
  free(View_GetNameOfView(view)) ;
}
