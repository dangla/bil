//#include <stdlib.h>
#include "Mry.h"
#include "Temperature.h"



Temperature_t* (Temperature_Create)(void)
{
  Temperature_t* temperature = (Temperature_t*) Mry_New(Temperature_t) ;
  
  {
    Temperature_GetRoomValue(temperature) = Temperature_DefaultValue ;
  }
  
  return(temperature) ;
}



void (Temperature_Delete)(void* self)
{
  Temperature_t* temperature = (Temperature_t*) self ;
}
