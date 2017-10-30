#include <stdlib.h>
#include "Message.h"
#include "Temperature.h"



static Temperature_t* instancetemperature = NULL ;
static Temperature_t* (Temperature_Create)(void) ;



Temperature_t* (Temperature_Create)(void)
{
  Temperature_t* temperature = (Temperature_t*) malloc(sizeof(Temperature_t)) ;
  
  if(!temperature) arret("Temperature_Create") ;
  
  {
    Temperature_GetRoomValue(temperature) = Temperature_DefaultValue ;
  }
  
  return(temperature) ;
}



Temperature_t* Temperature_GetInstance(void)
{
  if(!instancetemperature) {
    instancetemperature = Temperature_Create() ;
  }
  
  return(instancetemperature) ;
}
