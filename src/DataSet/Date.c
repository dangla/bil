#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "Date.h"






Date_t*  (Date_New)(void)
{
  Date_t* date = (Date_t*) Mry_New(Date_t) ;
  
  return(date) ;
}



void  (Date_Delete)(void* self)
{
  Date_t* date = (Date_t*) self ;
}
