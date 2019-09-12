#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "Date.h"






Date_t*  Date_Create(int n_dates)
{
  Date_t* date = (Date_t*) malloc(n_dates*sizeof(Date_t)) ;
  
  if(!date) arret("Date_Create") ;
  
  return(date) ;
}

