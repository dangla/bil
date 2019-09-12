#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Message.h"
#include "ICond.h"



ICond_t* ICond_Create(int n_iconds)
{
  ICond_t* icond = (ICond_t*) malloc(n_iconds*sizeof(ICond_t)) ;
    
  if(!icond) arret("ICond_Create(1)") ;
    
    
  /* Allocation of space for the name of unknowns */
  {
    size_t sz = n_iconds*ICond_MaxLengthOfKeyWord*sizeof(char) ;
    char* name_unk = (char*) malloc(sz) ;
    int i ;
    
    if(!name_unk) arret("IConds_Create(2)") ;
  
    for(i = 0 ; i < n_iconds ; i++) {
      ICond_GetNameOfUnknown(icond + i) = name_unk + i*ICond_MaxLengthOfKeyWord ;
    }
  }
    
    
  /* Allocation of space for the name of files of nodal values */
  {
    size_t sz = n_iconds*ICond_MaxLengthOfFileName*sizeof(char) ;
    char* filename = (char*) malloc(sz) ;
    int i ;
    
    if(!filename) arret("IConds_Create(3)") ;
  
    for(i = 0 ; i < n_iconds ; i++) {
      ICond_GetFileNameOfNodalValues(icond + i) = filename + i*ICond_MaxLengthOfFileName ;
      ICond_GetFileNameOfNodalValues(icond + i)[0] = '\0' ;
    }
  }
  
  return(icond) ;
}
