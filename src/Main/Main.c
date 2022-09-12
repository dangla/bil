#include <stdlib.h>
#include "Bil.h"


int main(int argc,char** argv)
{
  Bil_t* bil = Bil_Create(argc,argv) ;
  
  Bil_Execute(bil) ;
  
  Bil_Delete(bil) ;
  
  return(0) ;
}
