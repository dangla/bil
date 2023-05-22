#include "Entry.h"


int main(int argc,char** argv)
{
  Entry_t* entry = Entry_Create(argc,argv) ;
  
  Entry_Execute(entry) ;
  
  Entry_Delete(entry) ;
  
  return(0) ;
}
