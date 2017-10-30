#ifndef PREDEFINEDMETHODS_H
#define PREDEFINEDMETHODS_H

static Module_SetModuleProp_t           SetModuleProp ;

#include "BaseName_SetModuleProp.h"

extern Module_SetModuleProp_t BaseName_SetModuleProp ;

int BaseName_SetModuleProp(Module_t* module)
{
  return(SetModuleProp(module));
}

#endif
