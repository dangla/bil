#ifndef PREDEFINEDMODULEMETHODS_H
#define PREDEFINEDMODULEMETHODS_H

static Module_SetModuleProp_t           SetModuleProp ;

#include "BaseName.h"

#define BaseName_SetModuleProp  BaseName(_SetModuleProp)

extern Module_SetModuleProp_t BaseName_SetModuleProp ;

int BaseName_SetModuleProp(Module_t* module)
{
  return(SetModuleProp(module));
}

#endif
