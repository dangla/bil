#ifndef PLASTICITYMODELS_H
#define PLASTICITYMODELS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


#include "ListOfPlasticityModels.h"


#define PlasticityModels_NbOfModels               (ListOfPlasticityModels_Nb)
#define PlasticityModels_ListOfNames               ListOfPlasticityModels_Names
#define PlasticityModels_ListOfSetModelProp        ListOfPlasticityModels_Methods(_SetModelProp)


#ifdef __CPLUSPLUS
}
#endif
#endif
