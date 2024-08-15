#ifndef LISTOFPLASTICITYMODELS_H
#define LISTOFPLASTICITYMODELS_H

#define ListOfPlasticityModels_Nb          8
#define ListOfPlasticityModels_Methods(m)  PlasticityACC##m,PlasticityACCPierre##m,PlasticityBBM##m,PlasticityBExM##m,PlasticityCamClay##m,PlasticityCamClayOffset##m,PlasticityDruckerPrager##m,PlasticityNSFS##m
#define ListOfPlasticityModels_Names       "PlasticityACC","PlasticityACCPierre","PlasticityBBM","PlasticityBExM","PlasticityCamClay","PlasticityCamClayOffset","PlasticityDruckerPrager","PlasticityNSFS"

#endif
