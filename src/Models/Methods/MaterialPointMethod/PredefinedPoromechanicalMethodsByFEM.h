#ifndef PREDEFINEDPOROMECHANICALMETHODSBYFEM_H
#define PREDEFINEDPOROMECHANICALMETHODSBYFEM_H

#include "PredefinedMethods.h"


#include "BaseName.h"
#include "CustomValues.h"
#include "MaterialPointModel.h"



#define ImplicitValues_t BaseName(_ImplicitValues_t)
#define ExplicitValues_t BaseName(_ExplicitValues_t)
#define ConstantValues_t BaseName(_ConstantValues_t)
#define OtherValues_t    BaseName(_OtherValues_t)


template<typename T>
struct ImplicitValues_t ;

template<typename T>
struct ExplicitValues_t;

template<typename T>
struct ConstantValues_t;

template<typename T>
struct OtherValues_t;



template<typename T>
using Values_t = CustomValues_t<T,ImplicitValues_t,ExplicitValues_t,ConstantValues_t,OtherValues_t> ;

using Values_d = Values_t<double> ;

#define Values_Index(V)  CustomValues_Index(Values_d,V,double)


#define MPM_t      BaseName(_MPM_t)


struct MPM_t: public MaterialPointModel_t<Values_d> {
  MaterialPointModel_SetInputs_t<Values_d> SetInputs;
  MaterialPointModel_Integrate_t<Values_d> Integrate;
  MaterialPointModel_Initialize_t<Values_d>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_d> SetTangentMatrix;
  MaterialPointModel_SetTransferMatrix_t<Values_d> SetTransferMatrix;
  MaterialPointModel_SetIndexes_t SetIndexes;
  MaterialPointModel_SetIncrements_t SetIncrements;
} ;

#endif
