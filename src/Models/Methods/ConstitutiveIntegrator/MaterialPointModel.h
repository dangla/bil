#ifndef MATERIALPOINTMODEL_H
#define MATERIALPOINTMODEL_H

#include "Element.h"
#include "Model.h"



template<typename V>
using MaterialPointModel_SetInputs_t = V* (Element_t*,double const&,int const&,double const* const*,V&);

template<typename V,typename U = V>
using MaterialPointModel_Integrate_t = U* (Element_t*,double const&,double const&,V const&,U&) ;

template<typename V>
using MaterialPointModel_Initialize_t = V* (Element_t*,double const&,V&);

template<typename V>
using MaterialPointModel_SetTangentMatrix_t = int (Element_t*,double const&,double const&,int const&,V const&,V const&,int const&,double*);

template<typename V>
using MaterialPointModel_SetTransferMatrix_t = int (Element_t*,double const&,int const&,V const&,double*);

template<typename V>
using MaterialPointModel_SetFluxes_t = V* (Element_t*,double const&,int const&,int const&,V const&,V*);

using MaterialPointModel_SetIndexes_t = void (Element_t*,int*);

using MaterialPointModel_SetIncrements_t = void (Element_t*,double*);



template<typename V,typename U = V>
struct MaterialPointModel_t {
  public:
  
  MaterialPointModel_t() {}
  virtual ~MaterialPointModel_t() {}
  
  virtual V* SetInputs(Element_t*,double const&,int const&,double const* const*,V&) {return(NULL);}

  virtual U* Integrate(Element_t*,double const&,double const&,V const&,U&) {return(NULL);}

  virtual V* Initialize(Element_t*,double const&,V&) {return(NULL);}

  virtual int SetTangentMatrix(Element_t*,double const&,double const&,int const&,V const&,V const&,int const&,double*) {return(-1);}

  virtual int SetTransferMatrix(Element_t*,double const&,int const&,V const&,double*) {return(-1);}

  virtual V* SetFluxes(Element_t*,double const&,int const&,int const&,V const&,V*) {return(NULL);}

  virtual void SetIndexes(Element_t*,int*) {}

  virtual void SetIncrements(Element_t*,double*) {}
};

#endif
