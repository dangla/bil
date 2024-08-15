#ifndef LOCALVARIABLES_H
#define LOCALVARIABLES_H


#include "LocalSpaceTime.h"
#include "Element.h"
#include "CustomValues.h"
#include "FEM.h"

template<typename T> struct LocalVariables_t;

#if 1
template<typename T>
using LocalVariables_SetInputs_t = T& (const LocalVariables_t<T>&,const int&,T&);
#endif

#if 0
template<typename T>
struct LocalVariables_SetInputs_t {
  T& operator()(const LocalVariables_t<T>&,const int&,T&);
} ;
#endif



template<typename T>
struct LocalVariables_t: virtual public LocalSpaceTime_t {
  private:
  const double* const* _u;
  const double* _f;
  /* This is a default internal memory */
  T  _savedlocalvalue[2];
  T* _localvalue;
  T* _localderivative;

  template<typename U>
  void _setvalues(const int& p,const double* f) {
    if(f) {
      U* val1 = (U*) f ;
      U& val2 = *_localvalue ;
      
      val2 = val1[p];
    }
  }
  
  template<typename U>
  void _storevalues(const int& p,double* f) {
    if(f) {
      U* val1 = (U*) f ;
      U& val2 = *_localvalue ;
      
      val1[p] = val2;
    }
  }
  
  public:
  /* Constructors */
  LocalVariables_t(const Element_t* el,const double& t,const double& dt,const double* const* u,const double* f):
  LocalSpaceTime_t(el,t,dt),_u(u),_f(f) {
    _localvalue = _savedlocalvalue;
    _localderivative = _savedlocalvalue + 1;
  }
  LocalVariables_t(const Element_t* el,const double& t,const double& dt,const double* const* u):
  LocalSpaceTime_t(el,t,dt),_u(u),_f(NULL) {
    _localvalue = _savedlocalvalue;
    _localderivative = _savedlocalvalue + 1;
  }
  //LocalVariables_t(void): LocalSpaceTime_t() {}
  
  /* Destructor */
  ~LocalVariables_t(void) {}
  
  /* Function call operator */
  //virtual T& operator()(const int&) {}
  
  /* Functions */
  template<typename F>
  T& SetInputs(F setin,const int& p) {
    return(setin(*this,p,*_localvalue));
  }

  void SetValues(const int& p) {
    {
      const double* f = GetImplicitTerm();
      using U = CustomValues_TypeOfImplicitValues(T);
      
      _setvalues<U>(p,f);
    }
  
    {
      const double* f = GetExplicitTerm();
      using U = CustomValues_TypeOfExplicitValues(T);
    
      _setvalues<U>(p,f);
    }
  
    {
      const double* f = GetConstantTerm();
      using U = CustomValues_TypeOfConstantValues(T);
      
      _setvalues<U>(p,f);
    }
  }
  
  void StoreImplicitTerms(const int& p,double* f) {
    using U = CustomValues_TypeOfImplicitValues(T);
      
    _storevalues<U>(p,f);
  }
  
  void StoreExplicitTerms(const int& p,double* f) {
    using U = CustomValues_TypeOfExplicitValues(T);
      
    _storevalues<U>(p,f);
  }

  void DisplacementVectorAndStrainFEM(const int& p,const int& u_mech,double* x) {
    const Element_t* el = GetElement();
    const double* const* u = _u; //GetPointerToNodalUnknowns();
    IntFct_t* intfct = Element_GetIntFct(el) ;
    FEM_t*    fem    = FEM_GetInstance(el) ;
    int dim = Element_GetDimensionOfSpace(el) ;
    
    /* Displacements */
    {
      for(int i = 0 ; i < dim ; i++) {
        x[i] = FEM_ComputeUnknown(fem,u,intfct,p,u_mech + i) ;
      }
    
      for(int i = dim ; i < 3 ; i++) {
        x[i] = 0 ;
      }
    }
    
    /* Strain */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,u_mech) ;
    
      for(int i = 0 ; i < 9 ; i++) {
        x[3 + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
  }

  void ValueAndGradientFEM(const int& p,const int& u_mass,double* x)
  {
    const Element_t* el = GetElement();
    const double* const* u = _u; //GetPointerToNodalUnknowns();
    IntFct_t* intfct = Element_GetIntFct(el) ;
    FEM_t*    fem    = FEM_GetInstance(el) ;
  
    /* Scalar */
    x[0] = FEM_ComputeUnknown(fem,u,intfct,p,u_mass) ;
    
    /* Scalar gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,u_mass) ;
    
      for(int i = 0 ; i < 3 ; i++) {
        x[1 + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  //void   SetPointerToNodalUnknowns(const double* const* u) {_u = u;}
  //void   SetImplicitTerm(const double* f) {_f = f;}
  void  SetLocalValue(T* v) {_localvalue = v;}
  void  ResetLocalValue(void) {_localvalue = _savedlocalvalue;}
  
  /* Accessors */
  const double* const*   GetPointerToNodalUnknowns(void) {return _u;}
  double*    GetImplicitTerm(void) {return _f;}
  double*    GetExplicitTerm(void) {
    if(Element_GetNbOfExplicitTerms(GetElement())) {
      return Element_GetExplicitTerm(GetElement());
    } else {
      return NULL;
    }
  }
  double*    GetConstantTerm(void) {
    if(Element_GetNbOfConstantTerms(GetElement())) {
      return Element_GetConstantTerm(GetElement());
    } else {
      return NULL;
    }
  }
  T&    GetLocalValue() {return *_localvalue;}
  T&    GetLocalDerivative() {return *_localderivative;}
  int   GetSizeOfLocalValues() {return (sizeof(T)/sizeof(double));}
} ;


#endif
