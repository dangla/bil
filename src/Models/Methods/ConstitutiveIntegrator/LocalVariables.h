#ifndef LOCALVARIABLES_H
#define LOCALVARIABLES_H


#include "LocalSpaceTime.h"
#include "Math_.h"
#include "Element.h"
#include "IntFct.h"
#include "CustomValues.h"
#include "FEM.h"
#include "FVM.h"


#define LocalVariables_MaxNbOfLocalValues   (Math_Max(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints)+1)



template<typename T>
struct LocalVariables_t: virtual public LocalSpaceTime_t {
  private:
  double const* const* _u;
  double const* _f;
  T  _localvalue[LocalVariables_MaxNbOfLocalValues];
  T  _localderivative[LocalVariables_MaxNbOfLocalValues];

  template<typename U>
  void _setvalues(int const& p,double const* f) {
    if(f) {
      U const* val1 = (U const*) f ;
      U* val2 = _localvalue + p ;
      
      val2[0] = val1[p];
    }
  }
  
  template<typename U>
  void _storevalues(int const& p,double* f) const {
    if(f) {
      U* val1 = (U*) f ;
      U* val2 = _localvalue + p ;
      
      val1[p] = val2[0];
    }
  }
  
  public:
  /* Constructors */
  LocalVariables_t(Element_t const* el,double const& t,double const& dt,double const* const* u,double const* f): LocalSpaceTime_t(el,t,dt),_u(u),_f(f) {}
  
  /* Destructor */
  //~LocalVariables_t(void) {}
  
  /* Function call operator */
  //virtual T& operator()(int const&) {}
  
  /* Functions */
  #if 0
  template<typename F>
  T* SetInputs(F setin,int const& p) {
    if(p < LocalVariables_MaxNbOfLocalValues) {
      T* v = setin(*this,p,_localvalue[p]);
      return(v);
    }
    return(NULL);
  }
  #endif

  T* SetValues(int const& p) {
    if(p < LocalVariables_MaxNbOfLocalValues) {
      {
        double const* f = GetImplicitTerm();
        using U = CustomValues_TypeOfImplicitValues(T);
      
        _setvalues<U>(p,f);
      }
  
      {
        double const* f = GetExplicitTerm();
        using U = CustomValues_TypeOfExplicitValues(T);
    
        _setvalues<U>(p,f);
      }
  
      {
        double const* f = GetConstantTerm();
        using U = CustomValues_TypeOfConstantValues(T);
      
        _setvalues<U>(p,f);
      }
    
      return(_localvalue + p);
    }
    
    return(NULL);
  }
  
  void StoreImplicitTerms(int const& p,double* f) const {
    if(p < LocalVariables_MaxNbOfLocalValues) {
      using U = CustomValues_TypeOfImplicitValues(T);
      
      _storevalues<U>(p,f);
    }
  }
  
  void StoreExplicitTerms(int const& p,double* f) const {
    if(p < LocalVariables_MaxNbOfLocalValues) {
      using U = CustomValues_TypeOfExplicitValues(T);
      
      _storevalues<U>(p,f);
    }
  }
  
  void StoreConstantTerms(int const& p,double* f) const {
    if(p < LocalVariables_MaxNbOfLocalValues) {
      using U = CustomValues_TypeOfConstantValues(T);
      
      _storevalues<U>(p,f);
    }
  }

  double* DisplacementVectorAndStrainFEM(int const& p,int const& u_mech,double* x) const {
    Element_t const* el = GetElement();
    double const* const* u = _u; //GetPointerToNodalUnknowns();
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
    
    return(x);
  }

  double* ValueAndGradientFEM(int const& p,int const& u_mass,double* x) const {
    Element_t const* el = GetElement();
    double const* const* u = _u; //GetPointerToNodalUnknowns();
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
    
    return(x);
  }

  #if 0
  double* ValueFVM(int const& n,int const& u_mass,double* x) const {
    Element_t const* el = GetElement();
    double const* const* u = _u;
    //IntFct_t* intfct = Element_GetIntFct(el) ;
  
    /* Scalar */
    x[0] = Element_GetValueOfNodalUnknown(el,u,n,u_mass);
    
    return(x);
  }
  #endif
  
  T* GradientFVM(int const& i,int const& j,T& grdv) const {
    Element_t const* el = GetElement();
    FVM_t* fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    int nn = Element_GetNbOfNodes(el) ;
    double dij  = dist[nn*i + j] ;
    T const* val = _localvalue;
    
    if(i < LocalVariables_MaxNbOfLocalValues
    && j < LocalVariables_MaxNbOfLocalValues) {
      grdv = (val[j] - val[i])/dij ;
      
      return(&grdv);
    }
    
    return(NULL);
  }
  
  //void   SetPointerToNodalUnknowns(double const* const* u) {_u = u;}
  //void   SetImplicitTerm(double const* f) {_f = f;}
  //void  SetLocalValue(T* v) {_localvalue = v;}
  //void  ResetLocalValue(void) {_localvalue = _savedlocalvalue;}
  
  /* Accessors */
  double const* const*   GetPointerToNodalUnknowns(void) const {return _u;}
  double const*    GetImplicitTerm(void) const {return _f;}
  double const*    GetExplicitTerm(void) const {
    if(Element_GetNbOfExplicitTerms(GetElement())) {
      return Element_GetExplicitTerm(GetElement());
    } else {
      return NULL;
    }
  }
  double const*    GetConstantTerm(void) const {
    if(Element_GetNbOfConstantTerms(GetElement())) {
      return Element_GetConstantTerm(GetElement());
    } else {
      return NULL;
    }
  }
  T*    GetLocalValue() {return _localvalue;}
  T*    GetLocalDerivative() {return _localderivative;}
  int   GetSizeOfLocalValues() const {return (sizeof(T)/sizeof(double));}
} ;


#endif
