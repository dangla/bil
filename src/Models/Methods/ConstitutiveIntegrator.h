#ifndef CONSTITUTIVEINTEGRATOR_H
#define CONSTITUTIVEINTEGRATOR_H

#include "autodiff.h"

#include "Element.h"
#include "LocalVariables.h"
#include "LocalSpaceTime.h"

template<typename T>
using ConstitutiveIntegrator_Integrate_t = T& (const Element_t*,const double&,const double&,const T&,T&) ;

#if 0
template<typename T>
struct ConstitutiveIntegrator_Integrate_t {
  T&  operator()(const Element_t*,const double&,const double&,const T&,T&) ;
} ;
#endif


template <typename SETIN, typename INTEG>
struct ConstitutiveIntegrator_base {
  private:
  SETIN _setin;
  INTEG _integ;
  
  public:
  /* Constructors */
  ConstitutiveIntegrator_base(SETIN& setin,INTEG& integ): _setin(setin),_integ(integ) {}
  
  /* Destructor */
  ~ConstitutiveIntegrator_base(void) {}
  
  /* Accessors */
  SETIN& GetSetInputs(void) {return _setin;}
  INTEG& GetIntegrate(void) {return _integ;}
} ;



template <typename T,typename SETIN = LocalVariables_SetInputs_t<T>*, typename INTEG = ConstitutiveIntegrator_Integrate_t<T>*>
struct ConstitutiveIntegrator_t: public ConstitutiveIntegrator_base<SETIN,INTEG>, public LocalSpaceTime_t {
  private:
  LocalVariables_t<T> _var;
  LocalVariables_t<T> _var_n;
  
  public:
  /* Constructors */
  ConstitutiveIntegrator_t(const Element_t* el,const double& t,const double& dt,const double* const* u,const double* const* u_n,const double* f_n,SETIN setin,INTEG integ):
  ConstitutiveIntegrator_base<SETIN,INTEG>(setin,integ),LocalSpaceTime_t(el,t,dt),_var(el,t,dt,u),_var_n(el,t,dt,u_n,f_n) {}
  
  /* Destructor */
  ~ConstitutiveIntegrator_t(void) {}
  
  /* Functor */
  //template<typename T>
  //T* operator()(const double*,T*) ;
  
  /* Function call operator */
  //double* operator()(const int p) {}
  
  T& Integrate(const int& p) {
    const Element_t* el = GetElement();
    const double& t = GetCurrentTime();
    const double& dt = GetTimeIncrement();
    T& v = _var.GetLocalValue() ;
    T& v_n = _var_n.GetLocalValue() ;
    SETIN& setin = GetSetInputs();
    INTEG& integ = GetIntegrate();

    _var_n.SetValues(p);
    _var_n.SetInputs(setin,p);
    
    _var.SetInputs(setin,p);
    
    return(integ(el,t,dt,v_n,v));
  }
  
  void StoreImplicitTerms(const int& p,double* f) {
    _var.StoreImplicitTerms(p,f) ;
  }
  
  void StoreExplicitTerms(const int& p,double* f) {
    _var.StoreExplicitTerms(p,f) ;
  }
  
  T& Differentiate(const double& dxi,const int& i)
  {
    Element_t* el = GetElement();
    const double& t = GetCurrentTime();
    const double& dt = GetTimeIncrement();
    T& v = _var.GetLocalValue() ;
    T& dv = _var.GetLocalDerivative() ;
    T& v_n = _var_n.GetLocalValue() ;
    double* x = (double*) &v ;
    double* dx = (double*) &dv ;
    int nbofvariables = _var.GetSizeOfLocalValues() ;
    INTEG& integ = GetIntegrate();
    
    dv = v;
  
    /* We increment the variable as (x + dx) */
    dx[i] += dxi ;
  
    integ(el,t,dt,v_n,dv) ;
  
    /* The numerical derivative as (f(x + dx) - f(x))/dx */
    for(int j = 0 ; j < nbofvariables ; j++) {
      dx[j] -= x[j] ;
      dx[j] /= dxi ;
    }

    return(dv) ;
  }
  #if defined HAVE_AUTODIFF
  template<typename U> // U = CustomValues_t<I,E,C,real>
  T& AutoDifferentiate(const int& i)
  {
    Element_t* el = GetElement();
    const double& t = GetCurrentTime();
    const double& dt = GetTimeIncrement();
    T& v = _var.GetLocalValue() ;
    T& dv = _var.GetLocalDerivative() ;
    T& v_n = _var_n.GetLocalValue() ;
    double* x = (double*) &v ;
    double* dx = (double*) &dv ;
    int nbofvariables = _var.GetSizeOfLocalValues() ;
    INTEG& integ = GetIntegrate();
    U w ;
    real* r = (real*) w ;
  
    // w = v; /* Can we do this? */
    for(int j = 0 ; j < nbofvariables ; j++) {
      r[j] = x[j] ;
    }
  
    /* This is the partial derivative wrt i */
    r[i][1] = 1 ;
  
    integ(el,t,dt,v_n,w) ;
  
    for(int j = 0 ; j < nbofvariables ; j++) {
      dx[j] = r[j][1] ;
    }

    return(dv) ;
  }
  #endif
  
  void  SetLocalValue(double* v) {_var.SetLocalValue((T*) v);}
  void  ResetLocalValue(void) {_var.ResetLocalValue;}
  
  /* Accessors */
} ;



#endif
