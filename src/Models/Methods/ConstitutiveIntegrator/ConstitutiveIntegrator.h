#ifndef CONSTITUTIVEINTEGRATOR_H
#define CONSTITUTIVEINTEGRATOR_H

#include "autodiff.h"

#include "Element.h"
#include "LocalVariables.h"
#include "LocalSpaceTime.h"

template<typename T>
using ConstitutiveIntegrator_SetInputs_t = T* (Element_t const*,double const&,double const&,double const* const*,int const&,T&);

template<typename T>
using ConstitutiveIntegrator_Integrate_t = T* (Element_t const*,double const&,double const&,const T&,T&) ;

template<typename T>
using ConstitutiveIntegrator_SetFluxes_t = T* (Element_t const*,int const&,int const&,T const&,T&,T&);

template<typename T>
using ConstitutiveIntegrator_Initialize_t = T* (Element_t const*,double const&,T&);



template <typename SETIN, typename INTEG>
struct ConstitutiveIntegrator_base {
  private:
  SETIN* _setin;
  INTEG* _integ;
  
  public:
  /* Constructors */
  ConstitutiveIntegrator_base(SETIN* setin,INTEG* integ): _setin(setin),_integ(integ) {}
  
  /* Destructor */
  //~ConstitutiveIntegrator_base(void) {}
  
  /* Accessors */
  SETIN& GetSetInputs(void) const {return *_setin;}
  INTEG& GetIntegrate(void) const {return *_integ;}
} ;



template <typename T,typename SETIN = ConstitutiveIntegrator_SetInputs_t<T>, typename INTEG = ConstitutiveIntegrator_Integrate_t<T>>
struct ConstitutiveIntegrator_t: public ConstitutiveIntegrator_base<SETIN,INTEG>, public LocalSpaceTime_t {
  private:
  LocalVariables_t<T> _var;
  LocalVariables_t<T> _var_n;
  
  public:
  /* Constructors */
  ConstitutiveIntegrator_t(Element_t const* el,double const& t,double const& dt,double const* const* u,double const* const* u_n,double const* f_n,SETIN* setin,INTEG* integ):
  ConstitutiveIntegrator_base<SETIN,INTEG>(setin,integ),LocalSpaceTime_t(el,t,dt),_var(el,t,dt,u,f_n),_var_n(el,t,dt,u_n,f_n) {}
  
  /* Destructor */
  //~ConstitutiveIntegrator_t(void) {}
  
  /* Functor */
  //template<typename T>
  //T* operator()(double const*,T*) ;
  
  /* Function call operator */
  //double* operator()(int const p) {}
  
  T* SetInputs(int const& p) {
    /* We must use "this->" because "GetSetInputs" should be looked-up here 
     * (it's a non-dependent name).
     */
    Element_t const* el = GetElement();
    double const& t = GetCurrentTime();
    double const& dt = GetTimeIncrement();
    double const* const* u = _var.GetPointerToNodalUnknowns();
    T* v = _var.GetLocalValue();
    SETIN& setin = this->GetSetInputs();
    
    return(setin(el,t,dt,u,p,v[p]));
  }
  
  T* SetValues(int const& p) {
    return(_var.SetValues(p));
  }
  
  T* Integrate(int const& p) {
    Element_t const* el = GetElement();
    double const& t = GetCurrentTime();
    double const& dt = GetTimeIncrement();
    double const* const* u = _var.GetPointerToNodalUnknowns();
    double const* const* u_n = _var_n.GetPointerToNodalUnknowns();
    T* v = _var.GetLocalValue();
    T* v_n = _var_n.GetLocalValue();
    SETIN& setin = this->GetSetInputs();
    INTEG& integ = this->GetIntegrate();

    _var_n.SetValues(p);
    
    setin(el,t,dt,u_n,p,v_n[p]);
    setin(el,t,dt,u,p,v[p]);
    
    return(integ(el,t,dt,v_n[p],v[p]));
  }
  
  void StoreImplicitTerms(int const& p,double* f) {
    _var.StoreImplicitTerms(p,f) ;
  }
  
  void StoreExplicitTerms(int const& p,double* f) {
    _var.StoreExplicitTerms(p,f) ;
  }
  
  void StoreConstantTerms(int const& p,double* f) {
    _var.StoreConstantTerms(p,f) ;
  }
  
  T* Differentiate(int const& p,double const& dxi,int const& i) {
    T* v = _var.GetLocalValue() ;
    T* dv = _var.GetLocalDerivative() ;
    double* dx = (double*) (dv + p) ;
    
    dv[p] = v[p];
  
    /* We increment the variable as (x + dx) */
    dx[i] += dxi ;
    
    {
      Element_t const* el = GetElement();
      double const& t = GetCurrentTime();
      double const& dt = GetTimeIncrement();
      T* v_n = _var_n.GetLocalValue() ;
      INTEG& integ = this->GetIntegrate();
      T* vo = integ(el,t,dt,v_n[p],dv[p]) ;
      
      if(!vo) return(NULL) ;
    }
    
    dv[p] -= v[p];
    dv[p] /= dxi;
    
    return(dv + p) ;
  }
  
  #if defined HAVE_AUTODIFF
  template<typename U> // U = CustomValues_t<I<real>,E<real>,C<real>>
  T* AutoDifferentiate(int const& p,int const& i) {
    T* v = _var.GetLocalValue() ;
    T* v_n = _var_n.GetLocalValue() ;
    double* x = (double*) (v + p) ;
    U w ;
    real* r = (real*) w ;
    int nbofvariables = _var.GetSizeOfLocalValues() ;
  
    // w = v; /* Can we do this? */
    {
      for(int j = 0 ; j < nbofvariables ; j++) {
        r[j] = x[j] ;
      }
    }
  
    /* This is the partial derivative wrt i */
    r[i][1] = 1 ;
  
    {
      Element_t* el = GetElement();
      double const& t = GetCurrentTime();
      double const& dt = GetTimeIncrement();
      INTEG& integ = this->GetIntegrate();
      T* vo = integ(el,t,dt,v_n[p],w) ;
      
      if(!vo) return(NULL) ;
    }
  
    {
      T* dv = _var.GetLocalDerivative() ;
      double* dx = (double*) (dv + p) ;
      
      for(int j = 0 ; j < nbofvariables ; j++) {
        dx[j] = r[j][1] ;
      }

      return(dv + p) ;
    }
  }
  #endif
  
  T* FiniteGradient(int const i,int const j) {
    T* dv = _var.GetLocalDerivative();
    
    return(_var.GradientFVM(i,j,dv[i]));
  }
  
  template<typename F>
  T* SetFluxes(F* fluxes,int const& i,int const& j,T const& grdv) {
    Element_t* el = GetElement();
    T* v = _var.GetLocalValue() ;
    return((*fluxes)(el,i,j,grdv,v[i],v[j]));
  }
  
  template<typename F>
  T* Initialize(F* init,int const& p) {
    T* v = SetInputs(p);
    
    if(v) {
      Element_t* el = GetElement();
      double const& t = GetCurrentTime();
      
      return((*init)(el,t,*v));
    } else {
      return(NULL);
    }
  }
  
  
  /* Accessors */
  //LocalVariables_t<T>& GetLocalVariables(void) {return(_var);}
  //T* GetLocalValue(void) {return(_var.GetLocalValue());}
} ;



#endif
