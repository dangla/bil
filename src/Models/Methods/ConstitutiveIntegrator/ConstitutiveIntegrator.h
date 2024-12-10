#ifndef CONSTITUTIVEINTEGRATOR_H
#define CONSTITUTIVEINTEGRATOR_H

#include "autodiff.h"

#include "Element.h"
#include "LocalVariables.h"
#include "LocalSpaceTime.h"
#include "MaterialPointModel.h"


template<typename T,typename MPM = MaterialPointModel_t<T>>
struct ConstitutiveIntegrator_t {
  public:
  /* Constructors */
  ConstitutiveIntegrator_t(MPM* mpm): _mpm(mpm) {}
  
  /* Destructor */
  ~ConstitutiveIntegrator_t(void) {}

  
  void Set(Element_t* el,double const& t,double const& dt,double const* const* u_n,double* f_n,double const* const* u,double* f) {
    _localspacetime = LocalSpaceTime_t(el,t,dt);
    _var   = LocalVariables_t<T>(u,f);
    _var_n = LocalVariables_t<T>(u_n,f_n);
  }

  T* ExtractInputs(int const& p) {
    Element_t* el = GetElement();
    double const& t = GetCurrentTime();
    double const* const* u = _var.GetPointerToNodalUnknowns();
    T* v = _var.GetLocalValue();
    T* vo = _mpm->SetInputs(el,t,p,u,v[p]);
    
    return(vo);
  }
  
  T* ExtractValues(int const& p) {
    Element_t* el = GetElement();
    T* vo = _var.Extract(el,p);
    
    return(vo) ;
  }
  
  T* IntegrateValues(int const& p) {
    Element_t* el = GetElement();
    double const& t = GetCurrentTime();
    double const& dt = GetTimeIncrement();
    double const* const* u = _var.GetPointerToNodalUnknowns();
    T* vp_n = _var_n.Extract(el,p);
    T* v    = _var.GetLocalValue();
    T* vp   = _mpm->SetInputs(el,t,p,u,v[p]);
    
    if(vp) {
      T* vo = _mpm->Integrate(el,t,dt,vp_n[0],vp[0]);
      
      return(vo);
    }
    
    return(NULL);
  }

  T* InitializeValues(int const& p) {
    Element_t* el = GetElement();
    double const& t = GetCurrentTime();
    double const* const* u = _var.GetPointerToNodalUnknowns();
    T* vp = _var.Extract(el,p);
    T* vo = _mpm->SetInputs(el,t,p,u,vp[0]);
    
    if(vo) {
      T* vi = _mpm->Initialize(el,t,vo[0]);
      
      return(vi);
    }
    
    return(NULL);
  }

  /*
   * Operations specific to FEM or FVM
   * ---------------------------------
   */
  int  ComputeImplicitTermsByFEM(void) {
    Element_t* el = GetElement();
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    
    for(int p = 0 ; p < NbOfIntPoints ; p++) {
      T* val = IntegrateValues(p) ;
      
      if(!val) return(1);
      
      StoreImplicitTerms(p) ;
    }

    return(0) ;
  }

  int  ComputeImplicitTermsByFVM(void) {
    Element_t* el = GetElement();
    int nn = Element_GetNbOfNodes(el) ;

    for(int i = 0 ; i < nn ; i++) {
      T* val = IntegrateValues(i) ;
      
      if(!val) return(1) ;
  
      if(Element_IsSubmanifold(el)) continue ;
      
      /* Flux */
      //if(setfl) 
      {
        for(int j = 0 ; j < i ; j++) {
          T* grdv = FiniteGradient(i,j);

          ComputeFluxes(i,j,grdv[0]) ;
        }
      }
    }

    {
      for(int i = 0 ; i < nn ; i++) {
        StoreImplicitTerms(i) ;
      }
    }

    return(0) ;
  }


  int  ComputeExplicitTermsByFEM(void) {
    Element_t* el = GetElement();
  
    if(Element_GetNbOfExplicitTerms(el)) {
      IntFct_t*  intfct = Element_GetIntFct(el) ;
      int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        T* val = IntegrateValues(p) ;
      
        if(!val) return(1);
    
        StoreExplicitTerms(p) ;
      }
    }
  
    return(0) ;
  }

  int  ComputeExplicitTermsByFVM(void) {
    Element_t* el = GetElement();
  
    if(Element_GetNbOfExplicitTerms(el)) {
      int nn = Element_GetNbOfNodes(el) ;
      
      for(int i = 0 ; i < nn ; i++) {
        T* val = IntegrateValues(i) ;
      
        if(!val) return(1) ;
        
        StoreExplicitTerms(i) ;
      }
    }

    return(0) ;
  }

  int ComputeInitialStateByFEM(void) {
    Element_t* el = GetElement();
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

    /* Pre-initialization */
    //if(initi) 
    {
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        T* val = InitializeValues(p);
      
        //if(!val) return(1);
        if(!val) break;
        
        StoreImplicitTerms(p) ;
        StoreConstantTerms(p) ;
      }
    }
    
    {
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        T* val = IntegrateValues(p) ;
      
        if(!val) return(1);
      
        StoreImplicitTerms(p) ;
        StoreExplicitTerms(p) ;
      }
    }

    return(0) ;
  }

  int ComputeInitialStateByFVM(void) {
    Element_t* el = GetElement();
    int nn = Element_GetNbOfNodes(el) ;  
  
    /* Pre-initialization */
    //if(initi) 
    {
      for(int i = 0 ; i < nn ; i++) {
        T* val = InitializeValues(i);
      
        //if(!val) return(1) ;
        if(!val) break ;
        
        StoreImplicitTerms(i) ;
        StoreConstantTerms(i) ;
      }
    }
  

    {
      for(int i = 0 ; i < nn ; i++) {
        T* val = IntegrateValues(i) ;
      
        if(!val) return(1) ;
        
        StoreExplicitTerms(i) ;
  
        if(Element_IsSubmanifold(el)) continue ;
      
        //if(setfl) 
        {
          for(int j = 0 ; j < i ; j++) {
            T* grdv = FiniteGradient(i,j);

            ComputeFluxes(i,j,grdv[0]);
          }
        }
      }
    }
  
    /* storages */
    {
      for(int i = 0 ; i < nn ; i++) {
        StoreImplicitTerms(i) ;
      }
    }
  
    return(0) ;
  }
  
  
  /* Compute matrices by FEM */
  double* ComputeTangentStiffnessMatrixByFEM(void) {
    Element_t* el = GetElement();
    int n = 9;
    double c[IntFct_MaxNbOfIntPoints*81] ;
    int dec = SetTangentMatrixForFEM(n,c);
      
    if(dec < 0) {
      return(NULL) ;
    } else {
      FEM_t* fem = FEM_GetInstance(el) ;
      IntFct_t*  intfct = Element_GetIntFct(el) ;
      double* kp = FEM_ComputeStiffnessMatrix(fem,intfct,c,dec) ;
    
      return(kp);
    }
  }
  
  double* ComputeConductionMatrixByFEM(void) {
    Element_t* el = GetElement();
    int neq = Element_GetNbOfEquations(el);
    int ndif = neq;
    int const n = 3*ndif;
    double c[IntFct_MaxNbOfIntPoints*n*n] ;    
    int dec = SetTransferMatrixForFEM(n,c);
      
    if(dec < 0) {
      return(NULL) ;
    } else {
      IntFct_t*  intfct = Element_GetIntFct(el) ;
      FEM_t* fem = FEM_GetInstance(el) ;
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec,ndif) ;
    
      return(kc);
    }
  }

  double* ComputePoromechanicalMatrixByFEM(int const& e_mech) {
    Element_t* el = GetElement();
    int dim = Element_GetDimensionOfSpace(el);
    int neq = Element_GetNbOfEquations(el);
    int ndif = neq - dim;
    int ncols = 9 + ndif;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* kp;
    double* kc;
    
    {
      int const n = ncols;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;
      int dec = SetTangentMatrixForFEM(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,ndif,e_mech) ;
      }
    }
    
    {
      int const n = 3*ndif;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;    
      int dec = SetTransferMatrixForFEM(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec,ndif) ;
      }
    }
    
    if(kp) {
      if(kc) {
        double* k = FEM_AssemblePoroelasticAndConductionMatrices(fem,ndif,e_mech,kc,kp);
        
        return(k);
      } else {
        return(kp);
      }
    } else {
      return(kc);
    }
  }

  #if defined HAVE_AUTODIFF
  template<typename U> // U = CustomValues_t<I<real>,E<real>,C<real>>
  double* ComputeAutodiffPoromechanicalMatrixByFEM(int const& e_mech) {
    Element_t* el = GetElement();
    int dim = Element_GetDimensionOfSpace(el);
    int neq = Element_GetNbOfEquations(el);
    int ndif = neq - dim;
    int ncols = 9 + ndif;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* kp;
    double* kc;
    
    {
      int const n = ncols;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;
      int dec = SetAutodiffTangentMatrixForFEM<U>(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,ndif,e_mech) ;
      }
    }
    
    {
      int const n = 3*ndif;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;    
      int dec = SetTransferMatrixForFEM(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec,ndif) ;
      }
    }
    
    if(kp) {
      if(kc) {
        double* k = FEM_AssemblePoroelasticAndConductionMatrices(fem,ndif,e_mech,kc,kp);
        
        return(k);
      } else {
        return(kp);
      }
    } else {
      return(kc);
    }
  }
  #endif

  double* ComputeMassConservationMatrixByFEM(void) {
    Element_t* el = GetElement();
    int neq = Element_GetNbOfEquations(el);
    int ndif = neq;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* kp;
    double* kc;
    
    {
      int const n = ndif;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;
      int dec = SetTangentMatrixForFEM(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kp = FEM_ComputeMassMatrix(fem,intfct,c,dec,ndif) ;
      }
    }
    
    {
      int const n = 3*ndif;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;    
      int dec = SetTransferMatrixForFEM(n,c);
      
      if(dec < 0) {
        return(NULL) ;
      } else {
        kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec,ndif) ;
      }
    }
    
    if(kp) {
      if(kc) {
        int nn  = Element_GetNbOfNodes(el);
        int ndof = nn*neq;
        
        for(int i = 0 ; i < ndof*ndof ; i++) {
          kp[i] += kc[i];
        }
        
        FEM_FreeBufferFrom(fem,kc);
      }
      
      return(kp);
    } else {
      return(kc);
    }
  }

  /* Compute matrices by FVM */
  double* ComputeMassConservationMatrixByFVM(void) {
    Element_t* el = GetElement();
    int nn  = Element_GetNbOfNodes(el);
    int neq = Element_GetNbOfEquations(el);
    int ndof = nn*neq;
    double c[ndof*ndof] ;
    int dec = SetTangentMatrixForFVM(neq,c);
    
    if(dec < 0) {
      return(NULL) ;
    } else {
      FVM_t* fvm = FVM_GetInstance(el) ;
      double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,neq) ;
      
      return(km);
    }
  }

  /* Compute residus by FEM */
  double* ComputeStrainWorkResiduByFEM(int const& index) {
    Element_t* el = GetElement();
    double* vim = Element_GetCurrentImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    double* x = (double*) val ;
    double* stress = x + index;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,stress,nvi) ;
    
    return(rw);
  }
  
  double* ComputeBodyForceResiduByFEM(int const& index) {
    Element_t* el = GetElement();
    double* vim = Element_GetCurrentImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    double* x = (double*) val ;
    double* bforce = x + index;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,bforce,nvi) ;
    
    return(rbf);
  }
  
  double* ComputeMechanicalEquilibiumResiduByFEM(int const& istress,int const& ibforce) {
    Element_t* el = GetElement();
    double* vim = Element_GetCurrentImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    double* x = (double*) val ;
    double* stress = x + istress;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* rw = NULL;
    
    if(istress >= 0) {
      rw = FEM_ComputeStrainWorkResidu(fem,intfct,stress,nvi) ;
    } else {
      int dim = Element_GetDimensionOfSpace(el) ;
      int nn = Element_GetNbOfNodes(el) ;
      int ndof = nn*dim ;
      size_t SizeNeeded = ndof*(sizeof(double)) ;
      rw = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
      
      for(int i = 0 ; i < ndof ; i++) {      
        rw[i] = 0;
      }
    }


    if(ibforce >= 0) {
      double* bforce = x + ibforce;
      int dim = Element_GetDimensionOfSpace(el);
      int nn = Element_GetNbOfNodes(el) ;
      
      for(int j = 0 ; j < dim ; j++) {        
        if(bforce[j] != 0) {
          double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,bforce + j,nvi) ;
          
          for(int i = 0 ; i < nn ; i++) {      
            rw[i*dim + j] += -rbf[i];
          }

          FEM_FreeBufferFrom(fem,rbf);
        }
      }
    }
    
    return(rw);
  }
  
  double* ComputeMassConservationResiduByFEM(int const& imass,int const& iflow) {
    Element_t* el = GetElement();
    double dt = GetTimeIncrement();
    double* vi = Element_GetCurrentImplicitTerm(el) ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vi ;
    TI* val_n  = (TI*) vi_n ;
    double* x = (double*) val ;
    double* x_n = (double*) val_n ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np = IntFct_GetNbOfFunctions(intfct);
    int nn = Element_GetNbOfNodes(el);
    double* rm = NULL;
    
    if(imass >= 0) {
      double* mass = x + imass;
      double* mass_n = x_n + imass;
      double g1[IntFct_MaxNbOfIntPoints] ;
    
      for(int i = 0 ; i < np ; i++) {
        g1[i] = mass[i*nvi] - mass_n[i*nvi] ;
      }
    
      rm = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    } else {
      size_t SizeNeeded = nn*(sizeof(double)) ;
      rm = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
      
      for(int i = 0 ; i < nn ; i++) {
        rm[i] = 0 ;
      }
    }
  
    if(iflow >= 0) {
      double* flow = x + iflow;
      double* rf = FEM_ComputeFluxResidu(fem,intfct,flow,nvi) ;
      
      for(int i = 0 ; i < nn ; i++) {
        rm[i] += -dt*rf[i] ;
      }

      FEM_FreeBufferFrom(fem,rf);
    }
    
    return(rm);
  }
  
  /* Compute residus by FVM */
  double* ComputeBodyForceResiduByFVM(int const& index) {
    Element_t* el = GetElement();
    double* vim = Element_GetCurrentImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    double* x = (double*) val ;
    FVM_t* fvm = FVM_GetInstance(el) ;
    int nn = Element_GetNbOfNodes(el);
    double* rbf;
    
    if(index >= 0) {
      double* bforce = x + index;
      
      rbf = FVM_ComputeBodyForceResidu(fvm,bforce,nvi) ;
    } else {
      size_t SizeNeeded = nn*(sizeof(double)) ;
      rbf = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
      
      for(int i = 0 ; i < nn ; i++) {
        rbf[i] = 0 ;
      }
    }
    
    return(rbf);
  }
  
  double*  ComputeFluxResiduByFVM(int const& iflow) {
    Element_t* el = GetElement();
    double* vi = Element_GetCurrentImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vi ;
    double* x = (double*) val ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    FVM_t* fvm = FVM_GetInstance(el) ;
    int nn = Element_GetNbOfNodes(el);
    double* rf;       
    
    if(iflow >= 0) {
      double* flow = x + iflow;
      
      rf = FVM_ComputeFluxResidu(fvm,flow,nvi) ;
    } else {
      size_t SizeNeeded = nn*(sizeof(double)) ;
      rf = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
      
      for(int i = 0 ; i < nn ; i++) {
        rf[i] = 0 ;
      }
    }
      
    return(rf);
  }


  double* ComputeMassConservationResiduByFVM(int const& imass,int const& iflow) {
    Element_t* el = GetElement();
    double dt = GetTimeIncrement();
    double* vi = Element_GetCurrentImplicitTerm(el) ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    using TI = CustomValues_TypeOfImplicitValues(T);
    TI* val  = (TI*) vi ;
    TI* val_n  = (TI*) vi_n ;
    double* x = (double*) val ;
    double* x_n = (double*) val_n ;
    int nvi = CustomValues_NbOfImplicitValues(T);
    FVM_t* fvm = FVM_GetInstance(el) ;
    int nn = Element_GetNbOfNodes(el);    
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(int i = 0 ; i < nn ; i++) {
      for(int j = 0 ; j < nn ; j++) {
        if(i == j) {
          if(imass >= 0) {
            double* mass = x + imass;
            double* mass_n = x_n + imass;
            
            g[i*nn + i] = mass[i*nvi] - mass_n[i*nvi] ;
          } else {
            g[i*nn + i] = 0 ;
          }
        } else {
          if(iflow >= 0) {
            double* flow = x + iflow;
            
            g[i*nn + j] = dt * flow[i*nvi + j] ;
          } else {
            g[i*nn + j] = 0;
          }
        }
      }
    }
    
    {
      double* rm = FVM_ComputeMassAndFluxResidu(fvm,g,nn) ;
      
      return(rm);
    }
  }
  
  
  private:
  LocalSpaceTime_t _localspacetime;
  LocalVariables_t<T> _var;
  LocalVariables_t<T> _var_n;
  MPM* _mpm;
  
  
  /* Accessors */
  Element_t*     GetElement(void) const {return _localspacetime.GetElement();}
  double const&  GetCurrentTime(void) const {
    //Element_t* el = GetElement();
    //double const& t = Element_GetTime(el);
    
    //return t;
    return _localspacetime.GetCurrentTime();
  }
  double const&  GetTimeIncrement(void) const {
    //Element_t* el = GetElement();
    //double const& dt = Element_GetTimeStep(el);
    
    //return dt;
    return _localspacetime.GetTimeIncrement();
  }

  #if 0
  T* SetInputs(Element_t* el,double const& t,double const* const* u,int const& p,T& v) {
    return(_mpm->SetInputs(el,t,p,u,v));
  }

  template<typename U>
  U* Integrate(Element_t* el,double const& t,double const& dt,T const& v_n,U& v) {
    return(_mpm->Integrate(el,t,dt,v_n,v));
  }

  T* Initialize(Element_t* el,double const& t,T& v) {
    return(_mpm->Initialize(el,t,v));
  }

  int SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,T const& v,T const& dv,int const& k,double* c) {
    return(_mpm->SetTangentMatrix(el,t,dt,p,v,dv,k,c));
  }

  int SetTransferMatrix(Element_t* el,double const& dt,int const& p,T const& v,double* c) {
    return(_mpm->SetTransferMatrix(el,dt,p,v,c));
  }

  T* SetFluxes(Element_t* el,double const& t,int const& i,int const& j,T const& grd,T* v) {
    return(_mpm->SetFluxes(el,t,i,j,grd,v));
  }

  int* SetIndexes(Element_t* el,int* ind) {
    _mpm->SetIndexes(el,ind);
    return(ind);
  }

  double* SetIncrements(Element_t* el,double* dui) {
    _mpm->SetIncrements(el,dui);
    return(dui);
  }
  #endif
  
  
  void StoreImplicitTerms(int const& p) {
    _var.StoreImplicitTerms(p) ;
  }
  
  void StoreExplicitTerms(int const& p) {
    Element_t* el = GetElement();
    
    _var.StoreExplicitTerms(el,p) ;
  }
  
  void StoreConstantTerms(int const& p) {
    Element_t* el = GetElement();
    
    _var.StoreConstantTerms(el,p) ;
  }
  
  
  T* DifferentiateValues(int const& p,double const& dxi,int const& i) {
    T* v = _var.GetLocalValue() ;
    T* dv = _var.GetLocalDerivative() ;
    double* dx = (double*) (dv + p) ;
    
    dv[p] = v[p];
  
    /* We increment the variable as (x + dx) */
    dx[i] += dxi ;
    
    {
      Element_t* el = GetElement();
      double const& t = GetCurrentTime();
      double const& dt = GetTimeIncrement();
      T* v_n = _var_n.GetLocalValue() ;
      T* vo = _mpm->Integrate(el,t,dt,v_n[p],dv[p]) ;
      
      if(!vo) return(NULL) ;
    }
    
    dv[p] -= v[p];
    dv[p] /= dxi;
    
    return(dv + p) ;
  }
  
  /* Below:
   * T = CustomValues_t<C<double>...>
   * U = CustomValues_t<C<real>...>
   */
  #if defined HAVE_AUTODIFF
  template<typename U>
  T* AutoDifferentiateValues(int const& p,int const& i) {
    T* v = _var.GetLocalValue() ;
    T* v_n = _var_n.GetLocalValue() ;
    U w ;
    real* r ;
    
    if(!CustomValues_IsValueType(U,real)) {
      Message_FatalError("AutoDifferentiateValues: the type is not \"real\"!");
    }
    
    r = (real*) &w ;
    
    #if 1
    w = v[p];
    #else
    {
      double* x = (double*) (v + p) ;
      int nbofvariables = _var.GetSizeOfLocalValues() ;
      
      for(int j = 0 ; j < nbofvariables ; j++) {
        r[j] = x[j] ;
      }
    }
    #endif
    
    /* This is the partial derivative wrt i */
    r[i][1] = 1 ;
  
    {
      Element_t* el = GetElement();
      double const& t = GetCurrentTime();
      double const& dt = GetTimeIncrement();
      U* vo = _mpm->Integrate(el,t,dt,v_n[p],w) ;
      
      if(!vo) return(NULL) ;
    }
  
    {
      T* dv = _var.GetLocalDerivative() ;
      double* dx = (double*) (dv + p) ;
      int nbofvariables = _var.GetSizeOfLocalValues() ;
      
      for(int j = 0 ; j < nbofvariables ; j++) {
        dx[j] = r[j][1] ;
      }

      return(dv + p) ;
    }
  }
  #endif
  
  T* FiniteGradient(int const i,int const j) {
    Element_t* el = GetElement();
    T* dv = _var.GetLocalDerivative();
    
    return(_var.GradientFVM(el,i,j,dv[i]));
  }
  
  T* ComputeFluxes(int const& i,int const& j,T const& grdv) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    T* v = _var.GetLocalValue() ;
    T* vo = _mpm->SetFluxes(el,t,i,j,grdv,v);
      
    return(vo);
  }
  
  T* ComputeFluxDerivatives(int const& i,int const& j,T const& grdv) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    T* dv = _var.GetLocalDerivative() ;
    T* dvo = _mpm->SetFluxes(el,t,i,j,grdv,dv);
      
    return(dvo);
  }


  int SetTangentMatrixForFEM(int const& ncols,double* c) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    double dt = GetTimeIncrement();
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int    ind[9+Model_MaxNbOfEquations];
    double dui[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(el,ind);
    _mpm->SetIncrements(el,dui);
    
    {
      int dec = ncols*ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }

    for(int p = 0 ; p < np ; p++) {
      //T* val = IntegrateValues(p) ;
      T* val   = _var.Extract(el,p);
      T* val_n = _var_n.Extract(el,p);
      
      if(!val) return(-1);
            
      for(int k = 0 ; k < ncols ; k++) {
        double dui_k = (dui) ? dui[k] : 0;
        int    ind_k = (ind) ? ind[k] : -1;
        
        if(ind_k >= 0) {
          T* dval = DifferentiateValues(p,dui_k,ind_k) ;
        
          shift = _mpm->SetTangentMatrix(el,t,dt,p,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(el,t,dt,p,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }


  #if defined HAVE_AUTODIFF
  template<typename U>
  int SetAutodiffTangentMatrixForFEM(int const& ncols,double* c) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    double dt = GetTimeIncrement();
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int    ind[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(el,ind);
    
    {
      int dec = ncols*ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }

    for(int p = 0 ; p < np ; p++) {
      //T* val = IntegrateValues(p) ;
      T* val   = _var.Extract(el,p);
      T* val_n = _var_n.Extract(el,p);
      
      if(!val) return(-1);
            
      for(int k = 0 ; k < ncols ; k++) {
        int    ind_k = (ind) ? ind[k] : -1;
        
        if(ind_k >= 0) {
          T* dval = AutoDifferentiateValues<U>(p,ind_k) ;
        
          shift = _mpm->SetTangentMatrix(el,t,dt,p,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(el,t,dt,p,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  #endif
  

  int SetTangentMatrixForFVM(int const& ncols,double* c) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    double dt = GetTimeIncrement();
    int nn = Element_GetNbOfNodes(el) ;
    int    ind[9+Model_MaxNbOfEquations];
    double dui[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(el,ind);
    _mpm->SetIncrements(el,dui);
    
    {
      int  dec = ncols*ncols ;
      
      for(int i = 0 ; i < nn*nn*dec ; i++) c[i] = 0 ;
    }
  
    for(int i = 0 ; i < nn ; i++) {
      //T* val = IntegrateValues(i) ;
      T* val   = _var.Extract(el,i);
      T* val_n = _var_n.Extract(el,i);
    
      if(!val) return(-1) ;
    
      for(int k = 0 ; k < ncols ; k++) {
        double  dui_k = dui[k] ;
        int     ind_k = ind[k];
        
        if(ind_k >= 0) {
          T* dval = DifferentiateValues(i,dui_k,ind_k) ;
      
          for(int j = 0 ; j < nn ; j++) {
            T* dv = ComputeFluxDerivatives(i,j,dval[0]) ;
          }
      
          shift = _mpm->SetTangentMatrix(el,t,dt,i,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(el,t,dt,i,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  

  #if defined HAVE_AUTODIFF
  template<typename U>
  int SetAutodiffTangentMatrixForFVM(int const& ncols,double* c) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    double dt = GetTimeIncrement();
    int nn = Element_GetNbOfNodes(el) ;
    int ind[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(el,ind);
    
    {
      int  dec = ncols*ncols ;
      
      for(int i = 0 ; i < nn*nn*dec ; i++) c[i] = 0 ;
    }
  
    for(int i = 0 ; i < nn ; i++) {
      //T* val = IntegrateValues(i) ;
      T* val   = _var.Extract(el,i);
      T* val_n = _var_n.Extract(el,i);
    
      if(!val) return(-1) ;
    
      for(int k = 0 ; k < ncols ; k++) {
        int     ind_k = ind[k];
        
        if(ind_k >= 0) {
          T* dval = AutoDifferentiateValues<U>(i,ind_k) ;
      
          for(int j = 0 ; j < nn ; j++) {
            T* dv = ComputeFluxDerivatives(i,j,dval[0]) ;
          }
      
          shift = _mpm->SetTangentMatrix(el,t,dt,i,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(el,t,dt,i,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  #endif
  
  
  int SetTransferMatrixForFEM(int const& ncols,double* c) {
    Element_t* el = GetElement();
    double t  = GetCurrentTime();
    double dt = GetTimeIncrement();
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int shift;
    
    /* initialization */
    {
      int dec = ncols * ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }
  
    for(int p = 0 ; p < np ; p++) {
      T* val = _var.Extract(el,p);
    
      shift = _mpm->SetTransferMatrix(el,dt,p,val[0],c);
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
} ;


#undef SETIN
#undef INTEG
#undef SETFL
#undef INITI
#undef SETGM
#undef SETRM

#endif
