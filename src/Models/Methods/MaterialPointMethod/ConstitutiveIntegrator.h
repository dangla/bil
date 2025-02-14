#ifndef CONSTITUTIVEINTEGRATOR_H
#define CONSTITUTIVEINTEGRATOR_H

#include "autodiff.h"

#include "Element.h"
#include "LocalVariables.h"


template<typename MPM>
struct ConstitutiveIntegrator_t {
  private:
  template<typename T> using V = typename MPM::template Value_type<T>;
  
  public:
  /* Constructors */
  ConstitutiveIntegrator_t(MPM* mpm): _mpm(mpm) {}
  
  /* Destructor */
  ~ConstitutiveIntegrator_t(void) {}

  
  void Set(Element_t* el,double const& t,double const& dt,double const* const* u_n,double* f_n,double const* const* u,double* f) {
    _el = el;
    _t  = t;
    _dt = dt;
    _var   = LocalVariables_t<V<double>>(u,f);
    _var_n = LocalVariables_t<V<double>>(u_n,f_n);
  }

  V<double>* ExtractInputs(int const& p) {
    double const* const* u = _var.GetPointerToNodalUnknowns();
    V<double>* v = _var.GetLocalValue();
    V<double>* vo = _mpm->SetInputs(_el,_t,p,u,v[p]);
    
    return(vo);
  }
  
  V<double>* ExtractValues(int const& p) {
    V<double>* vo = _var.Extract(_el,p);
    
    return(vo) ;
  }
  
  V<double>* IntegrateValues(int const& p) {
    double const* const* u = _var.GetPointerToNodalUnknowns();
    V<double>* vp_n = _var_n.Extract(_el,p);
    V<double>* v    = _var.GetLocalValue();
    V<double>* vp   = _mpm->SetInputs(_el,_t,p,u,v[p]);
    
    if(vp) {
      V<double>* vo = _mpm->Integrate(_el,_t,_dt,vp_n[0],vp[0]);
      
      return(vo);
    }
    
    return(NULL);
  }

  V<double>* InitializeValues(int const& p) {
    double const* const* u = _var.GetPointerToNodalUnknowns();
    V<double>* vp = _var.Extract(_el,p);
    V<double>* vo = _mpm->SetInputs(_el,_t,p,u,vp[0]);
    
    if(vo) {
      V<double>* vi = _mpm->Initialize(_el,_t,vo[0]);
      
      return(vi);
    }
    
    return(NULL);
  }

  /*
   * Operations specific to FEM or FVM
   * ---------------------------------
   */
  int  ComputeImplicitTermsByFEM(void) {
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    
    for(int p = 0 ; p < NbOfIntPoints ; p++) {
      V<double>* val = IntegrateValues(p) ;
      
      if(!val) return(1);
      
      StoreImplicitTerms(p) ;
    }

    return(0) ;
  }

  int  ComputeImplicitTermsByFVM(void) {
    int nn = Element_GetNbOfNodes(_el) ;

    for(int i = 0 ; i < nn ; i++) {
      V<double>* val = IntegrateValues(i) ;
      
      if(!val) return(1) ;
  
      if(Element_IsSubmanifold(_el)) continue ;
      
      /* Flux */
      {
        for(int j = 0 ; j < i ; j++) {
          V<double>* grdv = FiniteGradient(i,j);

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
    if(Element_GetNbOfExplicitTerms(_el)) {
      IntFct_t*  intfct = Element_GetIntFct(_el) ;
      int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        V<double>* val = IntegrateValues(p) ;
      
        if(!val) return(1);
    
        StoreExplicitTerms(p) ;
      }
    }
  
    return(0) ;
  }

  int  ComputeExplicitTermsByFVM(void) {  
    if(Element_GetNbOfExplicitTerms(_el)) {
      int nn = Element_GetNbOfNodes(_el) ;
      
      for(int i = 0 ; i < nn ; i++) {
        V<double>* val = IntegrateValues(i) ;
      
        if(!val) return(1) ;
        
        StoreExplicitTerms(i) ;
      }
    }

    return(0) ;
  }

  int ComputeInitialStateByFEM(void) {
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

    /* Pre-initialization */
    {
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        V<double>* val = InitializeValues(p);
      
        if(!val) break;
        
        StoreImplicitTerms(p) ;
        StoreConstantTerms(p) ;
      }
    }
    
    {
      for(int p = 0 ; p < NbOfIntPoints ; p++) {
        V<double>* val = IntegrateValues(p) ;
      
        if(!val) return(1);
      
        StoreImplicitTerms(p) ;
        StoreExplicitTerms(p) ;
      }
    }

    return(0) ;
  }

  int ComputeInitialStateByFVM(void) {
    int nn = Element_GetNbOfNodes(_el) ;  
  
    /* Pre-initialization */
    {
      for(int i = 0 ; i < nn ; i++) {
        V<double>* val = InitializeValues(i);
      
        if(!val) break ;
        
        StoreImplicitTerms(i) ;
        StoreConstantTerms(i) ;
      }
    }
  

    {
      for(int i = 0 ; i < nn ; i++) {
        V<double>* val = IntegrateValues(i) ;
      
        if(!val) return(1) ;
        
        StoreExplicitTerms(i) ;
  
        if(Element_IsSubmanifold(_el)) continue ;
      
        {
          for(int j = 0 ; j < i ; j++) {
            V<double>* grdv = FiniteGradient(i,j);

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
    int n = 9;
    double c[IntFct_MaxNbOfIntPoints*81] ;
    int dec = SetTangentMatrixForFEM(n,c);
      
    if(dec < 0) {
      return(NULL) ;
    } else {
      FEM_t* fem = FEM_GetInstance(_el) ;
      IntFct_t*  intfct = Element_GetIntFct(_el) ;
      double* kp = FEM_ComputeStiffnessMatrix(fem,intfct,c,dec) ;
    
      return(kp);
    }
  }
  
  double* ComputeConductionMatrixByFEM(void) {
    int neq = Element_GetNbOfEquations(_el);
    int ndif = neq;
    int const n = 3*ndif;
    double c[IntFct_MaxNbOfIntPoints*n*n] ;    
    int dec = SetTransferMatrixForFEM(n,c);
      
    if(dec < 0) {
      return(NULL) ;
    } else {
      IntFct_t*  intfct = Element_GetIntFct(_el) ;
      FEM_t* fem = FEM_GetInstance(_el) ;
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec,ndif) ;
    
      return(kc);
    }
  }

  double* ComputePoromechanicalMatrixByFEM(int const& e_mech) {
    int dim = Element_GetDimensionOfSpace(_el);
    int neq = Element_GetNbOfEquations(_el);
    int ndif = neq - dim;
    int ncols = 9 + ndif;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
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

  double* ComputeMassConservationMatrixByFEM(void) {
    int neq = Element_GetNbOfEquations(_el);
    int ndif = neq;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
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
        int nn  = Element_GetNbOfNodes(_el);
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

  #if defined HAVE_AUTODIFF
  double* ComputeAutodiffPoromechanicalMatrixByFEM(int const& e_mech) {
    int dim = Element_GetDimensionOfSpace(_el);
    int neq = Element_GetNbOfEquations(_el);
    int ndif = neq - dim;
    int ncols = 9 + ndif;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    double* kp;
    double* kc;
    
    {
      int const n = ncols;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;
      int dec = SetAutodiffTangentMatrixForFEM(n,c);
      
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

  double* ComputeAutodiffMassConservationMatrixByFEM(void) {
    int neq = Element_GetNbOfEquations(_el);
    int ndif = neq;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    double* kp;
    double* kc;
    
    {
      int const n = ndif;
      double c[IntFct_MaxNbOfIntPoints*n*n] ;
      int dec = SetAutodiffTangentMatrixForFEM(n,c);
      
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
        int nn  = Element_GetNbOfNodes(_el);
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
  #endif

  /* Compute matrices by FVM */
  double* ComputeMassConservationMatrixByFVM(void) {
    int nn  = Element_GetNbOfNodes(_el);
    int neq = Element_GetNbOfEquations(_el);
    int ndof = nn*neq;
    double c[ndof*ndof] ;
    int dec = SetTangentMatrixForFVM(neq,c);
    
    if(dec < 0) {
      return(NULL) ;
    } else {
      FVM_t* fvm = FVM_GetInstance(_el) ;
      double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,neq) ;
      
      return(km);
    }
  }
  
  #if defined HAVE_AUTODIFF
  double* ComputeAutodiffMassConservationMatrixByFVM(void) {
    int nn  = Element_GetNbOfNodes(_el);
    int neq = Element_GetNbOfEquations(_el);
    int ndof = nn*neq;
    double c[ndof*ndof] ;
    int dec = SetAutodiffTangentMatrixForFVM(neq,c);
    
    if(dec < 0) {
      return(NULL) ;
    } else {
      FVM_t* fvm = FVM_GetInstance(_el) ;
      double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,neq) ;
      
      return(km);
    }
  }
  #endif

  /* Compute residus by FEM */
  double* ComputeStrainWorkResiduByFEM(int const& index) {
    double* vim = Element_GetCurrentImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    double* x = (double*) val ;
    double* stress = x + index;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,stress,nvi) ;
    
    return(rw);
  }
  
  double* ComputeBodyForceResiduByFEM(int const& index) {
    double* vim = Element_GetCurrentImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    double* x = (double*) val ;
    double* bforce = x + index;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,bforce,nvi) ;
    
    return(rbf);
  }
  
  double* ComputeMechanicalEquilibiumResiduByFEM(int const& istress,int const& ibforce) {
    double* vim = Element_GetCurrentImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    double* x = (double*) val ;
    double* stress = x + istress;
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    double* rw = NULL;
    
    if(istress >= 0) {
      rw = FEM_ComputeStrainWorkResidu(fem,intfct,stress,nvi) ;
    } else {
      int dim = Element_GetDimensionOfSpace(_el) ;
      int nn = Element_GetNbOfNodes(_el) ;
      int ndof = nn*dim ;
      size_t SizeNeeded = ndof*(sizeof(double)) ;
      rw = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
      
      for(int i = 0 ; i < ndof ; i++) {      
        rw[i] = 0;
      }
    }


    if(ibforce >= 0) {
      double* bforce = x + ibforce;
      int dim = Element_GetDimensionOfSpace(_el);
      int nn = Element_GetNbOfNodes(_el) ;
      
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
    double* vi = Element_GetCurrentImplicitTerm(_el) ;
    double* vi_n = Element_GetPreviousImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vi ;
    TI* val_n  = (TI*) vi_n ;
    double* x = (double*) val ;
    double* x_n = (double*) val_n ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    FEM_t* fem = FEM_GetInstance(_el) ;
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int np = IntFct_GetNbOfFunctions(intfct);
    int nn = Element_GetNbOfNodes(_el);
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
        rm[i] += -_dt*rf[i] ;
      }

      FEM_FreeBufferFrom(fem,rf);
    }
    
    return(rm);
  }
  
  /* Compute residus by FVM */
  double* ComputeBodyForceResiduByFVM(int const& index) {
    double* vim = Element_GetCurrentImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vim ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    double* x = (double*) val ;
    FVM_t* fvm = FVM_GetInstance(_el) ;
    int nn = Element_GetNbOfNodes(_el);
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
    double* vi = Element_GetCurrentImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vi ;
    double* x = (double*) val ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    FVM_t* fvm = FVM_GetInstance(_el) ;
    int nn = Element_GetNbOfNodes(_el);
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
    double* vi = Element_GetCurrentImplicitTerm(_el) ;
    double* vi_n = Element_GetPreviousImplicitTerm(_el) ;
    using TI = CustomValues_TypeOfImplicitValues(V<double>);
    TI* val  = (TI*) vi ;
    TI* val_n  = (TI*) vi_n ;
    double* x = (double*) val ;
    double* x_n = (double*) val_n ;
    int nvi = CustomValues_NbOfImplicitValues(V<double>);
    FVM_t* fvm = FVM_GetInstance(_el) ;
    int nn = Element_GetNbOfNodes(_el);    
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
            
            g[i*nn + j] = _dt * flow[i*nvi + j] ;
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
  Element_t* _el;
  double _t;
  double _dt;
  LocalVariables_t<V<double>> _var;
  LocalVariables_t<V<double>> _var_n;
  MPM* _mpm;
  
  
  /* Accessors */
  #if 0
  V<double>* SetInputs(double const* const* u,int const& p,V<double>& v) {
    return(_mpm->SetInputs(_el,_t,p,u,v));
  }

  template<typename T>
  V<T>* Integrate(V<double> const& v_n,V<T>& v) {
    return(_mpm->Integrate(_el,_t,_dt,v_n,v));
  }

  V<double>* Initialize(V<double>& v) {
    return(_mpm->Initialize(_el,_t,v));
  }

  int SetTangentMatrix(int const& p,V<double> const& v,V<double> const& dv,int const& k,double* c) {
    return(_mpm->SetTangentMatrix(_el,_t,_dt,p,v,dv,k,c));
  }

  int SetTransferMatrix(int const& p,V<double> const& v,double* c) {
    return(_mpm->SetTransferMatrix(_el,_dt,p,v,c));
  }

  V<double>* SetFluxes(int const& i,int const& j,V<double> const& grd,V<double>* v) {
    return(_mpm->SetFluxes(_el,_t,i,j,grd,v));
  }

  void SetIndexes(int* ind) {
    _mpm->SetIndexes(_el,ind);
  }

  void SetIncrements(double* dui) {
    _mpm->SetIncrements(_el,dui);
  }
  #endif
  
  
  void StoreImplicitTerms(int const& p) {
    _var.StoreImplicitTerms(p) ;
  }
  
  void StoreExplicitTerms(int const& p) {    
    _var.StoreExplicitTerms(_el,p) ;
  }
  
  void StoreConstantTerms(int const& p) {    
    _var.StoreConstantTerms(_el,p) ;
  }
  
  
  V<double>* DifferentiateValues(int const& p,double const& dxi,int const& i) {
    V<double>* v = _var.GetLocalValue() ;
    V<double>* dv = _var.GetLocalDerivative() ;
    double* dx = (double*) (dv + p) ;
    
    dv[p] = v[p];
  
    /* We increment the variable as (x + dx) */
    dx[i] += dxi ;
    
    {
      V<double>* v_n = _var_n.GetLocalValue() ;
      V<double>* vo = _mpm->Integrate(_el,_t,_dt,v_n[p],dv[p]) ;
      
      if(!vo) return(NULL) ;
    }
    
    dv[p] -= v[p];
    dv[p] /= dxi;
    
    return(dv + p) ;
  }
  
  /* Below:
   * V<double> = CustomValues_t<double,C...>
   * V<real> = CustomValues_t<real,C...>
   */
  #if defined HAVE_AUTODIFF
  V<double>* AutoDifferentiateValues(int const& p,int const& i) {
    V<double>* v = _var.GetLocalValue() ;
    V<real> w ;
    real* r ;
    
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
      V<double>* v_n = _var_n.GetLocalValue() ;
      V<real>* vo = _mpm->Integrate(_el,_t,_dt,v_n[p],w) ;
      
      if(!vo) return(NULL) ;
    }
  
    {
      V<double>* dv = _var.GetLocalDerivative() ;
      double* dx = (double*) (dv + p) ;
      int nbofvariables = _var.GetSizeOfLocalValues() ;
      
      for(int j = 0 ; j < nbofvariables ; j++) {
        dx[j] = r[j][1] ;
      }

      return(dv + p) ;
    }
  }
  #endif
  
  V<double>* FiniteGradient(int const i,int const j) {
    V<double>* dv = _var.GetLocalDerivative();
    
    return(_var.GradientFVM(_el,i,j,dv[i]));
  }
  
  V<double>* ComputeFluxes(int const& i,int const& j,V<double> const& grdv) {
    V<double>* v = _var.GetLocalValue() ;
    V<double>* vo = _mpm->SetFluxes(_el,_t,i,j,grdv,v);
      
    return(vo);
  }
  
  V<double>* ComputeFluxDerivatives(int const& i,int const& j,V<double> const& grdv) {
    V<double>* dv = _var.GetLocalDerivative() ;
    V<double>* dvo = _mpm->SetFluxes(_el,_t,i,j,grdv,dv);
      
    return(dvo);
  }


  int SetTangentMatrixForFEM(int const& ncols,double* c) {
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int    ind[9+Model_MaxNbOfEquations];
    double dui[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(_el,ind);
    _mpm->SetIncrements(_el,dui);
    
    {
      int dec = ncols*ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }

    for(int p = 0 ; p < np ; p++) {
      //V<double>* val = IntegrateValues(p) ;
      V<double>* val   = _var.Extract(_el,p);
      V<double>* val_n = _var_n.Extract(_el,p);
      
      if(!val) return(-1);
            
      for(int k = 0 ; k < ncols ; k++) {
        double dui_k = dui[k];
        int    ind_k = ind[k];
        
        if(ind_k >= 0) {
          V<double>* dval = DifferentiateValues(p,dui_k,ind_k) ;
        
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,p,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,p,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }


  #if defined HAVE_AUTODIFF
  int SetAutodiffTangentMatrixForFEM(int const& ncols,double* c) {
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int    ind[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(_el,ind);
    
    {
      int dec = ncols*ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }

    for(int p = 0 ; p < np ; p++) {
      //V<double>* val = IntegrateValues(p) ;
      V<double>* val   = _var.Extract(_el,p);
      V<double>* val_n = _var_n.Extract(_el,p);
      
      if(!val) return(-1);
            
      for(int k = 0 ; k < ncols ; k++) {
        int    ind_k = ind[k];
        
        if(ind_k >= 0) {
          V<double>* dval = AutoDifferentiateValues(p,ind_k) ;
        
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,p,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,p,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  #endif
  

  int SetTangentMatrixForFVM(int const& ncols,double* c) {
    int nn = Element_GetNbOfNodes(_el) ;
    int    ind[9+Model_MaxNbOfEquations];
    double dui[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(_el,ind);
    _mpm->SetIncrements(_el,dui);
    
    {
      int  dec = ncols*ncols ;
      
      for(int i = 0 ; i < nn*nn*dec ; i++) c[i] = 0 ;
    }
  
    for(int i = 0 ; i < nn ; i++) {
      //V<double>* val = IntegrateValues(i) ;
      V<double>* val   = _var.Extract(_el,i);
      V<double>* val_n = _var_n.Extract(_el,i);
    
      if(!val) return(-1) ;
    
      for(int k = 0 ; k < ncols ; k++) {
        double  dui_k = dui[k] ;
        int     ind_k = ind[k];
        
        if(ind_k >= 0) {
          V<double>* dval = DifferentiateValues(i,dui_k,ind_k) ;
      
          for(int j = 0 ; j < nn ; j++) {
            V<double>* dv = ComputeFluxDerivatives(i,j,dval[0]) ;
          }
      
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,i,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,i,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  

  #if defined HAVE_AUTODIFF
  int SetAutodiffTangentMatrixForFVM(int const& ncols,double* c) {
    int nn = Element_GetNbOfNodes(_el) ;
    int ind[9+Model_MaxNbOfEquations];
    int shift;
    
    for(int i = 0 ; i < ncols ; i++) ind[i] = -1 ;
    
    _mpm->SetIndexes(_el,ind);
    
    {
      int  dec = ncols*ncols ;
      
      for(int i = 0 ; i < nn*nn*dec ; i++) c[i] = 0 ;
    }
  
    for(int i = 0 ; i < nn ; i++) {
      //V<double>* val = IntegrateValues(i) ;
      V<double>* val   = _var.Extract(_el,i);
      V<double>* val_n = _var_n.Extract(_el,i);
    
      if(!val) return(-1) ;
    
      for(int k = 0 ; k < ncols ; k++) {
        int     ind_k = ind[k];
        
        if(ind_k >= 0) {
          V<double>* dval = AutoDifferentiateValues(i,ind_k) ;
      
          for(int j = 0 ; j < nn ; j++) {
            V<double>* dv = ComputeFluxDerivatives(i,j,dval[0]) ;
          }
      
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,i,val[0],dval[0],k,c);
        } else {
          shift = _mpm->SetTangentMatrix(_el,_t,_dt,i,val[0],val[0],k,c);
        }
      }
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
  #endif
  
  
  int SetTransferMatrixForFEM(int const& ncols,double* c) {
    IntFct_t*  intfct = Element_GetIntFct(_el) ;
    int np  = IntFct_GetNbOfPoints(intfct) ;
    int shift;
    
    /* initialization */
    {
      int dec = ncols * ncols ;
      
      for(int i = 0 ; i < np*dec ; i++) c[i] = 0 ;
    }
  
    for(int p = 0 ; p < np ; p++) {
      V<double>* val = _var.Extract(_el,p);
    
      shift = _mpm->SetTransferMatrix(_el,_dt,p,val[0],c);
      
      if(shift == 0) break;
    }

    return(shift) ;
  }
} ;


#endif
