#ifndef MATERIALPOINTMODEL_H
#define MATERIALPOINTMODEL_H

#include "Element.h"

template<template<typename> typename V>
using MaterialPointModel_SetInputs_t = V<double>* (Element_t*,double const&,int const&,double const* const*,V<double>&);

template<template<typename> typename V,typename T = double>
using MaterialPointModel_Integrate_t = V<T>* (Element_t*,double const&,double const&,V<double> const&,V<T>&) ;

template<template<typename> typename V>
using MaterialPointModel_Initialize_t = V<double>* (Element_t*,double const&,V<double>&);

template<template<typename> typename V>
using MaterialPointModel_SetTangentMatrix_t = int (Element_t*,double const&,double const&,int const&,V<double> const&,V<double> const&,int const&,double*);

template<template<typename> typename V>
using MaterialPointModel_SetTransferMatrix_t = int (Element_t*,double const&,int const&,V<double> const&,double*);

template<template<typename> typename V>
using MaterialPointModel_SetFluxes_t = V<double>* (Element_t*,double const&,int const&,int const&,V<double> const&,V<double>*);

using MaterialPointModel_SetIndexes_t = void (Element_t*,int*);

using MaterialPointModel_SetIncrements_t = void (Element_t*,double*);



#define MaterialPointModel_TypeOfValue(MPM)    MPM::Value_type



template<template<typename> typename V>
struct MaterialPointModel_t {
  public:
  template<typename T>
  using Value_type = V<T>;
  
  //MaterialPointModel_t() {}
  virtual ~MaterialPointModel_t() {}
  
  virtual V<double>* SetInputs(Element_t*,double const&,int const&,double const* const*,V<double>&) {return(NULL);}

  virtual V<double>* Integrate(Element_t*,double const&,double const&,V<double> const&,V<double>&) {return(NULL);}

  virtual V<double>* Initialize(Element_t*,double const&,V<double>&) {return(NULL);}

  virtual int SetTangentMatrix(Element_t*,double const&,double const&,int const&,V<double> const&,V<double> const&,int const&,double*) {return(-1);}

  virtual int SetTransferMatrix(Element_t*,double const&,int const&,V<double> const&,double*) {return(-1);}

  virtual V<double>* SetFluxes(Element_t*,double const&,int const&,int const&,V<double> const&,V<double>*) {return(NULL);}

  virtual void SetIndexes(Element_t*,int*) {}

  virtual void SetIncrements(Element_t*,double*) {}
};






#if 0
/** Explanation of the parameters common to all functions.
 *  el    = pointer to element object
 *  t     = current time
 *  dt    = time increment (i.e. the previous time is t-dt)
 *  u     = pointer to pointer to primary nodal unknowns
 *  p     = p^th interpolation Gauss point of the FE
 *  val   = custom values at the current time
 *  val_n = custom values at the previous time (always an input)
 */
Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
/** On output:
 *  val is initialized with the primary nodal unknowns (strain,pressure,temperature,etc...)
 * 
 *  Return a pointer to val
 */
{
  return(&val) ;
}

Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
/** On output:
 *  val = initialize val.
 * 
 *  Return a pointer to val.
 */
{
  return(&val);
}

template <typename T>
Values_t<T>* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_t<T>& val)
/** On output:
 *  val = update from the integration of the constitutive law from t-dt to t.
 * 
 *  Return a pointer to val.
 **/
{
  return(&val) ;
}

void MPM_t::SetIndexes(Element_t* el,int* ind)
/** On ouput:
 *  ind = a pointer to an array of ncol integers
 *  ncol = nb of the tangent matrix columns
 *  
 *  ind[k] = index of the k^th unknown in the custom values struct
 */
{
}

void MPM_t::SetIncrements(Element_t* el,double* dui)
/** On ouputs:
 *  dui = a pointer to an array of ncol doubles
 *  ncol = nb of the tangent matrix columns
 *  
 *  dui[k] = arbitrary small increment of the k^th unknown used for 
 *           the numerical derivatives (see operator Differentiate)
 */
{
}

int MPM_t::SetTangentMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
/** On input:
 *  k = the k^th column of the tangent matrix to be filled.
 *  dval = the derivatives of val wrt the k^th primary unknown
 * 
 *  On output:
 *  c = pointer to the matrix to be partially filled (only the column k).
 * 
 *  Return the shift (the size of the matrix: ncols*nrows) if succeeds 
 *  or < 0 if fails.
 *  Exemples:
 *    - for an elastic matrix, shift = 81
 *    - for a poroelastic matrix, shift = 100
 */
{
  return(dec) ;
}

int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
/** On output:
 *  c = pointer to the transfer matrix to be filled
 * 
 *  Return the shift (the size of the matrix: ncols*nrows)
 *  Example: if there are ndif diffusion process, shift = 9*ndif*ndif
 */
{
  int dec = 9 * 9 ;

  return(dec) ;
}
#endif



#include "ConstitutiveIntegrator.h"

/* Macros */
#define MaterialPointModel_DefineNbOfInternalValues(MPM,el,N) \
        do {\
          int const nvi = CustomValues_NbOfImplicitValues(MaterialPointModel_TypeOfValue(MPM)<double>);\
          int const nve = CustomValues_NbOfExplicitValues(MaterialPointModel_TypeOfValue(MPM)<double>);\
          int const nv0 = CustomValues_NbOfConstantValues(MaterialPointModel_TypeOfValue(MPM)<double>);\
          Element_GetNbOfImplicitTerms(el) = N*nvi ;\
          Element_GetNbOfExplicitTerms(el) = N*nve ;\
          Element_GetNbOfConstantTerms(el) = N*nv0 ;\
        } while(0)



/* We use a C extension provided by GNU C:
 * A compound statement enclosed in parentheses may appear 
 * as an expression in GNU C.
 * (https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs) */
 
/*
ConstitutiveIntegrator_t<MaterialPointModel_TypeOfValue(std::remove_reference_t<decltype(mpm)>),std::remove_reference_t<decltype(mpm)>> MaterialPointModel_ci(&mpm);
*/

#define MaterialPointModel_InitializeValues(MPM,el,t,i) \
        ({\
          MaterialPointModel_TypeOfValue(MPM)<double> MaterialPointModel_val;\
          do {\
            double* vi = Element_GetImplicitTerm(el) ;\
            double** u = Element_ComputePointerToNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,0,u,vi,u,vi) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_val = *MaterialPointModel_ci.InitializeValues(i);\
          } while(0);\
          MaterialPointModel_val;\
        })


#define MaterialPointModel_OutputValues(MPM,el,t,i) \
        ({\
          MaterialPointModel_TypeOfValue(MPM)<double> MaterialPointModel_val;\
          do {\
            double* vi = Element_GetImplicitTerm(el) ;\
            double** u = Element_ComputePointerToNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,0,u,vi,u,vi) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_val = *MaterialPointModel_ci.IntegrateValues(i);\
          } while(0);\
          MaterialPointModel_val;\
        })

        
#define MaterialPointModel_ComputeInitialStateByFEM(MPM,el,t) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi = Element_GetImplicitTerm(el) ;\
            double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,0,u,vi,u,vi) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeInitialStateByFEM();\
          } while(0);\
          MaterialPointModel_i;\
        })
        
#define MaterialPointModel_ComputeInitialStateByFVM(MPM,el,t) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi = Element_GetImplicitTerm(el) ;\
            double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,0,u,vi,u,vi) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeInitialStateByFVM();\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeExplicitTermsByFEM(MPM,el,t) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,0,u_n,vi_n,u_n,vi_n) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeExplicitTermsByFEM();\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeExplicitTermsByFVM(MPM,el,t) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,0,u_n,vi_n,u_n,vi_n) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeExplicitTermsByFVM();\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeImplicitTermsByFEM(MPM,el,t,dt) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi    = Element_GetCurrentImplicitTerm(el);\
            double* vi_n  = Element_GetPreviousImplicitTerm(el);\
            double** u   = Element_ComputePointerToCurrentNodalUnknowns(el);\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi);\
            Element_ComputeMaterialProperties(el,t);\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeImplicitTermsByFEM();\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeImplicitTermsByFVM(MPM,el,t,dt) \
        ({\
          int MaterialPointModel_i;\
          do {\
            double* vi    = Element_GetCurrentImplicitTerm(el);\
            double* vi_n  = Element_GetPreviousImplicitTerm(el);\
            double** u   = Element_ComputePointerToCurrentNodalUnknowns(el);\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi);\
            Element_ComputeMaterialProperties(el,t);\
            MaterialPointModel_i = MaterialPointModel_ci.ComputeImplicitTermsByFVM();\
          } while(0);\
          MaterialPointModel_i;\
        })



#define MaterialPointModel_ComputeTangentStifnessMatrixByFEM(MPM,el,t,dt,k)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double*  vi   = Element_GetCurrentImplicitTerm(el) ;\
            double*  vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u    = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n  = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              double* kp = MaterialPointModel_ci.ComputeTangentStiffnessMatrixByFEM();\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = kp[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })


#ifdef USE_AUTODIFF
#define MaterialPointModel_ComputePoromechanicalMatrixByFEM(MPM,el,t,dt,k,E_MECH) \
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double*  vi   = Element_GetCurrentImplicitTerm(el) ;\
            double*  vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u    = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n  = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              double* kp = MaterialPointModel_ci.ComputeAutodiffPoromechanicalMatrixByFEM(E_MECH);\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = kp[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#else
#define MaterialPointModel_ComputePoromechanicalMatrixByFEM(MPM,el,t,dt,k,E_MECH) \
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double*  vi   = Element_GetCurrentImplicitTerm(el) ;\
            double*  vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u    = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n  = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              double* kp = MaterialPointModel_ci.ComputePoromechanicalMatrixByFEM(E_MECH);\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = kp[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#endif
        
#ifdef USE_AUTODIFF
#define MaterialPointModel_ComputeMassConservationMatrixByFEM(MPM,el,t,dt,k)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u   = Element_ComputePointerToNodalUnknowns(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            {\
              double* km = MaterialPointModel_ci.ComputeAutodiffMassConservationMatrixByFEM();\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#else
#define MaterialPointModel_ComputeMassConservationMatrixByFEM(MPM,el,t,dt,k)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u   = Element_ComputePointerToNodalUnknowns(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            {\
              double* km = MaterialPointModel_ci.ComputeMassConservationMatrixByFEM();\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#endif
         
#ifdef USE_AUTODIFF
#define MaterialPointModel_ComputeMassConservationMatrixByFVM(MPM,el,t,dt,k)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u   = Element_ComputePointerToNodalUnknowns(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            {\
              double* km = MaterialPointModel_ci.ComputeAutodiffMassConservationMatrixByFVM();\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#else
#define MaterialPointModel_ComputeMassConservationMatrixByFVM(MPM,el,t,dt,k)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u   = Element_ComputePointerToNodalUnknowns(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            int ndof = Element_GetNbOfDOF(el);\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0;\
            if(Element_IsSubmanifold(el)) return(0) ;\
            Element_ComputeMaterialProperties(el,t) ;\
            {\
              double* km = MaterialPointModel_ci.ComputeMassConservationMatrixByFVM();\
              for(int i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })
#endif

#define MaterialPointModel_ComputeMechanicalEquilibriumResiduByFEM(MPM,el,t,dt,r,E_MECH,Stress,BodyForce)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              int istress = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,Stress[0],double);\
              int ibforce = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,BodyForce[0],double);\
              double* rw = MaterialPointModel_ci.ComputeMechanicalEquilibiumResiduByFEM(istress,ibforce);\
              int nn = Element_GetNbOfNodes(el) ;\
              int dim = Element_GetDimensionOfSpace(el) ;\
              int neq = Element_GetNbOfEquations(el);\
              for(int i = 0 ; i < nn ; i++) {\
                for(int j = 0 ; j < dim ; j++) {\
                  r[(i)*neq + (E_MECH + j)] -= rw[i*dim + j] ;\
                }\
              }\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeMassConservationResiduByFEM(MPM,el,t,dt,r,E_MASS,Mass,MassFlow)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              int imass = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,Mass,double);\
              int iflow = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,MassFlow[0],double);\
              double* ra =  MaterialPointModel_ci.ComputeMassConservationResiduByFEM(imass,iflow);\
              int nn = Element_GetNbOfNodes(el) ;\
              int neq = Element_GetNbOfEquations(el);\
              for(int i = 0 ; i < nn ; i++) {\
                r[(i)*neq + (E_MASS)] -= ra[i] ;\
              }\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeMassConservationResiduByFVM(MPM,el,t,dt,r,E_MASS,Mass,MassFlow)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              int imass = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,Mass,double);\
              int iflow = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,MassFlow[0],double);\
              double* ra =  MaterialPointModel_ci.ComputeMassConservationResiduByFVM(imass,iflow);\
              int nn = Element_GetNbOfNodes(el) ;\
              int neq = Element_GetNbOfEquations(el);\
              for(int i = 0 ; i < nn ; i++) {\
                r[(i)*neq + (E_MASS)] -= ra[i] ;\
              }\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeFluxResiduByFVM(MPM,el,t,dt,r,E_MASS,MassFlow)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              int iflow = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,MassFlow[0],double);\
              double* ra =  MaterialPointModel_ci.ComputeFluxResiduByFVM(iflow);\
              int nn = Element_GetNbOfNodes(el) ;\
              int neq = Element_GetNbOfEquations(el);\
              for(int i = 0 ; i < nn ; i++) {\
                r[(i)*neq + (E_MASS)] -= ra[i] ;\
              }\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })


#define MaterialPointModel_ComputeBodyForceResiduByFVM(MPM,el,t,dt,r,E_MASS,Mass)\
        ({\
          int MaterialPointModel_i = 0;\
          do {\
            double* vi   = Element_GetCurrentImplicitTerm(el) ;\
            double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
            double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
            double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
            MPM mpm;\
            ConstitutiveIntegrator_t<MPM> MaterialPointModel_ci(&mpm);\
            if(Element_IsSubmanifold(el)) return(0) ;\
            MaterialPointModel_ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
            {\
              int imass = CustomValues_Index(MaterialPointModel_TypeOfValue(MPM)<double>,Mass,double);\
              double* ra =  MaterialPointModel_ci.ComputeBodyForceResiduByFVM(imass);\
              int nn = Element_GetNbOfNodes(el) ;\
              int neq = Element_GetNbOfEquations(el);\
              for(int i = 0 ; i < nn ; i++) {\
                r[(i)*neq + (E_MASS)] -= ra[i] ;\
              }\
            }\
          } while(0);\
          MaterialPointModel_i;\
        })

#endif
