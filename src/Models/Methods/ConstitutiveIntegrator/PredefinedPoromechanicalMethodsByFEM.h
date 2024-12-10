#ifndef PREDEFINEDPOROMECHANICALMETHODSBYFEM_H
#define PREDEFINEDPOROMECHANICALMETHODSBYFEM_H


#include "FEM.h"

#include <BaseName.h>
#include <CustomValues.h>



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

#define Values_Index(V)  CustomValues_Index(Values_t,V)


#define MPM_t      BaseName(_MPM_t)


#include "MaterialPointModel.h"
struct MPM_t: public MaterialPointModel_t<Values_d> {
  MaterialPointModel_SetInputs_t<Values_d> SetInputs;
  template<typename T>
  MaterialPointModel_Integrate_t<Values_d,Values_t<T>> Integrate;
  MaterialPointModel_Initialize_t<Values_d>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_d> SetTangentMatrix;
  MaterialPointModel_SetTransferMatrix_t<Values_d> SetTransferMatrix;
  MaterialPointModel_SetIndexes_t SetIndexes;
  MaterialPointModel_SetIncrements_t SetIncrements;
} ;



#include <ConstitutiveIntegrator.h>
using CI_t = ConstitutiveIntegrator_t<Values_d,MPM_t>;



#define ConstitutiveIntegrator_DefineElementProp(el,intfcts) \
        do {\
          IntFct_t* intfct = Element_GetIntFct(el) ;\
          int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ; \
          int const nvi = ((int) sizeof(ImplicitValues_t<char>));\
          int const nve = ((int) sizeof(ExplicitValues_t<char>));\
          int const nv0 = ((int) sizeof(ConstantValues_t<char>));\
          Element_GetNbOfImplicitTerms(el) = NbOfIntPoints*nvi ;\
          Element_GetNbOfExplicitTerms(el) = NbOfIntPoints*nve ;\
          Element_GetNbOfConstantTerms(el) = NbOfIntPoints*nv0 ;\
          return(0) ;\
        } while(0)



#define ConstitutiveIntegrator_ComputeLoads(el,t,dt,cg,r) \
        do {\
          IntFct_t* fi = Element_GetIntFct(el) ;\
          int nn = Element_GetNbOfNodes(el) ;\
          int ndof = nn*NEQ ;\
          FEM_t* fem = FEM_GetInstance(el) ;\
          double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;\
          for(int i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;\
          return(0) ;\
        } while(0)



#define ConstitutiveIntegrator_ComputeInitialState(el,t) \
        do {\
          double* vi0 = Element_GetImplicitTerm(el) ;\
          double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
          MPM_t mpm;\
          CI_t ci(&mpm) ;\
          if(Element_IsSubmanifold(el)) return(0) ;\
          ci.Set(el,t,0,u,vi0,u,vi0) ;\
          Element_ComputeMaterialProperties(el) ;\
          return(ci.ComputeInitialStateByFEM());\
        } while(0)


#define ConstitutiveIntegrator_ComputeExplicitTerms(el,t) \
        do {\
          double* vi_n = Element_GetPreviousImplicitTerm(el) ;\
          double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
          MPM_t mpm;\
          CI_t ci(&mpm) ;\
          if(Element_IsSubmanifold(el)) return(0) ;\
          ci.Set(el,t,0,u,vi_n,u,vi_n) ;\
          Element_ComputeMaterialProperties(el) ;\
          return(ci.ComputeExplicitTermsByFEM());\
        } while(0)


#define ConstitutiveIntegrator_ComputeImplicitTerms(el,t,dt) \
        do {\
          double* vi    = Element_GetCurrentImplicitTerm(el) ;\
          double* vi_n  = Element_GetPreviousImplicitTerm(el) ;\
          double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
          double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
          MPM_t mpm;\
          CI_t ci(&mpm) ;\
          if(Element_IsSubmanifold(el)) return(0) ;\
          ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
          Element_ComputeMaterialProperties(el) ;\
          return(ci.ComputeImplicitTermsByFEM()) ;\
        } while(0)


#define ConstitutiveIntegrator_ComputeMatrix(el,t,dt,k,E_MECH)\
        do {\
          double*  vi   = Element_GetCurrentImplicitTerm(el) ;\
          double*  vi_n = Element_GetPreviousImplicitTerm(el) ;\
          double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;\
          double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;\
          int nn = Element_GetNbOfNodes(el) ;\
          int ndof = nn*NEQ ;\
          MPM_t mpm;\
          CI_t ci(&mpm) ;\
          for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;\
          if(Element_IsSubmanifold(el)) return(0) ;\
          Element_ComputeMaterialProperties(el) ;\
          ci.Set(el,t,dt,u_n,vi_n,u,vi) ;\
          {\
            double* kp = ci.ComputePoromechanicalMatrixByFEM(E_MECH);\
            for(int i = 0 ; i < ndof*ndof ; i++) {\
              k[i] = kp[i] ;\
            }\
          }\
          return(0) ;\
        } while(0)




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
 * Return 0 if succeeds or < 0 if fails.
 */
{
  return(dec) ;
}

int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
/** On output:
 *  c = pointer to the transfer matrix to be filled
 * 
 *  Return the shift (i.e. the size of the matrix, i.e col*row)
 */
{
  int dec = 9 * 9 ;

  return(dec) ;
}
#endif

#endif
