#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "Poroplasticity with hardening (2019)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (1+dim)

/* Equation index */
#define E_MASS   (0)
#define E_MECH   (1)

/* Indices of unknowns (generic indices) */
#define U_MECH   E_MECH
#define U_MASS   E_MASS

/* Unknown index */
#define U_P_L   U_MASS
#define U_DISP  U_MECH




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




struct MPM_t: public MaterialPointModel_t<Values_t> {
  MaterialPointModel_SetInputs_t<Values_t> SetInputs;
  MaterialPointModel_Integrate_t<Values_t> Integrate;
  MaterialPointModel_Initialize_t<Values_t>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_t> SetTangentMatrix;
  MaterialPointModel_SetTransferMatrix_t<Values_t> SetTransferMatrix;
} ;





/* We define some names for implicit terms */
template<typename T = double>
struct ImplicitValues_t {
  T Displacement[3];
  T Strain[9];
  T Pressure;
  T GradPressure[3];
  T Mass_liquid;
  T MassFlow_liquid[3];
  T Stress[9];
  T BodyForce[3];
  T PlasticStrain[9];
  T HardeningVariable;
  T PlasticMultiplier;
  T YieldFunctionValue;
} ;



/* We define some names for explicit terms */
template<typename T = double>
struct ExplicitValues_t {
  T Permeability_liquid;
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
template<typename T = double>
struct ConstantValues_t {};



template<typename T = double>
struct OtherValues_t {
  T MassDensity_liquid;
  T Porosity;
};



/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*,double) ;


#define ComputeTangentStiffnessTensor(...)  Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...)             Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)         Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define UpdateElastoplasticTensor(...) Plasticity_UpdateElastoplasticTensor(plasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...)     Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)


/* Parameters */
static double  gravite ;
static double  rho_s ;
static double* sig0 ;
static double  hardv0 ;
static double  rho_l0 ;
static double  p_l0 ;
static double  phi0 ;
static double  biot,N ;
static double  k_l ;
static double  k_int,mu_l ;
static double  beta ;
static double* cijkl ;
static Plasticity_t* plasty ;
static int     plasticmodel ;

#define  SetPlasticModel(I) \
         do { \
           if(plasticmodel < 0) { \
             plasticmodel = I ; \
           } else if(plasticmodel != I) { \
             Message_FatalError("Incompatible model") ; \
           } \
         } while(0)


#define GetProperty(a)      Element_GetPropertyValue(el,a)



int pm(const char *s)
{
         if(!strcmp(s,"gravity"))    { return (0) ;
  } else if(!strcmp(s,"young"))      { return (1) ;
  } else if(!strcmp(s,"poisson"))    { return (2) ;
  } else if(!strcmp(s,"porosity"))   { return (3) ;
  } else if(!strcmp(s,"rho_l"))      { return (4) ;
  } else if(!strcmp(s,"k_int"))      { return (5) ;
  } else if(!strcmp(s,"mu_l"))       { return (6) ;
  } else if(!strcmp(s,"b"))          { return (7) ;
  } else if(!strcmp(s,"N"))          { return (8) ;
  } else if(!strcmp(s,"rho_s"))      { return (9) ;
  } else if(!strcmp(s,"beta"))       { return (10) ;
  } else if(!strcmp(s,"p_l0"))       { return (11) ;
  } else if(!strcmp(s,"k_l"))        { return (12) ;
  } else if(!strcmp(s,"sig0"))       {
    return(13) ;
  } else if(!strncmp(s,"sig0_",5))   {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(13 + 3*i + j) ;
    
  } else if(!strcmp(s,"harv0")) {
    return(22) ;
    
    /* Model 1: Drucker-Prager */
  } else if(!strcmp(s,"initial_cumulative_plastic_shear_strain")) {
    SetPlasticModel(1) ;
    return(22) ;
  } else if(!strcmp(s,"cohesion"))   { 
    SetPlasticModel(1) ;
    return (23) ;
  } else if(!strcmp(s,"friction"))   { 
    SetPlasticModel(1) ;
    return (24) ;
  } else if(!strcmp(s,"dilatancy"))  { 
    SetPlasticModel(1) ;
    return (25) ;
    
    /* Model 2: Cam-clay */
  } else if(!strcmp(s,"initial_pre-consolidation_pressure")) {
    SetPlasticModel(2) ;
    return(22) ;
  } else if(!strcmp(s,"slope_of_swelling_line")) {
    SetPlasticModel(2) ;
    return(23) ;
  } else if(!strcmp(s,"slope_of_virgin_consolidation_line")) {
    SetPlasticModel(2) ;
    return(24) ;
  } else if(!strcmp(s,"slope_of_critical_state_line"))  {
    SetPlasticModel(2) ;
    return(25) ;
  } else if(!strcmp(s,"initial_void_ratio")) {
    SetPlasticModel(2) ;
    return(26) ;
  } else if(!strcmp(s,"initial_plastic_void_ratio")) {
    SetPlasticModel(2) ;
    
  } else return(-1) ;
}


void GetProperties(Element_t* el,double t)
{
  gravite = GetProperty("gravity") ;
  phi0    = GetProperty("porosity") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l0  = GetProperty("rho_l") ;
  k_l     = GetProperty("k_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  biot    = GetProperty("b") ;
  N       = GetProperty("N") ;
  beta    = GetProperty("beta") ;
  sig0    = &GetProperty("sig0") ;
  //hardv0  = GetProperty("hardv0") ;
  
  plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
    hardv0  = Plasticity_GetHardeningVariable(plasty)[0] ;
  }
}



int SetModelProp(Model_t* model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_MASS,"liq") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_DISP + i,name_unk[i]) ;
  }
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
    
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 28 ;
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;
  
  plasticmodel = -1 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  
  /* Plasticity */
  {
    plasty = Plasticity_Create() ;
      
    Material_AppendData(mat,1,plasty,Plasticity_t,"Plasticity") ;
  }
  
  /* Elastic and plastic properties */
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    {
      /* Elasticity */
      {
        double young    = Material_GetPropertyValue(mat,"young") ;
        double poisson  = Material_GetPropertyValue(mat,"poisson") ;
        
        Elasticity_SetToIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson) ;
      
        {
          double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
          Elasticity_ComputeStiffnessTensor(elasty,c) ;
        }
      }
      
      /* Drucker-Prager */
      if(plasticmodel == 1) {
        double cohesion = Material_GetPropertyValue(mat,"cohesion") ;
        double af       = Material_GetPropertyValue(mat,"friction")*M_PI/180. ;
        double ad       = Material_GetPropertyValue(mat,"dilatancy")*M_PI/180. ;
        
        Plasticity_SetTo(plasty,DruckerPrager) ;
        Plasticity_SetParameters(plasty,af,ad,cohesion,NULL) ;
      
      /* Cam-Clay */
      } else if(plasticmodel == 2) {
        double kappa  = Material_GetPropertyValue(mat,"slope_of_swelling_line") ;
        double lambda = Material_GetPropertyValue(mat,"slope_of_virgin_consolidation_line") ;
        double M      = Material_GetPropertyValue(mat,"slope_of_critical_state_line") ;
        double pc0    = Material_GetPropertyValue(mat,"initial_pre-consolidation_pressure") ;
        double e0     = Material_GetPropertyValue(mat,"initial_void_ratio") ;
        
        Plasticity_SetTo(plasty,CamClay) ;
        Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0) ;
        
      } else {
        Message_FatalError("Unknown model") ;
      }
    }

#if 0
    {      
      Elasticity_PrintStiffnessTensor(elasty) ;
    }
#endif
  }
  
  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;

  printf("\n\
The system consists in (1 + dim) equations\n\
\t 1. The mass balance equation for water (liq)\n\
\t 2. The equilibrium equation (meca_1,meca_2,meca_3)\n") ;

  printf("\n\
The primary unknowns are\n\
\t 1. The liquid pressure (p_l)\n\
\t 2. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # mass density of solid skeleton\n") ;
  fprintf(ficd,"young = 5.8e+09   # Young's modulus\n") ;
  fprintf(ficd,"poisson = 0.3     # Poisson's ratio\n") ;
  fprintf(ficd,"porosity = 0.15   # porosity\n") ;
  fprintf(ficd,"rho_l = 1000      # mass density of fluid\n") ;
  fprintf(ficd,"p_l0 = 4.7e+06    # initial pressure of fluid\n") ;
  fprintf(ficd,"p_g = 0           # gas pressure\n") ;
  fprintf(ficd,"k_l = 2e+09       # compression modulus of fluid\n") ;
  fprintf(ficd,"k_int = 1e-19     # intrinsic permeability\n") ;
  fprintf(ficd,"mu_l = 0.001      # viscosity of liquid\n") ;
  fprintf(ficd,"b = 0.8           # Biot's coefficient\n") ;
  fprintf(ficd,"N = 4.e-11        # compressibility of pores\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion\n") ;
  fprintf(ficd,"friction = 25     # friction angle\n") ;
  fprintf(ficd,"dilatancy = 25    # dilatancy angle \n") ;
  fprintf(ficd,"beta = 0.8        # plastic Biot's coefficient\n") ;
  fprintf(ficd,"sig0_ij = -11.5e6 # initial stress sig0_ij\n") ;
  fprintf(ficd,"Curves = my_file  # file name: p_c S_l k_rl\n") ;
  
  return(0) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  
  /* Continuity of pressure across zero-thickness element */
  {
    if(Element_HasZeroThickness(el)) {
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"p_l") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"liq") ;
    }
  }
  
  {
    IntFct_t* intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;
  
    MaterialPointModel_DefineNbOfInternalValues(MPM_t,el,NbOfIntPoints);
  }
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
  
    /* hydraulic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_MASS) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
      
    /* other */
    } else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
{
  int i = MaterialPointModel_ComputeInitialStateByFEM(MPM_t,el,t);
  
  return(i);
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  int i = MaterialPointModel_ComputeExplicitTermsByFEM(MPM_t,el,t);
  
  return(i);
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  int i = MaterialPointModel_ComputeImplicitTermsByFEM(MPM_t,el,t,dt);
  
  return(i);
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
  int i = MaterialPointModel_ComputePoromechanicalMatrixByFEM(MPM_t,el,t,dt,k,E_MECH);
  
  return(i);
}




int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
  /* Initialization */
  {
    int ndof = Element_GetNbOfDOF(el) ;
    
    for(int i = 0 ; i < ndof ; i++) r[i] = 0 ;
  }
  /* 1. Mechanics */
  MaterialPointModel_ComputeMechanicalEquilibriumResiduByFEM(MPM_t,el,t,dt,r,E_MECH,Stress,BodyForce);
  /* 2. Conservation of total mass */
  MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_MASS,Mass_liquid,MassFlow_liquid);
  
  return(0);
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 9 ;
  double* vi  = Element_GetCurrentImplicitTerm(el) ;
  double* ve  = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /* Initialization */
  {
    int    i ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }
  
  /*
    Input data
  */
  GetProperties(el,t) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Pressure */
    double p_l = Element_ComputeUnknown(el,u,intfct,p,U_P_L) ;
    /* Displacement */
    double* dis = Element_ComputeDisplacementVector(el,u,intfct,p,U_DISP) ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double effsig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double crit = 0 ;
    double k_h = 0 ;
    CustomValues_t<double,ImplicitValues_t>* val1 = (CustomValues_t<double,ImplicitValues_t>*) vi ;
    CustomValues_t<double,ExplicitValues_t>* val2 = (CustomValues_t<double,ExplicitValues_t>*) ve ;
    int i;
    
    /* Averaging */
    for(i = 0 ; i < np ; i++) {
      double p_l1 = val1[i].Pressure ;
      
      for(int j = 0 ; j < 3 ; j++) w_l[j]  += val1[i].MassFlow_liquid[j]/np ;
      for(int j = 0 ; j < 9 ; j++) sig[j]  += val1[i].Stress[j]/np ;
      for(int j = 0 ; j < 9 ; j++) effsig[j] += val1[i].Stress[j]/np ;
      
      effsig[0] += beta * p_l1 / np ;
      effsig[4] += beta * p_l1 / np ;
      effsig[8] += beta * p_l1 / np ;
      
      for(int j = 0 ; j < 9 ; j++) eps_p[j] += val1[i].PlasticStrain[j]/np ;
      
      hardv += val1[i].HardeningVariable/np ;
      
      k_h += val2[i].Permeability_liquid/np ;
      
      crit += val1[i].YieldFunctionValue/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Pore pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_l      ,"Fluid mass flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,eps_p    ,"Plastic strains",9) ;
    Result_Store(r + i++,&hardv   ,"Hardening variable",1) ;
    Result_Store(r + i++,&k_h     ,"Permeability",1) ;
    Result_Store(r + i++,&crit    ,"Yield function",1) ;
    Result_Store(r + i++,effsig   ,"Effective stresses",9) ;
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }
  
  return(NbOfOutputs) ;
}




Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
{
  LocalVariables_t<Values_d> var(u,NULL);
  
  /* Displacements and strains */
  var.DisplacementVectorAndStrainFEM(el,p,U_MECH,val.Displacement) ;
    
  /* Pressure and pressure gradient */
  var.ValueAndGradientFEM(el,p,U_MASS,&val.Pressure) ;
  
  return(&val) ;
}



int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define B1(i,j)        ((c1)[(i)*3+(j)])
  
  int    dec = 100 ;
  
  if(k == 0) {
    double* c0 = c + p*dec ;
    
    /* Pressure */
    double pl  = val.Pressure ;
    
    /* Yield and potential function gradients */
    double roted = 0 ;
    
    /* Criterion */
    double crit = val.YieldFunctionValue ;


    /* initialization */
    {
      int i ;
      
      for(i = 0 ; i < dec ; i++) c0[i] = 0 ;
    }
    

    /* Mechanics */
    {
      double sig[9] ;
      int i ;
    
      for(i = 0 ; i < 9 ; i++) sig[i] = val.Stress[i] ;
    
      /* Elastic effective stresses */
      sig[0] += beta*pl ;
      sig[4] += beta*pl ;
      sig[8] += beta*pl ;
      
      /* Tangent stiffness matrix */
      {
        double* c1 = c0 ;
      
        {
          /* Criterion */
          if(crit >= 0.) {
            double hardv = val.HardeningVariable;
            double dlambda = val.PlasticMultiplier;
            /* Continuum tangent stiffness matrix */
            //ComputeTangentStiffnessTensor(sig,&HARDV) ;
            /* Consistent tangent stiffness matrix */
            ComputeTangentStiffnessTensor(sig,&hardv,&dlambda) ;
            
            CopyTangentStiffnessTensor(c1) ;
            
            {
              double fcg = Plasticity_GetFjiCijklGlk(plasty)[0] ;
              double hm  = Plasticity_GetHardeningModulus(plasty)[0] ;
              
              roted = 1/(hm + fcg) ;
            }
          
          } else {
      
            CopyElasticTensor(c1) ;
          }
        }
      }
      
      
      /* Coupling matrix */
      {
        double* c1 = c0 + 81 ;
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = - biot ;
      
        if(crit >= 0.) {
          double* dfsds = Plasticity_GetYieldFunctionGradient(plasty) ;
          double* cg = Plasticity_GetCijklGlk(plasty) ;
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          for(i = 0 ; i < 9 ; i++) {
            c1[i] -= cg[i]*(beta - biot)*trf*roted ;
          }
        }
      }
    }
    
    
    /* Hydraulics */
    {
      /* Fluid mass density */
      double rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;
    
    
      /* Coupling matrix */
      {
        double* c1 = c0 + 81 + 9 ;
        int i ;
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*biot ;
      
        if(crit >= 0.) {
          double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty) ;
          double* fc = Plasticity_GetFjiCijkl(plasty) ;
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
      
          for(i = 0 ; i < 9 ; i++) {
            c1[i] += rho_l*fc[i]*(beta - biot)*trg*roted ;
          }
        }
      }
      
      
      /* Storage matrix */
      {
        double* c1 = c0 + 81 + 9 + 9 ;
        double* eps  = val.Strain ;
        double* eps_p  = val.PlasticStrain ;
        double tre   = eps[0] + eps[4] + eps[8] ;
        double tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
        double phi_p = beta*tre_p ;
        double phi   = phi0 + biot*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
        
        c1[0] = rho_l*N + rho_l0*phi/k_l ;
      
        if(crit >= 0.) {
          double* dfsds = Plasticity_GetYieldFunctionGradient(plasty) ;
          double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty) ;
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          c1[0] += rho_l*(beta - biot)*(beta - biot)*trf*trg*roted ;
        }
      }
    }
  }
  
  return(dec) ;
#undef B1
}



int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
{
  int    dec = 9 ;

  {
    double* c1 = c + p*dec ;
    
    /* initialization */
    for(int i = 0 ; i < dec ; i++) c1[i] = 0 ;
    
    {
      /* Permeability tensor */
      c1[0] = dt*val.Permeability_liquid ;
      c1[4] = c1[0] ;
      c1[8] = c1[0] ;
    }
  }
  
  return(dec) ;
}




Values_d* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_d& val)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  val.Strain ;
  double* eps_n =  val_n.Strain ;
  /* Plastic strains */
  double* eps_p  = val.PlasticStrain ;
  double* eps_pn = val_n.PlasticStrain ;
  /* Pressure */
  double  pl   = val.Pressure ;
  double  pl_n = val_n.Pressure ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = val.Stress ;
    double* sig_n = val_n.Stress ;
    double  hardv = val_n.HardeningVariable ;
    double  dpl   = pl - pl_n ;
    
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
      
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        for(j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*deps[j] ;
        }
      }
      #undef C
      
      sig[0] += - biot*dpl ;
      sig[4] += - biot*dpl ;
      sig[8] += - biot*dpl ;
    
      /* Elastic trial effective stresses (with beta coefficient) */
      sig[0] += beta*pl ;
      sig[4] += beta*pl ;
      sig[8] += beta*pl ;
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      //Plasticity_FreeBuffer(plasty) ;
      
      /* Projection */
      {
        double* crit = ReturnMapping(sig,eps_p,&hardv) ;
        double* dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        val.YieldFunctionValue  = crit[0] ;
        val.HardeningVariable = hardv ;
        val.PlasticMultiplier = dlambda[0] ;
      }
    
      /* Total stresses */
      sig[0] -= beta*pl ;
      sig[4] -= beta*pl ;
      sig[8] -= beta*pl ;
    }
  }
  
  
  /* Backup mass flow */
  {
    /* Porosity */
    double tre   = eps[0] + eps[4] + eps[8] ;
    double tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
    double phi_p = beta*tre_p ;
    double phi   = phi0 + biot*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
    
    /* Fluid mass density */
    double rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;
    
    /* Fluid mass flow */
    {
      /* Transfer coefficient */
      double k_h = val_n.Permeability_liquid ;
    
      /* Pressure gradient */
      double* gpl = val.GradPressure ;
    
      /* Mass flow */
      double* w_l = val.MassFlow_liquid ;
      int i ;
    
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_h*gpl[i] ;
      w_l[dim - 1] += k_h*rho_l*gravite ;

      /* permeability */
      val.Permeability_liquid = rho_l*k_int/mu_l ;
    }
    
    /* Liquid mass content, body force */
    {
      double m_l = rho_l*phi ;
      double* f_mass = val.BodyForce ;
    
      val.Mass_liquid = m_l ;
      val.MassDensity_liquid = rho_l ;
      val.Porosity = phi ;
      
      for(int i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_l)*gravite ;
    }
  }
  
  return(&val);
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  DataFile_t* datafile = Element_GetDataFile(el) ;
  
  {
    {
      /* Initial stresses */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
      } else {
        for(int i = 0 ; i < 9 ; i++) val.Stress[i]  = sig0[i] ;
        val.HardeningVariable = hardv0 ;
      }
      
      for(int i = 0 ; i < 9 ; i++) val.PlasticStrain[i]  = 0 ;
    }
  }
  
  return(&val);
}
