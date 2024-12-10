#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "Plasticity with hardening (2017)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (dim)

/* Equation index */
#define E_MECH   (0)

/* Indices of unknowns (generic indices) */
#define U_MECH     E_MECH

/* Unknown index */
#define U_DISP     U_MECH



#include "BaseName.h"
#include "CustomValues.h"
#include "ConstitutiveIntegrator.h"
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



using Values_d = CustomValues_t<double,ImplicitValues_t,ExplicitValues_t,ConstantValues_t,OtherValues_t> ;

#define Values_Index(V)  CustomValues_Index(Values_d,V,double)


#define MPM_t      BaseName(_MPM_t)




struct MPM_t: public MaterialPointModel_t<Values_d> {
  MaterialPointModel_SetInputs_t<Values_d> SetInputs;
  MaterialPointModel_Integrate_t<Values_d> Integrate;
  MaterialPointModel_Initialize_t<Values_d>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_d> SetTangentMatrix;
} ;



//using CI_t = ConstitutiveIntegrator_t<Values_d,MPM_t>;
using CI_t = ConstitutiveIntegrator_t<Values_d>;




/* We define some names for implicit terms */
template<typename T = double>
struct ImplicitValues_t {
  T Displacement[3];
  T Strain[9];
  T Stress[9];
  T BodyForce[3];
  T PlasticStrain[9];
  T HardeningVariable;
  T CriterionValue;
  T PlasticMultiplier;
} ;


template<typename T = double>
struct OtherValues_t {
};



/* We define some names for explicit terms */
template<typename T = double>
struct ExplicitValues_t {
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
template<typename T = double>
struct ConstantValues_t {
  T InitialStress[9];
};

static MPM_t mpm1;
//static MPM_t* mpm = &mpm1;
static MaterialPointModel_t<Values_d>* mpm = &mpm1;

static CI_t ci(mpm) ;





/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void   GetProperties(Element_t*) ;

static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;


#define ComputeTangentStiffnessTensor(...)  Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...)                  Plasticity_ReturnMapping(plasty,__VA_ARGS__)
//#define ComputeTangentStiffnessTensor(...)  Plasticity_GenericTangentStiffnessTensor(plasty,__VA_ARGS__)
//#define ReturnMapping(...)                  Plasticity_GenericReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)              Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...)     Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)


/* Material parameters */
static double  gravity ;
static double  rho_s ;
static double* sig0 ;
static double  hardv0 ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;
static double* cijkl ;
static Plasticity_t* plasty ;
static int     plasticmodel ;
#if SharedMS_APIis(OpenMP)
  #pragma omp threadprivate(gravity,rho_s,sig0,hardv0,macrogradient,macrostrain,cijkl,plasty,plasticmodel)
#endif

#define  SetPlasticModel(I) \
         do { \
           if(plasticmodel < 0) { \
             plasticmodel = I ; \
           } else if(plasticmodel != I) { \
             Message_FatalError("Incompatible model") ; \
           } \
         } while(0)


#define GetProperty(a)      Element_GetPropertyValue(el,a)

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))



/* To treat later */
#if SharedMS_APIis(OpenMP)
  #error "Not available yet" 
  #pragma omp threadprivate(Variable,Variable_n)
#endif



int pm(const char* s)
{
         if(!strcmp(s,"gravity")) {
    return(0) ;
  } else if(!strcmp(s,"rho_s")) {
    return(1) ;
  } else if(!strcmp(s,"poisson")) {
    return(2) ;
  } else if(!strcmp(s,"young")) {
    return(3) ;
  } else if(!strcmp(s,"sig0"))  {
    return(4) ;
  } else if(!strncmp(s,"sig0_",5)) {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(4 + 3*i + j) ;
    
    /* Model 1: Drucker-Prager */
  } else if(!strcmp(s,"initial_cumulative_plastic_shear_strain")) {
    SetPlasticModel(1) ;
    return(13) ;
  } else if(!strcmp(s,"cohesion")) {
    SetPlasticModel(1) ;
    return(14) ;
  } else if(!strcmp(s,"friction")) {
    SetPlasticModel(1) ;
    return(15) ;
  } else if(!strcmp(s,"dilatancy")) {
    SetPlasticModel(1) ;
    return(16) ;
    
    /* Model 2: Cam-clay */
  } else if(!strcmp(s,"initial_pre-consolidation_pressure")) {
    SetPlasticModel(2) ;
    return(13) ;
  } else if(!strcmp(s,"slope_of_swelling_line")) {
    SetPlasticModel(2) ;
    return(14) ;
  } else if(!strcmp(s,"slope_of_virgin_consolidation_line")) {
    SetPlasticModel(2) ;
    return(15) ;
  } else if(!strcmp(s,"slope_of_critical_state_line"))  {
    SetPlasticModel(2) ;
    return(16) ;
  } else if(!strcmp(s,"initial_void_ratio")) {
    SetPlasticModel(2) ;
    return(17) ;
    
  } else if(!strcmp(s,"macro-gradient")) {
    return(18) ;
  } else if(!strncmp(s,"macro-gradient_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(18 + 3*i + j) ;
    
  } else if(!strcmp(s,"macro-fctindex")) {
    return(27) ;
  } else if(!strncmp(s,"macro-fctindex_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(27 + 3*i + j) ;
    
  } else return(-1) ;
}



double* MacroGradient(Element_t* el,double t)
{
  double* gradient = macrogradient ;
  double  f[9] = {0,0,0,0,0,0,0,0,0} ;
  
  {
    Functions_t* fcts = Material_GetFunctions(Element_GetMaterial(el)) ;
    Function_t*  fct = Functions_GetFunction(fcts) ;
    int nf = Functions_GetNbOfFunctions(fcts) ;
    double* fctindex = &GetProperty("macro-fctindex") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      int idx = floor(fctindex[i] + 0.5) ;
      
      if(0 < idx && idx < nf + 1) {
        Function_t* macrogradfct = fct + idx - 1 ;
      
        f[i] = Function_ComputeValue(macrogradfct,t) ;
      }
    }
  }
  
  {
    double* g = &GetProperty("macro-gradient") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      gradient[i] = g[i] * f[i] ;
    }
  }
  
  return(gradient) ;
}



double* MacroStrain(Element_t* el,double t)
{
  double* strain = macrostrain ;
  double* grd = MacroGradient(el,t) ;
  int i ;
    
  for(i = 0 ; i < 3 ; i++) {
    int j ;
      
    for(j = 0 ; j < 3 ; j++) {
      strain[3*i + j] = 0.5*(grd[3*i + j] + grd[3*j + i]) ;
    }
  }
  
  return(strain) ;
}



void GetProperties(Element_t* el)
{
  gravity = GetProperty("gravity") ;
  rho_s   = GetProperty("rho_s") ;
  sig0    = &GetProperty("sig0") ;
  //hardv0  = GetProperty("hardv0") ;
  
  {
    int id = SharedMS_CurrentThreadId ;
    GenericData_t* gdat = Element_FindMaterialGenericData(el,Plasticity_t,"Plasticity") ;
    int n = (gdat) ? GenericData_GetNbOfData(gdat) : 0 ;
    
    if(id < n) {
      plasty = ((Plasticity_t*) GenericData_GetData(gdat)) + id ;
      //plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") + id ;
    } else {
      arret("GetProperties") ;
    }
  
    {
      Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
      cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
      hardv0  = Plasticity_GetHardeningVariable(plasty)[0] ;
    }
  }
}



int SetModelProp(Model_t* model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  for(i = 0 ; i < dim ; i++) {
    char name_eqn[7] ;
    sprintf(name_eqn,"meca_%d",i + 1) ;
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%d",i + 1) ;
    Model_CopyNameOfUnknown(model,U_DISP + i,name_unk) ;
  }
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int NbOfProp = 36 ;
  #if SharedMS_APIis(None)
    int nthreads = 1 ;
  #else
    int nthreads = SharedMS_MaxNbOfThreads ;
  #endif
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;
  
  plasticmodel = -1 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  
  
  /* Plasticity */
  {
    plasty = Mry_Create(Plasticity_t,nthreads,Plasticity_Create()) ;

    Material_AppendData(mat,nthreads,plasty,Plasticity_t,"Plasticity") ;
  }
  
  /* Elastic and plastic properties */
  for(i = 0 ; i < nthreads ; i++) {
    Plasticity_t* plastyi = plasty + i ;
    Elasticity_t* elasty = Plasticity_GetElasticity(plastyi) ;
    
    {
      /* Elasticity */
      {
        double young   = Material_GetPropertyValue(mat,"young") ;
        double poisson = Material_GetPropertyValue(mat,"poisson") ;
        
        Elasticity_SetToIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson) ;
        
        Elasticity_UpdateElasticTensors(elasty) ;
      }
      
      /* Drucker-Prager */
      if(plasticmodel == 1) {
        double cohesion = Material_GetPropertyValue(mat,"cohesion") ;
        double af       = Material_GetPropertyValue(mat,"friction")*M_PI/180. ;
        double ad       = Material_GetPropertyValue(mat,"dilatancy")*M_PI/180. ;
        
        Plasticity_SetTo(plastyi,DruckerPrager) ;
        Plasticity_SetParameters(plastyi,af,ad,cohesion,NULL) ;
      
      /* Cam-Clay with linear elasticity */
      } else if(plasticmodel == 2) {
        double kappa  = Material_GetPropertyValue(mat,"slope_of_swelling_line") ;
        double lambda = Material_GetPropertyValue(mat,"slope_of_virgin_consolidation_line") ;
        double M      = Material_GetPropertyValue(mat,"slope_of_critical_state_line") ;
        double pc0    = Material_GetPropertyValue(mat,"initial_pre-consolidation_pressure") ;
        double e0     = Material_GetPropertyValue(mat,"initial_void_ratio") ;
        
        Plasticity_SetTo(plastyi,CamClay) ;
        Plasticity_SetParameters(plastyi,kappa,lambda,M,pc0,e0) ;
        
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
The system consists in (dim) equations\n\
\t 1. The equilibrium equation (meca_1,meca_2,meca_3)\n") ;

  printf("\n\
The primary unknowns are\n\
\t 1. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 5.8e+09   # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.3     # coefficient de Poisson\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion (Drucker-Prager model) \n") ;
  fprintf(ficd,"friction = 25     # friction angle (Drucker-Prager model) \n") ;
  fprintf(ficd,"dilatancy = 25    # dilatancy angle (Drucker-Prager model) \n") ;
  fprintf(ficd,"sig0_ij = -11.5e6 # contrainte initiale sig0_ij\n") ;
  
  return(0) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int const nvi = ((int) sizeof(ImplicitValues_t<char>));
  int const nve = ((int) sizeof(ExplicitValues_t<char>));
  int const nv0 = ((int) sizeof(ConstantValues_t<char>));

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = nvi*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = nve*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = nv0*NbOfIntPoints ;
  
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

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
  
    {
      int i ;
      
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el,double t)
{
  double* vim0 = Element_GetImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
      
  ci.Set(el,t,0,u,vim0,u,vim0) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  Element_ComputeMaterialProperties(el) ;
  
  {
    int i = ci.ComputeInitialStateByFEM();
  
    return(i);
  }
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
    
  ci.Set(el,t,dt,u_n,vim_n,u,vim0) ;

  /*
    Input data
  */
  Element_ComputeMaterialProperties(el) ;
  
  {
    int i = ci.ComputeImplicitTermsByFEM();
    
    return(i);
  }
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])

  double*  vi   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  
  int    dec = 81 ;
      
  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;
  

  /* Initialization */
  for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0 ;


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  Element_ComputeMaterialProperties(el) ;


  /*
  ** Elastoplastic matrix
  */

  {
    double* kp = ci.ComputeTangentStiffnessMatrixByFEM() ;

    for(int i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  
  return(0) ;
#undef K
}




int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double*  vi   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;

  /* Initialization */
  for(int i = 0 ; i < ndof ; i++) r[i] = 0 ;

  if(Element_IsSubmanifold(el)) return(0) ;
      
  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;

  {
    int istress = Values_Index(Stress[0]);
    int ibforce = Values_Index(BodyForce[0]);
    double* rw = ci.ComputeMechanicalEquilibiumResiduByFEM(istress,ibforce);
    
    for(int i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_MECH + j) -= rw[i*dim + j] ;
    }
  }
  
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 6 ;
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;

  //if(Element_IsSubmanifold(el)) return(0) ;

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
  Element_ComputeMaterialProperties(el) ;
  
  
  //CI_t ci(&SetInputs,&Integrate,el,t,0,u,vim0,u,vim0) ;
  
  ci.Set(el,t,0,u,vim0,u,vim0) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Displacement */
    double* pdis = Element_ComputeDisplacementVector(el,u,intfct,p,U_DISP) ;
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double crit = 0 ;
    
    for(int i = 0 ; i < dim ; i++) {
      dis[i] = pdis[i] ;
    }

    if(ItIsPeriodic) {
      for(int i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          dis[i] += MacroGradient(el,t)[3*i + j] * s[j] ;
        }
      }
    }
    
    /* Averaging */
    for(p = 0 ; p < np ; p++) {
      Values_d& val1 = *ci.ExtractValues(p);

      for(int j = 0 ; j < 9 ; j++) sig[j] += val1.Stress[j]/np ;
      
      for(int j = 0 ; j < 9 ; j++) eps_p[j] += val1.PlasticStrain[j]/np ;
      
      hardv += val1.HardeningVariable/np ;
      
      crit += val1.CriterionValue/np ;
    }
    
    {
      int i = 0 ;
      
      Result_Store(r + i++,dis   ,"Displacements",3) ;
      Result_Store(r + i++,sig   ,"Stresses",9) ;
      Result_Store(r + i++,pdis  ,"Perturbated-displacements",3) ;
      Result_Store(r + i++,eps_p ,"Plastic-strains",9) ;
      Result_Store(r + i++,&hardv,"Hardening variable",1) ;
      Result_Store(r + i++,&crit ,"Yield function",1) ;
    
      if(i != NbOfOutputs) {
        Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
      }
    }
  }
  
  return(NbOfOutputs) ;
}


int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
{
  int    dec = 81 ;
  double* c0 = c + p*dec ;
  
    /* Mechanics */
    if(k == 0) {
      double* c1 = c0 ;
      
      
      /* Tangent stiffness matrix */
      {
        /* Criterion */
        double crit = val.CriterionValue ;
        
        if(crit >= 0.) {
          double* sig = val.Stress;
          double hardv = val.HardeningVariable;
          double dlambda = val.PlasticMultiplier;
          //Plasticity_FreeBuffer(plasty) ;
          
          /* Continuum tangent stiffness matrix */
          //ComputeTangentStiffnessTensor(SIG,&HARDV) ;
          /* Consistent tangent stiffness matrix */
          ComputeTangentStiffnessTensor(sig,&hardv,&dlambda) ;
      
          CopyTangentStiffnessTensor(c1) ;
          
        } else {
      
          CopyElasticTensor(c1) ;
        }
      }
    }
  
  return(dec) ;
}



Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
{
  LocalVariables_t<Values_d> var(u,NULL);
  
  /* Displacements and strains */
  var.DisplacementVectorAndStrainFEM(el,p,U_MECH,val.Displacement) ;
  
  if(ItIsPeriodic) {
    for(int i = 0 ; i < 9 ; i++) {
      val.Strain[i]   += MacroStrain(el,t)[i] ;
    }
  }
  
  return(&val) ;
}



Values_d*  MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_d& val)
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps = val.Strain ;
  double* eps_n = val_n.Strain ;
  /* Plastic strains */
  double* eps_p  = val.PlasticStrain ;
  double const* eps_pn = val_n.PlasticStrain ;

  /* Backup stresses, plastic strains */
  {
    double* sig   = val.Stress ;
    double const* sig_n = val_n.Stress ;
    double  hardv = val_n.HardeningVariable ;
    
    
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
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      //Plasticity_FreeBuffer(plasty) ;
    
      /* Projection */
      {
        double* crit = ReturnMapping(sig,eps_p,&hardv) ;
        double* dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        val.CriterionValue  = crit[0] ;
        val.HardeningVariable = hardv ;
        val.PlasticMultiplier = dlambda[0] ;
      }
    }
  }
  
  
  /* Backup body force */
  {
    {
      double* fmass = val.BodyForce ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) fmass[i] = 0 ;
      fmass[dim - 1] = (rho_s)*gravity ;
    }
  }
  
  return(&val) ;
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  DataFile_t* datafile = Element_GetDataFile(el) ;
  
    {
      /* Initial stresses, hardening variable */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
        for(int i = 0 ; i < 9 ; i++) val.InitialStress[i] = val.Stress[i] ;
      } else {
        for(int i = 0 ; i < 9 ; i++) val.InitialStress[i] = sig0[i] ;
        for(int i = 0 ; i < 9 ; i++) val.Stress[i]  = sig0[i] ;
        val.HardeningVariable = hardv0 ;
      }
      
      /* Initial plastic strains */
      for(int i = 0 ; i < 9 ; i++) val.PlasticStrain[i]  = 0 ;
    }
  
  return(&val) ;
}
