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
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (24)
#define NVE     (0)
#define NV0     (9)

/* Equation index */
#define E_mec   (0)

/* Unknown index */
#define U_u     (0)

/* We define some names for implicit terms */
#define SIG           (vim + 0)
#define F_MASS        (vim + 9)
#define EPS_P         (vim + 12)
#define HARDV         (vim + 21)[0]
#define CRIT          (vim + 22)[0]
#define DLAMBDA       (vim + 23)[0]

#define SIG_n         (vim_n + 0)
#define EPS_P_n       (vim_n + 12)
#define HARDV_n       (vim_n + 21)[0]

/* We define some names for explicit terms */

/* We define some names for constant terms */
#define SIG0          (v0  + 0)


/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void   GetProperties(Element_t*,double) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
//static double* ComputeVariablesDerivatives(Element_t*,double,double,double*,double,int) ;

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


/* We define some indices for the local variables */
enum {
I_U      = 0,
I_U2     = I_U + 2,
I_EPS,
I_EPS8   = I_EPS + 8,
I_SIG,
I_SIG8   = I_SIG + 8,
I_EPS_P,
I_EPS_P8 = I_EPS_P + 8,
I_FMASS,
I_FMASS2 = I_FMASS + 2,
I_HARDV,
I_CRIT,
I_DLAMBDA,
I_Last
} ;


#define NbOfVariables     (I_Last)
static double  Variable[NbOfVariables] ;
static double  Variable_n[NbOfVariables] ;
//static double dVariable[NbOfVariables] ;
#if SharedMS_APIis(OpenMP)
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



void GetProperties(Element_t* el,double t)
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
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%d",i + 1) ;
    Model_CopyNameOfUnknown(model,U_u + i,name_unk) ;
  }
  
  Model_GetComputePropertyIndex(model) = pm ;
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
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
        Plasticity_SetParameters(plastyi,af,ad,cohesion) ;
      
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

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
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
  double* v00  = Element_GetConstantTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  DataFile_t* datafile = Element_GetDataFile(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el,t) ;


  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      double* v0   = v00 + p*NV0 ;
      int    i ;
      
      /* Initial stresses, hardening variable */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = SIG[i] ;
      } else {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = sig0[i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]  = sig0[i] ;
        HARDV = hardv0 ;
      }
      
      /* Initial plastic strains */
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
    }
  }
  

  /* If there are initial displacements */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_FMASS + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
      DLAMBDA = x[I_DLAMBDA] ;
    }
  }
  
  return(0) ;
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
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el,t) ;
  
    
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_FMASS + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
      DLAMBDA = x[I_DLAMBDA] ;
    }
  }
  
  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;


  /* Initialization */
  {
    double zero = 0. ;
    int    i ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;
  }


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el,t) ;


  /*
  ** Elastoplastic matrix
  */
  {
    double c[IntFct_MaxNbOfIntPoints*81] ;
    int dec = ComputeTangentCoefficients(fem,dt,c) ;
    double* kp = FEM_ComputeElasticMatrix(fem,intfct,c,dec) ;
    
    {
      int i ;
      
      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    }
  }
  
  return(0) ;
#undef K
}




int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* vim = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialization */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;


  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  {
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  {
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_mec + dim - 1) -= -rbf[i] ;
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
  FEM_t* fem = FEM_GetInstance(el) ;

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
  GetProperties(el,t) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Displacement */
    double pdis[3] = {0,0,0} ;
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double crit = 0 ;
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      pdis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
      dis[i] = pdis[i] ;
    }

    if(ItIsPeriodic) {
      
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          dis[i] += MacroGradient(el,t)[3*i + j] * s[j] ;
        }
      }
    }
    
    /* Averaging */
    for(p = 0 ; p < np ; p++) {
      double* vim  = vim0 + p*NVI ;
      int j ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;
      
      hardv += HARDV/np ;
      
      crit += CRIT/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,dis   ,"Displacements",3) ;
    Result_Store(r + i++,sig   ,"Stresses",9) ;
    Result_Store(r + i++,pdis  ,"Perturbated-displacements",3) ;
    Result_Store(r + i++,eps_p ,"Plastic-strains",9) ;
    Result_Store(r + i++,&hardv,"Hardening variable",1) ;
    Result_Store(r + i++,&crit ,"Yield function",1) ;
  }
  
  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(FEM_t* fem,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
  Element_t* el  = FEM_GetElement(fem) ;
  double*  vim0   = Element_GetCurrentImplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 81 ;
  int    p ;
  double zero = 0. ;
  
  /*
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double dxi[Model_MaxNbOfEquations] ;
  
  {
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }
  */

  
  for(p = 0 ; p < np ; p++) {
    double* vim  = vim0 + p*NVI ;
    double* c0 = c + p*dec ;
    /* Variables */
    //double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    
    /* Initialization */
    {
      int i ;
      
      for(i = 0 ; i < dec ; i++) c0[i] = zero ;
    }
    

    /* Mechanics */
    {
      double* c1 = c0 ;
      
      
      /* Tangent stiffness matrix */
      {
        /* Criterion */
        double crit = CRIT ;
        
        if(crit >= 0.) {
          //Plasticity_FreeBuffer(plasty) ;
          
          /* Continuum tangent stiffness matrix */
          //ComputeTangentStiffnessTensor(SIG,&HARDV) ;
          /* Consistent tangent stiffness matrix */
          ComputeTangentStiffnessTensor(SIG,&HARDV,&DLAMBDA) ;
      
          CopyTangentStiffnessTensor(c1) ;
          
        } else {
      
          CopyElasticTensor(c1) ;
        }
      }
    }
  }
  
  return(dec) ;
#undef C1
#undef T4
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
//  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
//  double*   x      = Model_GetVariable(model,p) ;
  double*   x      = Variable ;
  double*   x_n    = Variable_n ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
  }
    
  /* Strains */
  {
    double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
    int    i ;
      
    if(ItIsPeriodic) {
      
      for(i = 0 ; i < 9 ; i++) {
        eps[i] += MacroStrain(el,t)[i] ;
      }
    }
    
    for(i = 0 ; i < 9 ; i++) {
      x[I_EPS + i] = eps[i] ;
    }
      
    FEM_FreeBufferFrom(fem,eps) ;
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_u) ;
      double* vim_n = f_n + p*NVI ;
      
      if(ItIsPeriodic) {
        
        for(i = 0 ; i < 9 ; i++) {
          eps_n[i] += MacroStrain(el,t-dt)[i] ;
        }
      }
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_EPS   + i] = eps_n[i] ;
        x_n[I_SIG   + i] = SIG_n[i] ;
        x_n[I_EPS_P + i] = EPS_P_n[i] ;
      }
      
      x_n[I_HARDV] = HARDV_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Plastic strains */
  double* eps_p  = x + I_EPS_P ;
  double* eps_pn = x_n + I_EPS_P ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = x + I_SIG ;
    double* sig_n = x_n + I_SIG ;
    double  hardv = x_n[I_HARDV] ;
    
    
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
        
        x[I_CRIT]  = crit[0] ;
        x[I_HARDV] = hardv ;
        x[I_DLAMBDA] = dlambda[0] ;
      }
    }
  }
  
  
  /* Backup body force */
  {
    {
      double* fmass = x + I_FMASS ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) fmass[i] = 0 ;
      fmass[dim - 1] = (rho_s)*gravity ;
    }
  }
}





#if 0
double* ComputeVariablesDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double* dx = dVariable ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,t,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}
#endif
