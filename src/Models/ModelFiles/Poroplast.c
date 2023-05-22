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
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (28)
#define NVE     (1)
#define NV0     (0)

/* Equation index */
#define E_liq   (0)
#define E_mec   (1)

/* Unknown index */
#define U_p_l   (0)
#define U_u     (1)

/* We define some names for implicit terms */
#define M_L           (vim)[0]
#define W_L           (vim + 1)
#define SIG           (vim + 4)
#define F_MASS        (vim + 13)
#define EPS_P         (vim + 16)
#define HARDV         (vim + 25)[0]
#define CRIT          (vim + 26)[0]
#define DLAMBDA       (vim + 27)[0]

#define M_L_n         (vim_n)[0]
#define SIG_n         (vim_n + 4)
#define EPS_P_n       (vim_n + 16)
#define HARDV_n       (vim_n + 25)[0]

/* We define some names for explicit terms */
#define K_L           (vex)[0]

/* We define some names for constant terms */


/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double*) ;
static int    ComputeTransferCoefficients(FEM_t*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static int     ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
//static double* ComputeVariablesDerivatives(Element_t*,double,double,double*,double,int) ;


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
static double  b,N ;
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


/* We define some indices for the local variables */
enum {
I_U      = 0,
I_U2     = I_U + 2,
I_P_L,
I_EPS,
I_EPS8   = I_EPS + 8,
I_SIG,
I_SIG8   = I_SIG + 8,
I_EPS_P,
I_EPS_P8 = I_EPS_P + 8,
I_Fmass,
I_Fmass2 = I_Fmass + 2,
I_M_L,
I_W_L,
I_W_L2 = I_W_L + 2,
I_HARDV,
I_CRIT,
I_RHO_L,
I_PHI,
I_K_H,
I_GRD_P_L,
I_GRD_P_L2 = I_GRD_P_L + 2,
I_DLAMBDA,
I_Last
} ;

#define NbOfVariables     (I_Last)
static double  Variable[NbOfVariables] ;
static double  Variable_n[NbOfVariables] ;
//static double dVariable[NbOfVariables] ;



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


void GetProperties(Element_t* el)
{
  gravite = GetProperty("gravity") ;
  phi0    = GetProperty("porosity") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l0  = GetProperty("rho_l") ;
  k_l     = GetProperty("k_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  b       = GetProperty("b") ;
  N       = GetProperty("N") ;
  beta    = GetProperty("beta") ;
  sig0    = &GetProperty("sig0") ;
  hardv0  = GetProperty("hardv0") ;
  
  plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
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
  Model_CopyNameOfEquation(model,E_liq,"liq") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_p_l,"p_l") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_u + i,name_unk[i]) ;
  }
  
  Model_GetComputePropertyIndex(model) = pm ;
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
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
        Plasticity_SetParameters(plasty,af,ad,cohesion) ;
      
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
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

    /** Define the length of tables */
    Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
    Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_liq) {
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
  double* vim0  = Element_GetImplicitTerm(el) ;
  double* vex0  = Element_GetExplicitTerm(el) ;
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
  GetProperties(el) ;
  
  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;

      
      /* Initial stresses */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
      } else {
        for(i = 0 ; i < 9 ; i++) SIG[i]  = sig0[i] ;
        HARDV = hardv0 ;
      }
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
    }
  }
  
    
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,0,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      M_L = x[I_M_L] ;
    
      for(i = 0 ; i < 3 ; i++) W_L[i] = x[I_W_L + i] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
      DLAMBDA = x[I_DLAMBDA] ;
    }
    
    
    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      double rho_l = x[I_RHO_L] ;
      double k_h = rho_l*k_int/mu_l ;
    
      K_L = k_h ;
    }
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  double* vim0 = Element_GetPreviousImplicitTerm(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    
    /* fluid mass density */
    double rho_l = x[I_RHO_L] ;
    
    /* permeability */
    double k_h = rho_l*k_int/mu_l ;
    
    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      
      K_L = k_h ;
    }
  }
  
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
  GetProperties(el) ;
  
    
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      M_L = x[I_M_L] ;
    
      for(i = 0 ; i < 3 ; i++) W_L[i] = x[I_W_L + i] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
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
  GetProperties(el) ;


  /*
  ** Poromechanics matrix
  */
  {
    double c[IntFct_MaxNbOfIntPoints*100] ;
    int dec = ComputeTangentCoefficients(fem,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1,U_u) ;
    
    {
      int i ;
      
      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    }
  }
  
  /*
  ** Conduction Matrix
  */
  {
    double c[IntFct_MaxNbOfIntPoints*100] ;
    int dec = ComputeTransferCoefficients(fem,dt,c) ;
    double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    int    i ;
  
    for(i = 0 ; i < nn ; i++) {
      int    j ;
      
      for(j = 0 ; j < nn ; j++) {
        K(E_liq + i*NEQ,U_p_l + j*NEQ) += dt*kc[i*nn + j] ;
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
  double* vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
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
    double* vim = vim_1 ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  {
    double* vim = vim_1 ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_mec + dim - 1) -= -rbf[i] ;
    }
    
  }
  
  
  /* 2. Hydraulics */
  
  /* 2.1 Accumulation Terms */
  {
    double* vim = vim_1 ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_L - M_L_n ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_liq) -= ra[i] ;
    }
  }
  
  /* 2.2 Transport Terms */
  {
    double* vim = vim_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_L,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= -dt*rf[i] ;
  }
  
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 7 ;
  double* vex  = Element_GetExplicitTerm(el) ;
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t* fem = FEM_GetInstance(el) ;

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
  GetProperties(el) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Pressure */
    double p_l = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    /* Displacement */
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double k_h = 0 ;
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_l[j] += W_L[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;
      
      hardv += HARDV/np ;
      
      k_h += K_L/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Pore pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_l      ,"Fluid mass flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,eps_p    ,"Plastic strains",9) ;
    Result_Store(r + i++,&hardv   ,"Hardening variable",1) ;
    Result_Store(r + i++,&k_h     ,"Permeability",1) ;
  }
  
  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(FEM_t* fem,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define T2(a,i,j)      ((a)[(i)*3+(j)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  Element_t* el  = FEM_GetElement(fem) ;
  double*  vim0  = Element_GetCurrentImplicitTerm(el) ;
//  double*  vim_n = Element_GetPreviousImplicitTerm(el) ;
//  double*  vex0  = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
//  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 100 ;
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
    
    /* Pressure */
    double pl  = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    //double pl  = x[I_P_L] ;
    
    /* Yield and potential function gradients */
    double roted = 0 ;
    
    /* Criterion */
    double crit = CRIT ;


    /* initialization */
    {
      int i ;
      
      for(i = 0 ; i < dec ; i++) c0[i] = zero ;
    }
    

    /* Mechanics */
    {
      double sig[9] ;
      int i ;
    
      for(i = 0 ; i < 9 ; i++) sig[i] = SIG[i] ;
    
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
            /* Continuum tangent stiffness matrix */
            //ComputeTangentStiffnessTensor(sig,&HARDV) ;
            /* Consistent tangent stiffness matrix */
            ComputeTangentStiffnessTensor(sig,&HARDV,DLAMBDA) ;
            
            CopyTangentStiffnessTensor(c1) ;
            
            {
              double fcg = Plasticity_GetFjiCijklGlk(plasty) ;
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
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = - b ;
      
        if(crit >= 0.) {
          double* dfsds = Plasticity_GetYieldFunctionGradient(plasty) ;
          double* cg = Plasticity_GetCijklGlk(plasty) ;
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          for(i = 0 ; i < 9 ; i++) {
            c1[i] -= cg[i]*(beta - b)*trf*roted ;
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
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*b ;
      
        if(crit >= 0.) {
          double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty) ;
          double* fc = Plasticity_GetFjiCijkl(plasty) ;
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
      
          for(i = 0 ; i < 9 ; i++) {
            c1[i] += rho_l*fc[i]*(beta - b)*trg*roted ;
          }
        }
      }
      
      
      /* Storage matrix */
      {
        double* c1 = c0 + 81 + 9 + 9 ;
        //double dxk   = dxi[U_p_l] ;
        //int    k     = I_P_L ;
        //double* dx   = ComputeVariablesDerivatives(el,t,dt,x,dxk,k) ;
        /* Porosity */
        double* eps  = FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
        double tre   = eps[0] + eps[4] + eps[8] ;
        double tre_p = EPS_P[0] + EPS_P[4] + EPS_P[8] ;
        double phi_p = beta*tre_p ;
        double phi   = phi0 + b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
        
        c1[0] = rho_l*N + rho_l0*phi/k_l ;
      
        if(crit >= 0.) {
          double* dfsds = Plasticity_GetYieldFunctionGradient(plasty) ;
          double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty) ;
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          c1[0] += rho_l*(beta - b)*(beta - b)*trf*trg*roted ;
        }
        //c1[0] = dx[I_M_L] ;
      }
    }
  }
  
  return(dec) ;
#undef C1
#undef B1
#undef T2
#undef T4
}



int ComputeTransferCoefficients(FEM_t* fem,double dt,double* c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  Element_t* el = FEM_GetElement(fem) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  double zero = 0. ;
  

  for(p = 0 ; p < np ; p++) {
    int i ;
    double* c1 = c + p*dec ;
    
    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    {
      double* vex  = vex0 + p*NVE ;
      
      /* Permeability tensor */
      c1[0] = K_L ;
      c1[4] = K_L ;
      c1[8] = K_L ;
    }
  }
  
  return(dec) ;
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
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Pressure */
    x[I_P_L] = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    
    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_p_l) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P_L + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_u) ;
      double* vim_n = f_n + p*NVI ;
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_EPS   + i] = eps_n[i] ;
        x_n[I_SIG   + i] = SIG_n[i] ;
        x_n[I_EPS_P + i] = EPS_P_n[i] ;
      }
      
      x_n[I_HARDV] = HARDV_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
    
    /* Pressure at previous time step */
    x_n[I_P_L] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l) ;
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;
      
      x[I_K_H]  = K_L ;
    }
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



int  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Plastic strains */
  double* eps_p  = x + I_EPS_P ;
  double* eps_pn = x_n + I_EPS_P ;
  /* Pressure */
  double  pl   = x[I_P_L] ;
  double  pl_n = x_n[I_P_L] ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = x + I_SIG ;
    double* sig_n = x_n + I_SIG ;
    double  hardv = x_n[I_HARDV] ;
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
      
      sig[0] += - b*dpl ;
      sig[4] += - b*dpl ;
      sig[8] += - b*dpl ;
    
      /* Elastic trial effective stresses (with beta coefficient) */
      sig[0] += beta*pl ;
      sig[4] += beta*pl ;
      sig[8] += beta*pl ;
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Projection */
      {
        double crit = ReturnMapping(sig,eps_p,&hardv) ;
        double dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        x[I_CRIT]  = crit ;
        x[I_HARDV] = hardv ;
        x[I_DLAMBDA] = dlambda ;
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
    double phi   = phi0 + b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
    
    /* Fluid mass density */
    double rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;
    
    /* Fluid mass flow */
    {
      /* Transfer coefficient */
      double k_h = x[I_K_H] ;
    
      /* Pressure gradient */
      double* gpl = x + I_GRD_P_L ;
    
      /* Mass flow */
      double* w_l = x + I_W_L ;
      int i ;
    
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_h*gpl[i] ;
      w_l[dim - 1] += k_h*rho_l*gravite ;
    }
    
    /* Liquid mass content, body force */
    {
      double m_l = rho_l*phi ;
      double* f_mass = x + I_Fmass ;
      int i ;
    
      x[I_M_L] = m_l ;
      x[I_RHO_L] = rho_l ;
      x[I_PHI] = phi ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_l)*gravite ;
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
