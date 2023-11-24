#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "Non stationary flow surface model for claystone (2023)"
#define AUTHORS "Hao-Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (1+dim)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (29)
#define NVE     (1)
#define NV0     (0)

/* Equation index */
#define E_liq   (0)
#define E_mec   (1)

/* Unknown index */
#define U_p_l   (0)
#define U_u     (1)

/* We define some names for implicit terms */
#define M_L           (vim   + 0)[0]
#define M_L_n         (vim_n + 0)[0]
#define W_L           (vim   + 1)

#define SIG           (vim   + 4)
#define SIG_n         (vim_n + 4)

#define F_MASS        (vim   + 13)

#define EPS_P         (vim   + 16)
#define EPS_P_n       (vim_n + 16)

#define HARDV         (vim   + 25)[0]
#define HARDV_n       (vim_n + 25)[0]

#define CRIT          (vim   + 26)[0]
#define DLAMBDA       (vim   + 27)[0]

#define EPSVR_P       (vim   + 28)[0]
#define EPSVR_P_n     (vim_n + 28)[0]


/* We define some names for explicit terms */
#define K_L           (vex + 0)[0]

/* We define some names for constant terms */


/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double,double*) ;
static int    ComputeTransferCoefficients(FEM_t*,double,double*) ;

//static Model_ComputeVariables_t             ComputeVariables ;
static double* ComputeVariables(Element_t*,void*,void*,void*,const double,const double,const int);
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static void  ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;



#define ComputeTangentStiffnessTensor(...)  Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...)                  Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)              Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define UpdateElastoplasticTensor(...)      Plasticity_UpdateElastoplasticTensor(plasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...)     Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)


/* Material properties */
#define SATURATION_CURVE        (saturationcurve)
#define SaturationDegree(pc)    (Curve_ComputeValue(SATURATION_CURVE,pc))
#define dSaturationDegree(pc)   (Curve_ComputeDerivative(SATURATION_CURVE,pc))

#define RELATIVEPERM_CURVE                (relativepermcurve)
#define RelativePermeabilityToLiquid(pc)  (Curve_ComputeValue(RELATIVEPERM_CURVE,pc))





/* Parameters */
static double  gravity ;
static double  rho_s ;
static double* sig0 ;
static double  rho_l0 ;
static double  p_g = 0 ;
static double  k_int ;
static double  mu_l ;
static double  e0 ;
static double  phi0 ;
static double  hardv0 ;
static Elasticity_t* elasty ;
static Plasticity_t* plasty ;
static Curve_t* saturationcurve ;
static Curve_t* relativepermcurve ;
static double* cijkl ;



/* Variable indexes */
enum {
  I_U  = 0,
  I_U2 = I_U + 2,
  I_P_L,
  I_EPS,
  I_EPS8 = I_EPS + 8,
  I_SIG,
  I_SIG8 = I_SIG + 8,
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
  I_K_L,
  I_GRD_P_L,
  I_GRD_P_L2 = I_GRD_P_L + 2,
  I_DLAMBDA,
  I_EPSVR_P,
  I_Last
} ;

#define NbOfVariables     (I_Last)
static double  Variable[NbOfVariables] ;
static double  Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;



int pm(const char *s)
{
         if(!strcmp(s,"gravity"))    { 
    return (0) ;
  } else if(!strcmp(s,"rho_s"))      { 
    return (1) ;
  } else if(!strcmp(s,"young")) { 
    return (2) ;
  } else if(!strcmp(s,"rho_l"))      { 
    return (3) ;
  } else if(!strcmp(s,"k_int"))      { 
    return (4) ;
  } else if(!strcmp(s,"mu_l"))       { 
    return (5) ;
  } else if(!strcmp(s,"p_l0"))       {
    return (6) ;
  } else if(!strcmp(s,"initial_stress"))       {
    return(7) ;
  } else if(!strncmp(s,"initial_stress_",15))   {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(7 + 3*i + j) ;
    
    /* BBM */
  } else if(!strcmp(s,"slope_of_swelling_line")) {
    return(16) ;
  } else if(!strcmp(s,"slope_of_virgin_consolidation_line")) {
    return(17) ;
  } else if(!strcmp(s,"slope_of_critical_state_line"))  {
    return(18) ;
  } else if(!strcmp(s,"initial_pre-consolidation_pressure")) {
    return(19) ;
  } else if(!strcmp(s,"initial_porosity")) {
    return(20) ;
  } else if(!strcmp(s,"suction_cohesion_coefficient")) {
    return(22) ;
  } else if(!strcmp(s,"reference_consolidation_pressure")) {
    return(23) ;
  } else if(!strcmp(s,"poisson")) {
    return(24) ;
  } else if(!strcmp(s,"reference_strain_rate")) {
    return(25) ;
  } else if(!strcmp(s,"viscous_exponent")) {
    return(26) ;
  } else return(-1) ;
}


void GetProperties(Element_t* el)
{
  gravity = Element_GetPropertyValue(el,"gravity") ;
  rho_s   = Element_GetPropertyValue(el,"rho_s") ;
  k_int   = Element_GetPropertyValue(el,"k_int") ;
  mu_l    = Element_GetPropertyValue(el,"mu_l") ;
  rho_l0  = Element_GetPropertyValue(el,"rho_l") ;
  sig0    = &Element_GetPropertyValue(el,"initial_stress") ;
  phi0    = Element_GetPropertyValue(el,"initial_porosity") ;
  e0      = phi0/(1 - phi0) ;
  
  plasty  = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  elasty  = Plasticity_GetElasticity(plasty) ;

  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
  
  hardv0  = Plasticity_GetHardeningVariable(plasty)[0] ;
  
  saturationcurve = Element_FindCurve(el,"sl") ;
  relativepermcurve = Element_FindCurve(el,"kl") ;
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
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 27 ;
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  
  /* Plasticity */
  {
    plasty = Plasticity_Create() ;
      
    Material_AppendData(mat,1,plasty,Plasticity_t,"Plasticity") ;
  }
  
  /* Elastic and plastic properties */
  {
    elasty = Plasticity_GetElasticity(plasty) ;
    
    {
      /* Elasticity */
      {
        double young    = Material_GetPropertyValue(mat,"young") ;
        double poisson  = Material_GetPropertyValue(mat,"poisson") ;
        
        Elasticity_SetToIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson) ;
        
        Elasticity_UpdateElasticTensors(elasty) ;
      }
    }
    {
      /* NSFS model */
      {
        double kappa  = Material_GetPropertyValue(mat,"slope_of_swelling_line") ;
        double lambda = Material_GetPropertyValue(mat,"slope_of_virgin_consolidation_line") ;
        double M      = Material_GetPropertyValue(mat,"slope_of_critical_state_line") ;
        double pc0    = Material_GetPropertyValue(mat,"initial_pre-consolidation_pressure") ;
        double coh    = Material_GetPropertyValue(mat,"suction_cohesion_coefficient") ;
        double p_ref  = Material_GetPropertyValue(mat,"reference_consolidation_pressure") ;
        double refstrainrate = Material_GetPropertyValue(mat,"reference_strain_rate") ;
        double viscexp = Material_GetPropertyValue(mat,"viscous_exponent") ;
        Curve_t* lc   = Material_FindCurve(mat,"lc") ;
        Curve_t* sl   = Material_FindCurve(mat,"sl") ;
        
        phi0    = Material_GetPropertyValue(mat,"initial_porosity") ;
        e0      = phi0/(1 - phi0) ;
        
        Plasticity_SetTo(plasty,NSFS) ;
        Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0,coh,p_ref,refstrainrate,viscexp,lc,sl) ;
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
  fprintf(ficd,"p_g = 0           # gas pressure\n") ;
  fprintf(ficd,"k_int = 1e-19     # intrinsic permeability\n") ;
  fprintf(ficd,"mu_l = 0.001      # viscosity of liquid\n") ;
  fprintf(ficd,"sig0_ij = -11.5e6 # initial stress sig0_ij\n") ;
  fprintf(ficd,"Curves = my_file  # file name: p_c S_l k_rl\n") ;
  
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

      
      /* Initial stresses, hardening variables */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
      } else {
        for(i = 0 ; i < 9 ; i++) SIG[i] = sig0[i] ;
        HARDV = hardv0 ;
      }
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
      EPSVR_P = 0 ;
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
      EPSVR_P = x[I_EPSVR_P] ;
    }
    
    
    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      double rho_l = x[I_RHO_L] ;
      double p_l = x[I_P_L] ;
      double pc = p_g - p_l ;
      double k_h = rho_l*k_int/mu_l*RelativePermeabilityToLiquid(pc) ;
    
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
    
    /* pressures */
    double p_l = x[I_P_L] ;
    double pc = p_g - p_l ;
    
    /* permeability */
    double k_h = rho_l*k_int/mu_l*RelativePermeabilityToLiquid(pc) ;
    
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
      EPSVR_P = x[I_EPSVR_P] ;
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
  int dim = Element_GetDimensionOfSpace(el) ;
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
    int dec = ComputeTangentCoefficients(fem,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    /* The matrix kp is stored as (u for displacement, p for pressure)
     * | Kuu Kup |
     * | Kpu Kpp |
     * i.e. the displacements u are in the positions 0 to dim-1 and
     * the pressure p is in the position dim.
     * So we need to store the matrix by accounting for the right indexes.
     */
    #define KP(i,j)   (kp[(i)*ndof + (j)])
    
    if(E_mec == 0) {
      int i ;
      
      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    } else {
      int n ;
      
      for(n = 0 ; n < nn ; n++) {
        int m ;
        
        for(m = 0 ; m < nn ; m++) {
          
          /* Mechanics */
          {
            int i ;
      
            /* Stiffness matrix */
            for(i = 0 ; i < dim ; i++) {
              int j ;
            
              for(j = 0 ; j < dim ; j++) {
                K(E_mec + i + n*NEQ,U_u + j + m*NEQ) = KP(0 + i + n*NEQ,0 + j + m*NEQ) ;
              }
            }
          
            /* Coupling matrix */
            for(i = 0 ; i < dim ; i++) {
              K(E_mec + i + n*NEQ,U_p_l + m*NEQ) = KP(0 + i + n*NEQ,dim + m*NEQ) ;
            }
          }
          
          /* Hydraulics */
          {
            int j ;
            
            /* Coupling matrix */
            for(j = 0 ; j < dim ; j++) {
              K(E_liq + n*NEQ,U_u + j + m*NEQ) = KP(dim + n*NEQ,0 + j + m*NEQ) ;
            }
          
            /* Storage matrix */
            K(E_liq + n*NEQ,U_p_l + m*NEQ) = KP(dim + n*NEQ,dim + m*NEQ) ;
          }
        }
      }
    }
    #undef KP
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
  int dim = Element_GetDimensionOfSpace(el) ;
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
  int NbOfOutputs = 10 ;
  double* vex0  = Element_GetExplicitTerm(el) ;
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Element_GetDimensionOfSpace(el) ;
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
    double pc  = p_g - p_l ;
    /* saturation */
    double sl = SaturationDegree(pc) ;
    /* Displacement */
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps[9] = {0,0,0,0,0,0,0,0,0} ;
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double tre,e ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double k_h = 0 ;
    double strainrate = 0 ;
    double crit = 0 ;
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++) {
      double* vim  = vim0 + i*NVI ;
      double* vex  = vex0 + i*NVE ;
      double* def =  FEM_ComputeLinearStrainTensor(fem,u,intfct,i,U_u) ;
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_l[j] += W_L[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps[j] += def[j]/np ;
      
      hardv += HARDV/np ;
      
      k_h += K_L/np ;
      
      strainrate += EPSVR_P/np ;
      
      crit += CRIT/np ;
    }
    
    tre = eps[0] + eps[4] + eps[8] ;
    e   = (1 + e0) * tre ;
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Liquid_pore_pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_l      ,"Fluid_mass_flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,&sl      ,"Saturation_degree",1) ;
    Result_Store(r + i++,&e       ,"Void_ratio_variation",1) ;
    Result_Store(r + i++,eps_p    ,"Plastic_strains",9) ;
    Result_Store(r + i++,&hardv   ,"Hardening_variable",1) ;
    Result_Store(r + i++,&strainrate   ,"Plastic_strain_rate",1) ;
    Result_Store(r + i++,&crit    ,"Yield_function_value",1) ;
      
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }
  
  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(FEM_t* fem,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define T2(a,i,j)      ((a)[(i)*3+(j)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  Element_t* el  = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double*  vim0  = Element_GetCurrentImplicitTerm(el) ;
  double*  vim0_n = Element_GetPreviousImplicitTerm(el) ;
//  double*  vex0  = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 100 ;
  int    p ;
  double zero = 0. ;
  
  double dxi[Model_MaxNbOfEquations] ;
  
  {
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }


  
  for(p = 0 ; p < np ; p++) {
    double* vim  = vim0 + p*NVI ;
    double* vim_n  = vim0_n + p*NVI ;
    double* c0 = c + p*dec ;
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim0_n,t,dt,p) ;
    
    /* Pressure */
    double p_l = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    //double p_l = x[I_P_L] ;
    double pc = p_g - p_l ;
    double sl = SaturationDegree(pc) ;


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
    
      /* Effective stresses */
      {
        double pcsl = pc*sl ;
    
        sig[0] += - pcsl ;
        sig[4] += - pcsl ;
        sig[8] += - pcsl ;
      }
      
      /* Tangent stiffness matrix */
      {
        double* c1 = c0 ;
        double crit = CRIT ;

        Elasticity_CopyStiffnessTensor(elasty,c1) ;
      
        {
          /* Criterion */
          if(crit >= 0.) {
            double p_co  = HARDV ;
            //double hardv[4] = {p_co,pc,dt,EPSVR_P} ;
            double hardv[4] = {p_co,pc,dt,EPSVR_P_n} ;
            
          /* Continuum tangent stiffness matrix */
            //ComputeTangentStiffnessTensor(sig,hardv) ;
          /* Consistent tangent stiffness matrix */
            ComputeTangentStiffnessTensor(sig,hardv,&DLAMBDA) ;

            CopyTangentStiffnessTensor(c1) ;
          }
        }
      }
      
      
      /* Coupling matrix */
      {
        double  dp_l = dxi[U_p_l] ;
        double* dxdp_l = ComputeVariableDerivatives(el,t,dt,x,dp_l,I_P_L) ;
        double* dsigdp_l = dxdp_l + I_SIG ;
        double* c1 = c0 + 81 ;

        for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp_l[i] ;
      }
    }
    
    
    /* Hydraulics */
    {
      /* Fluid mass density */
      double rho_l = rho_l0 ;
    
    
      /* Coupling matrix */
      {
        double* c1 = c0 + 81 + 9 ;
        int i ;
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*sl ;
      }
      
      
      /* Storage matrix */
      {
        double* c1 = c0 + 81 + 9 + 9 ;
        //double dxk   = dxi[U_p_l] ;
        //int    k     = I_P_L ;
        //double* dx   = ComputeVariableDerivatives(el,t,dt,x,dxk,k) ;
        /* Porosity */
        double* eps  = FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
        double tre   = eps[0] + eps[4] + eps[8] ;
        double phi   = phi0 + tre ;
        
        c1[0] = - rho_l*phi*dSaturationDegree(pc) ;
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





double* ComputeVariables(Element_t* el,void* vu,void* vu_n,void* vf_n,const double t,const double dt,const int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
//  Model_t*  model  = Element_GetModel(el) ;
//  double*   x      = Model_GetVariable(model,p) ;
  /* cast when type "const void*" is used 
  const double* const* u   = (const double* const*) vu ;
  const double* const* u_n = (const double* const*) vu_n ;
  const double*        f_n = (const double*) vf_n ;
  */
  double** u   = (double**) vu ;
  double** u_n = (double**) vu_n ;
  double*  f_n = (double*)  vf_n ;
  double*  x   = Variable ;
  double*  x_n = Variable_n ;
  
    
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
      x_n[I_EPSVR_P] = EPSVR_P_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
    
    /* Pressure at previous time step */
    x_n[I_P_L] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l) ;
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;
      
      x_n[I_K_L]  = K_L ;
    }
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  /* Inputs 
   * ------*/
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x   + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Stresses */
  double* sig_n =  x_n + I_SIG ;
  /* Plastic strains */
  double* eps_pn = x_n + I_EPS_P ;
  /* Pressure */
  double  p_l   = x[I_P_L] ;
  double  p_ln  = x_n[I_P_L] ;
  double  pc   = p_g - p_l ;
  double  pc_n = p_g - p_ln ;
  double  sl   = SaturationDegree(pc) ;
  double  sl_n = SaturationDegree(pc_n) ;
    

  /* Outputs 
   * ------*/

  /* Backup stresses, plastic strains */
  {
    double* sig   = x + I_SIG ;
    double* eps_p = x + I_EPS_P ;
    
    {
      double deps[9] ;
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
      
      /* Elastic trial effective stresses */
      {
        double pcsl_n = pc_n*sl_n ;
      
        sig[0] += - pcsl_n ;
        sig[4] += - pcsl_n ;
        sig[8] += - pcsl_n ;
      }

    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Return mapping */
      {
        double p_con = x_n[I_HARDV] ; /* pre-consolidation pressure at 0 suction at the previous time step */
        double epsvr_pn = x_n[I_EPSVR_P] ;
        double hardv[4] = {p_con,pc,dt,epsvr_pn} ;
        double* crit  = ReturnMapping(sig,eps_p,hardv) ;
        double* dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        x[I_CRIT]  = crit[0] ;
        x[I_HARDV] = hardv[0] ;
        x[I_DLAMBDA] = dlambda[0] ;
        
        {
          double epsv_pn = eps_pn[0] + eps_pn[4] + eps_pn[8] ;
          double epsv_p = eps_p[0] + eps_p[4] + eps_p[8] ;
          
          x[I_EPSVR_P] = (dt > 0) ? (epsv_p - epsv_pn)/dt : 0 ;
        }
      }
    
      /* Total stresses */
      {
        double pcsl = pc*sl ;
        
        sig[0] += pcsl ;
        sig[4] += pcsl ;
        sig[8] += pcsl ;
      }
    }
  }
  
  
  /* Backup mass flow */
  {
    /* Porosity */
    double tre   = eps[0] + eps[4] + eps[8] ;
    double phi   = phi0 + tre ;
    
    /* Fluid mass density */
    double rho_l = rho_l0 ;
    
    /* Fluid mass flow */
    {
      /* Transfer coefficient */
      double k_h = x_n[I_K_L] ;
    
      /* Pressure gradient */
      double* gpl = x + I_GRD_P_L ;
    
      /* Mass flow */
      double* w_l = x + I_W_L ;
      int i ;
    
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_h*gpl[i] ;
      w_l[dim - 1] += k_h*rho_l*gravity ;
    
      /* permeability */
      x[I_K_L] = rho_l*k_int/mu_l*RelativePermeabilityToLiquid(pc) ;
    }
    
    /* Liquid mass content, body force */
    {
      double  m_l = rho_l*phi*sl ;
      double* f_mass = x + I_Fmass ;
      int i ;
    
      x[I_M_L]   = m_l ;
      x[I_RHO_L] = rho_l ;
      x[I_PHI]   = phi ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_l)*gravity ;
    }
  }
}






#if 1
double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double* dx = dVariable ;
  double* x_n = Variable_n ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,t,dt,x_n,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}
#endif
