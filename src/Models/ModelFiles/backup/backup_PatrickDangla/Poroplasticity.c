#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"

#define TITLE "Poroplasticity with hardening (2016)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (1+dim)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (27)
#define NVE     (1)
#define NV0     (0)

/* Equation index */
#define E_liq   (0)
#define E_mec   (1)

/* Unknown index */
#define U_p_l   (0)
#define U_u     (1)

/* We define some names for implicit terms */
#define M_L           (vim[0])
#define W_L           (vim + 1)
#define SIG           (vim + 4)
#define F_MASS        (vim + 13)
#define EPS_P         (vim + 16)
#define GAM_P         (vim[25])
#define CRIT          (vim[26])

#define M_L_n         (vim_n[0])
#define SIG_n         (vim_n + 4)
#define EPS_P_n       (vim_n + 16)
#define GAM_P_n       (vim_n[25])

/* We define some names for explicit terms */
#define K_L           (vex[0])

/* We define some names for constant terms */


/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double*) ;
static int    ComputeTransferCoefficients(FEM_t*,double,double*) ;

static double ReturnMapping_DruckerPrager(double*,double*,double*) ;
static double Criterion_DruckerPrager(const double*,const double,double*,double*,double*) ;
static double UpdateElastoplasticTensor(double*,double*,double*,double*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,int) ;
static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static double* ComputeVariablesDerivatives(Element_t*,double,double*,double,int) ;


/* Parameters */
static double  gravite ;
static double  rho_s ;
static double  young,poisson ;
static double  sig0_11,sig0_22,sig0_33 ;
static double  rho_l0 ;
static double  p_l0 ;
static double  phi0 ;
static double  b,N ;
static double  k_l ;
static double  k_int,mu_l ;
static double  a_int ;
static double  beta ;
static double  cohesion,af,ad ;
static double  alpha,gam_R ;


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


#define NbOfVariables     (76)
static double  Variable[NbOfVariables] ;
static double dVariable[NbOfVariables] ;

#define I_U            (0)
#define I_P_L          (3)

#define I_EPS          (4)

#define I_SIG          (13)
#define I_EPS_P        (22)
#define I_EPS_n        (31)
#define I_EPS_P_n      (40)
#define I_Fmass        (49)
#define I_M_L          (52)
#define I_W_L          (53)
#define I_GAM_P        (56)
#define I_CRIT         (57)
#define I_RHO_L        (58)
#define I_PHI          (59)
#define I_P_Ln         (60)
#define I_K_H          (61)
#define I_GRD_P_L      (62)
#define I_SIG_n        (65)
#define I_GAM_P_n      (75)


int pm(const char *s)
{
       if(!strcmp(s,"gravite"))    return (0) ;
  else if(!strcmp(s,"young"))      return (1) ;
  else if(!strcmp(s,"poisson"))    return (2) ;
  else if(!strcmp(s,"porosite"))   return (3) ;
  else if(!strcmp(s,"rho_l"))      return (4) ;
  else if(!strcmp(s,"k_int"))      return (5) ;
  else if(!strcmp(s,"mu_l"))       return (6) ;
  else if(!strcmp(s,"b"))          return (7) ;
  else if(!strcmp(s,"N"))          return (8) ;
  else if(!strcmp(s,"rho_s"))      return (9) ;
  else if(!strcmp(s,"cohesion"))   return (10) ;
  else if(!strcmp(s,"frottement")) return (11) ;
  else if(!strcmp(s,"dilatance"))  return (12) ;
  else if(!strcmp(s,"beta"))       return (13) ;
  else if(!strcmp(s,"p_l0"))       return (14) ;
  else if(!strcmp(s,"k_l"))        return (15) ;
  else if(!strcmp(s,"sig0_11"))    return (16) ;
  else if(!strcmp(s,"sig0_22"))    return (17) ;
  else if(!strcmp(s,"sig0_33"))    return (18) ;
  else if(!strcmp(s,"alpha"))      return (19) ;
  else if(!strcmp(s,"gamma_R"))    return (20) ;
  else if(!strcmp(s,"a_int"))      return (21) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  gravite = GetProperty("gravite") ;
  young   = GetProperty("young") ;
  poisson = GetProperty("poisson") ;
  phi0    = GetProperty("porosite") ;
  k_int   = GetProperty("k_int") ;
  a_int   = GetProperty("a_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l0  = GetProperty("rho_l") ;
  k_l     = GetProperty("k_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  b       = GetProperty("b") ;
  N       = GetProperty("N") ;
  beta    = GetProperty("beta") ;
  sig0_11 = GetProperty("sig0_11") ;
  sig0_22 = GetProperty("sig0_22") ;
  sig0_33 = GetProperty("sig0_33") ;
  
  cohesion       = GetProperty("cohesion") ;
  af      = GetProperty("frottement")*M_PI/180. ;
  ad      = GetProperty("dilatance")*M_PI/180. ;
  alpha   = GetProperty("alpha") ;
  gam_R   = GetProperty("gamma_R") ;
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
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 22 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
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
  

  fprintf(ficd,"gravity = 0       # gravite\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 5.8e+09   # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.3     # coefficient de Poisson\n") ;
  fprintf(ficd,"porosite = 0.15   # porosite\n") ;
  fprintf(ficd,"rho_l = 1000      # masse volumique du fluide\n") ;
  fprintf(ficd,"p_l0 = 4.7e+06    # pression initiale du fluide\n") ;
  fprintf(ficd,"p_g = 0           # pression du gaz\n") ;
  fprintf(ficd,"k_l = 2e+09       # module de compression du fluide\n") ;
  fprintf(ficd,"k_int = 1e-19     # permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001      # viscosite du liquide\n") ;
  fprintf(ficd,"b = 0.8           # coefficient de Biot\n") ;
  fprintf(ficd,"N = 4.e-11        # compressibilite des pores\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion\n") ;
  fprintf(ficd,"frottement = 25   # frottement\n") ;
  fprintf(ficd,"dilatance = 25    # dilatance \n") ;
  fprintf(ficd,"beta = 0.8        # coefficient beta\n") ;
  fprintf(ficd,"sig0_11 = -11.5e6 # contrainte initiale sig0_11\n") ;
  fprintf(ficd,"sig0_22 = -11.5e6 # contrainte initiale sig0_22\n") ;
  fprintf(ficd,"sig0_33 = -11.5e6 # contrainte initiale sig0_33\n") ;
  fprintf(ficd,"alpha = 0.5       # coefficient alpha\n") ;
  fprintf(ficd,"gamma_R = 0.015   # coefficient gamma_R\n") ;
  fprintf(ficd,"a_int = 2.e10     # coefficient a_int\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;
  
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

      //if(!Context_IsPartialInitialization) {
      if(1) {
        for(i = 0 ; i < 9 ; i++) SIG[i] = 0 ;
        SIG[0] = sig0_11 ;
        SIG[4] = sig0_22 ;
        SIG[8] = sig0_33 ;
      }
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
    
      GAM_P = 0 ;
    }
  }
  
    
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,0,p) ;
    
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
      GAM_P = x[I_GAM_P] ;
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
  double* vim = Element_GetPreviousImplicitTerm(el) ;
  double* vex = Element_GetExplicitTerm(el) ;
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
  for(p = 0 ; p < NbOfIntPoints ; p++ , vex += NVE) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim,0,p) ;
    
    /* fluid mass density */
    double rho_l = x[I_RHO_L] ;
    
    /* permeability */
    double k_h = rho_l*k_int/mu_l ;
    
    double phi = x[I_PHI] ;
    double dphi = phi - phi0 ;
    
    if(dphi > 0) {
      if(dphi < 1.e-2) {
        k_h *= 1. + a_int*dphi*dphi*dphi ;
      } else k_h *= 1. + a_int*1.e-6 ;
    }
    
    /* storage in vex */
    {
      K_L = k_h ;
    }
  }
  
  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
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
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,dt,p) ;
    
    /* storage in vim */
    {
      int    i ;
      
      M_L = x[I_M_L] ;
    
      for(i = 0 ; i < 3 ; i++) W_L[i] = x[I_W_L + i] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      GAM_P = x[I_GAM_P] ;
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
    double gam_p = 0 ;
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
      
      gam_p += GAM_P/np ;
      
      k_h += K_L/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Pore pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_l      ,"Fluid mass flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,eps_p    ,"Plastic strains",9) ;
    Result_Store(r + i++,&gam_p   ,"Cumulative plastic shear strain",1) ;
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
  double*  vim   = Element_GetCurrentImplicitTerm(el) ;
  double*  vim_n = Element_GetPreviousImplicitTerm(el) ;
  double*  vex   = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 100 ;
  int    p ;
  double zero = 0. ;
  
  double dmu,lame,mu ;
  double dxi[Model_MaxNbOfEquations] ;
  
  /*
    Elastic parameters
  */
  dmu     = young/(1+poisson) ;
  mu      = 0.5*dmu ;
  lame    = dmu*poisson/(1-2*poisson) ;
  
  
  {
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }

  
  for(p = 0 ; p < np ; p++ , vim += NVI , vex += NVE) {
    double* c0 = c + p*dec ;
    /* Variables */
    //double* x = ComputeVariables(el,u,u_n,vim_n,dt,p) ;
    
    /* Pressure */
    double pl  = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    //double pl  = x[I_P_L] ;
    
    /* Yield and potential function gradients */
    double dfsds[9] ;
    double dgsds[9] ;
    double fc[9] ;
    double cg[9] ;
    double fcg = 0 ;
    
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
        
        for(i = 0 ; i < 3 ; i++) {
          int j ;
        
          for(j = 0 ; j < 3 ; j++) {
            C1(i,i,j,j) += lame ;
            C1(i,j,i,j) += mu ;
            C1(i,j,j,i) += mu ;
          }
        }
      
        /* Criterion */
        if(crit >= 0.) {
          double hm ;
          double gam_p = GAM_P ;
          double crit_1 = Criterion_DruckerPrager(sig,gam_p,dfsds,dgsds,&hm) ;
        
          fcg = UpdateElastoplasticTensor(dfsds,dgsds,fc,cg,hm,c1) ;
          
          if(fcg < 0) return(-1) ;
        }
      }
      
      
      /* Coupling matrix */
      {
        double* c1 = c0 + 81 ;
        
        for(i = 0 ; i < 3 ; i++) B1(i,i) = - b ;
      
        if(crit >= 0.) {
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          for(i = 0 ; i < 9 ; i++) {
            c1[i] -= cg[i]*(beta - b)*trf*fcg ;
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
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
      
          for(i = 0 ; i < 9 ; i++) {
            c1[i] += rho_l*fc[i]*(beta - b)*trg*fcg ;
          }
        }
      }
      
      
      /* Storage matrix */
      {
        double* c1 = c0 + 81 + 9 + 9 ;
        //double dxk   = dxi[U_p_l] ;
        //int    k     = I_P_L ;
        //double* dx   = ComputeVariablesDerivatives(el,dt,x,dxk,k) ;
        /* Porosity */
        double* eps  = FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
        double tre   = eps[0] + eps[4] + eps[8] ;
        double tre_p = EPS_P[0] + EPS_P[4] + EPS_P[8] ;
        double phi_p = beta*tre_p ;
        double phi   = phi0 + b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
        
        c1[0] = rho_l*N + rho_l0*phi/k_l ;
      
        if(crit >= 0.) {
          double trg = dgsds[0] + dgsds[4] + dgsds[8] ;
          double trf = dfsds[0] + dfsds[4] + dfsds[8] ;
        
          c1[0] += rho_l*(beta - b)*(beta - b)*trf*trg*fcg ;
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
  double* vex = Element_GetExplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  double zero = 0. ;
  

  for(p = 0 ; p < np ; p++ , vex += NVE) {
    int i ;
    double* c1 = c + p*dec ;
    
    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */
    c1[0] = K_L ;
    c1[4] = K_L ;
    c1[8] = K_L ;
  }
  
  return(dec) ;
}



double Criterion_DruckerPrager(const double* sig,const double gam_p,double* dfsds,double* dgsds,double* hm)
/** Drucker-Prager criterion. Inputs are: 
 *  the friction angle (af), the dilatancy angle (ad) and the cohesion.
 *  Return the value of the yield function. */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9] ;
  double crit ;
  double ff,dd,cc,cc0 ;
  int    i ;
  
  /*
    Input data
  */
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc      = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*c ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  /*
    Cohesion
  */
  {
    double c1 = (gam_p < gam_R) ? 1 - (1 - alpha)*gam_p/gam_R : alpha ;
    
    cc = cc0*c1*c1 ;
  }
  
  /*
    Criterion
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q + ff*p - cc ;
  
  /*
    Deviatoric stresses
  */
  for(i = 0 ; i < 9 ; i++) {
    dev[i] = sig[i] - p*id[i] ;
  }
  
  /*
    Yield function gradient
  */
  if(q > 0.) {
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = 1.5*dev[i]/q + id[i]*ff/3. ;
    }
    
  } else {
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = id[i]*ff/3. ;
    }
  }
  
  /*
    Potential function gradient
  */
  
  /* Elastic case */
  if(crit <= 0.) {
    if(q > 0.) {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    } else {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = id[i]*dd/3. ;
      }
    }
  }
  
  /* Plastic case */
  if(crit > 0.) {
    double k   = young/(1. - 2.*poisson)/3. ;
    double dmu = young/(1.+poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    /* Flow regime at the notch apex */
    } else {
      double dl = (ff*p - cc)/(k*ff*dd) ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = dev[i]/(dmu*dl) + id[i]*dd/3. ;
      }
    }
  }
  
  /* Hardening modulus */
  *hm = 0. ;
  if(gam_p < gam_R) {
    *hm = -2.*(1.-alpha)/gam_R*(1.-(1.-alpha)*gam_p/gam_R)*cc0 ;
    *hm *= sqrt(2*j2(dgsds)) ;
  }
  
  return(crit) ;
}



double ReturnMapping_DruckerPrager(double* sig,double* eps_p,double* gam_p)
/** Drucker-Prager return mapping.
 *  Return the criterion. */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double tol = 1.e-8 ;
  double p_t,q_t,sdev_t[9] ;
  double crit ;
  double ff,dd,k,dmu,mu,cc0 ;
  int    i ;
  
  /*
    Data
  */
  dmu     = young/(1.+poisson) ;
  mu      = dmu/2. ;
  k       = young/(1. - 2.*poisson)/3. ;
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc0     = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*cohesion ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*j2(sig)) ;
  for(i = 0 ; i < 9 ; i++) {
    sdev_t[i] = sig[i] - p_t*id[i] ;
  }
  
  /*
    Criterion
  */ 
  {
    double c1 = ((*gam_p) < gam_R) ? 1 - (1 - alpha)*(*gam_p)/gam_R : alpha ;
    double cc = cc0*c1*c1 ;
    
    crit = q_t + ff*p_t - cc ;
  }
  
  /*
    Return mapping: update plastic strains and stresses
  */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double p = p_t ;
    double q = q_t ;
    double gam_pn = *gam_p ;
    double gam_p1 = *gam_p ;
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      double fcrit = crit ;
      int    nf    = 0 ;
      double dl    = 0 ;
      
      while(fabs(fcrit) > tol*cc0) {
        double dqsdl ;
        double dpsdl ;
        double dccsdl ;
        double c1 ;
        double cc ;
        
	      /* Plastic strain increments */
	      for(i = 0 ; i < 9 ; i++) {
	        deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
	      }
        
	      /* Cumulative plastic shear strain */
	      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
        
	      /* p, q, cc */
	      c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
	      cc = cc0*c1*c1 ;
	      q  = q_t - dl*3*mu ;
	      p  = p_t - dl*k*dd ;
        
	      /* dqsdl, dpsdl, dccsdl */
	      dqsdl = -3*mu ;
	      dpsdl = -k*dd ;
	      dccsdl = 0. ;
	      if(gam_p1 < gam_R) {
	        dccsdl = -2*(1 - alpha)/gam_R*c1*cc0 ;
	        dccsdl *= 1.5*sqrt(2*j2(sdev_t))/q_t ;
	      }
        
        /* Criterion */
        fcrit = q + ff*p - cc ;
        
	      /* dl */
        {
          double df = dqsdl + ff*dpsdl - dccsdl ;
          
	        dl   -= fcrit/df ;
        }
        
	      if(nf++ > 20) {
	        printf("No convergence (ReturnMapping_DruckerPrager)") ;
	        exit(0) ;
	      }
      }
      
      /* Stresses */
      for(i = 0 ; i < 9 ; i++) {
	      sig[i] = sdev_t[i]*q/q_t + p*id[i] ;
      }
    }
    
    /* Flow regime at the notch apex */
    if(q <= 0.) {
      double c1 ;
      double cc ;
      double dl ;
      
      /* Deviatoric plastic strain increments */
      for(i = 0 ; i < 9 ; i++) deps_p[i]  = sdev_t[i]/dmu ;
        
      /* Cumulative plastic shear strain */
      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
      
      /* p, q, cc */
      c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
      cc = cc0*c1*c1 ;
      p  = cc/ff ;
      q  = 0. ;
      
      /* dl */
      dl   = (ff*p_t - cc)/(k*ff*dd) ;
      
      /* Plastic strain increments and stresses */
      for(i = 0 ; i < 9 ; i++) {
	      deps_p[i] = sdev_t[i]/dmu + dl*id[i]*dd/3. ;
	      sig[i]    = p*id[i] ;
      }
    }
    
      
    /* Total plastic strains */
    for(i = 0 ; i < 9 ; i++) eps_p[i] += deps_p[i] ;
    (*gam_p) = gam_p1 ;
  }
  
  return(crit) ;
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
//  double*   x      = Model_GetVariable(model,p) ;
  double*   x      = Variable ;
  
    
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
        x[I_EPS_n   + i] = eps_n[i] ;
        x[I_SIG_n   + i] = SIG_n[i] ;
        x[I_EPS_P_n + i] = EPS_P_n[i] ;
      }
      
      x[I_GAM_P_n] = GAM_P_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
    
    /* Pressure at previous time step */
    x[I_P_Ln] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l) ;
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;
      
      x[I_K_H]  = K_L ;
    }
  }
    
  ComputeSecondaryVariables(el,dt,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x + I_EPS_n ;
  /* Plastic strains */
  double* eps_p  = x + I_EPS_P ;
  double* eps_pn = x + I_EPS_P_n ;
  /* Pressure */
  double  pl   = x[I_P_L] ;
  double  pl_n = x[I_P_Ln] ;
    


  /* Backup stresses, plastic strains */
  {
    double  dmu   = young/(1 + poisson) ;
    double  lame  = dmu*poisson/(1 - 2*poisson) ;
    double* sig   = x + I_SIG ;
    double* sig_n = x + I_SIG_n ;
    double  gam_p = x[I_GAM_P_n] ;
    double  dpl   = pl - pl_n ;
    
    
    {
      double  deps[9],trde ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
      trde = deps[0] + deps[4] + deps[8] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] + dmu*deps[i] ;
      sig[0] += lame*trde - b*dpl ;
      sig[4] += lame*trde - b*dpl ;
      sig[8] += lame*trde - b*dpl ;
    
      /* Elastic trial effective stresses (with beta coefficient) */
      sig[0] += beta*pl ;
      sig[4] += beta*pl ;
      sig[8] += beta*pl ;
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Projection */
      {
        double crit = ReturnMapping_DruckerPrager(sig,eps_p,&gam_p) ;
        
        x[I_CRIT]  = crit ;
        x[I_GAM_P] = gam_p ;
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






double* ComputeVariablesDerivatives(Element_t* el,double dt,double* x,double dxi,int i)
{
  double* dx = dVariable ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}



double UpdateElastoplasticTensor(double* dfsds,double* dgsds,double* fc,double* cg,double fcg0,double* c)
/** Update the 4th rank elastoplastic tensor in c.
 *  Inputs are: 
 *  dfsds is the yield function gradient
 *  dgsds is the potential function gradient
 *  Outputs are:
 *  fc(k,l) = dfsds(j,i) * C(i,j,k,l)
 *  cg(i,j) = C(i,j,k,l) * dgsds(l,k)
 *  fcg is updated as fcg += dfsds(j,i) * C(i,j,k,l) * dgsds(l,k))
 *  Tensor c is then updated as C(i,j,k,l) += cg(i,j) * fc(k,l) / fcg
 *  Return the inverse of fcg: 1/fcg */
{
#define C(i,j)  ((c)[(i)*9+(j)])
  double fcg = fcg0 ;
  
  /* Tangent elastoplastic tensor */
  {
      
    /* Criterion */
    {
      int i ;
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        fc[i] = 0. ;
        cg[i] = 0. ;
          
        for(j = 0 ; j < 9 ; j++) {
              
          fc[i] += dfsds[j]*C(j,i) ;
          cg[i] += C(i,j)*dgsds[j] ;
          fcg   += dfsds[i]*C(i,j)*dgsds[j] ;
        }
      }
        
      if(fcg > 0.) {
        fcg = 1./fcg ;
      } else {
            
        printf("\n") ;
        printf("dfsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dfsds[i]) ;
        }
        printf("\n") ;
        printf("dgsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dgsds[i]) ;
        }
        printf("\n") ;
        printf("fcg = %e\n",fcg) ;
        printf("\n") ;
          
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          C(i,j) -= cg[i]*fc[j]*fcg ;
        }
      }
    }
  }
  
  return(fcg) ;
#undef C
}




#if 0
double Criterion_CamClay(double *sig,double pc,double *dfsds,double *dgsds,double *hm,Element_t *el)
/* Critere de Cam-Clay */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit,m2,v ;
  int    i ;
  /*
    Donnees
  */
  kappa   = GetProperty("kappa") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  phi0    = GetProperty("phi") ;
  m2      = m*m ;
  v       = 1./(lambda - kappa) ;
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Les gradients
  */
  for(i = 0 ; i < 9 ; i++) {
    double dev = sig[i] - p*id[i] ;
    
    dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    dgsds[i] = dfsds[i] ;
  }
  
  /* Le module d'ecrouissage */
  *hm = v/(1 - phi0)*p*(2*p + pc)*pc ;
  return(crit) ;
}


double ReturnMapping_CamClay(double *sig,double *sig_n,double *p_co,double *eps_p,Element_t *el)
/* Critere de Cam-Clay : return mapping algorithm (sig,p_co,eps_p)
   (d apres Borja & Lee 1990, modifie par Dangla)
*/
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double tol = 1.e-8 ;
  double p_t,q_t,p,q,pc,crit,m2,v,a ;
  double dl ;
  int    i ;
  /*
    Donnees
  */
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  phi0    = GetProperty("phi") ;
  m2      = m*m ;
  v       = 1./(lambda - kappa) ;
  
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  pc   = *p_co ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Algorithme de projection (closest point projection)
    Une seule boucle iterative pour le calcul de p, racine de
    q*q/m2 + p*(p + pc) = 0
    Les autres variables (pc,q,dl) sont donnees explicitement par p.
  */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double dlsdp  = ((1 - phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*dlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Les variables (pc,dl,q) sont explicites en p */
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
	      printf("pas de convergence (ReturnMapping_CamClay)") ;
	      exit(0) ;
      }
    }
  }
  
  /*
    Les contraintes et deformations plastiques
  */
  a = 1./(1 + 6*mu/m2*dl) ;
  
  for(i = 0 ; i < 9 ; i++) {
    double dev      = a*(sig[i] - p_t*id[i]) ;
    double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
    sig[i]   = p*id[i] + dev ;
    eps_p[i] = dl*dfsds ;
  }
  
  /* La pression de consolidation */
  *p_co = pc ;
  return(crit) ;
}
#endif
