#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

/* The Finite Volume Method */
#include "FVM.h"


#define TEMPERATURE   (293)

#define TITLE   "External sulfate attack of concrete (2016)" 
#define AUTHORS "Gu-Dangla"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ     (1)

/* Nb of (im/ex)plicit and constant terms */
#define NVE     (1)
#define NVI     (13)
#define NV0     (2)


/* Equation Indexes */
#define E_S     (0)

/* Primary Unknown Indexes */
#define U_C_SO4 (0)


/* Compiling options */
#define NOLOG_U     1
#define LOG_U       2
#define Ln10        (2.302585093)
#define U_SO4       LOG_U


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWN_n(n,i)   Element_GetValueOfNodalUnknown(el,u_n,n,i)


#if (U_SO4 == LOG_U)
  #define LogC_SO4(n)   (UNKNOWN(n,U_C_SO4))
  #define LogC_SO4n(n)  (UNKNOWN_n(n,U_C_SO4))
  #define C_SO4(n)      (pow(10,UNKNOWN(n,U_C_SO4)))
  #define C_SO4n(n)     (pow(10,UNKNOWN_n(n,U_C_SO4)))
#else
  #define C_SO4(n)      (UNKNOWN(n,U_C_SO4))
  #define C_SO4n(n)     (UNKNOWN_n(n,U_C_SO4))
  #define LogC_SO4(n)   (log10(UNKNOWN(n,U_C_SO4)))
  #define LogC_SO4n(n)  (log10(UNKNOWN_n(n,U_C_SO4)))
#endif


#define N_S(n)     (f[(n)])
#define W_S        (f[2])
#define N_CH(n)    (f[(3+(n))])
#define PoreRadius(n)    (f[(5+(n))])
#define PHI(n)     (f[(7+(n))])
#define N_AFt(n)   (f[(9+(n))])
#define N_C3AH6(n) (f[(11+(n))])

#define N_Sn(n)     (f_n[(n)])
#define PoreRadiusn(n)    (f_n[(5+(n))])
#define PHIn(n)     (f_n[(7+(n))])
#define N_AFtn(n)   (f_n[(9+(n))])
#define N_C3AH6n(n) (f_n[(11+(n))])


#define KF_SO4     (va[(0)])


#define V_Cem0(n)  (v0[(0+n)])


/* Water properties */
/* Molar mass (kg/mol) */
#define M_H2O     (18.e-3)
/* Molar volume of liquid water (m3/mol) */
#define V_H2O     (18.e-6)
/* Autoprotolysis constant (=1e-14*1.e6) */
#define LogK_H2O     (-8)



/* Material properties */
#define SATURATION_CURVE           (Element_GetCurve(el) + 0)
#define LiquidSaturationDegree(r)  (Curve_ComputeValue(SATURATION_CURVE,r))
#define dLiquidSaturationDegree(r) (Curve_ComputeDerivative(SATURATION_CURVE,r))
#define PoreEntryRadiusMax         (Curve_GetXRange(SATURATION_CURVE)[1])



/* CH Properties */
/* Molar volume of CH solid (m3/mole) */
#define V_CH       (33.e-6)      /* (33.e-3) */
/* Equilibrium constant (logK = -5.14 + (3*3)) */
#define LogK_CH        (3.86)


/* C3AH6 Properties */
/* Molar volume of C3AH6 solid (m3/mole) */
#define V_C3AH6       (149.52e-6)
/* Equilibrium constant (logK = -89.75 + (3*17)) */
#define LogK_C3AH6       (-38.75)


/* AFt Properties ((Ca6).(Al2).3(SO4).12(OH).26(H2O)) */
/* Molar volume of AFt solid (m3/mole) */
#define V_AFt      (710.32e-6)      /* Thermochimie (ANDRA) */
/* Equilibrium constant (logK = -112 + (3*23)) */
#define LogK_AFt       (-43)
/* Surface tension (N/m) */
#define Gamma_AFt   (0.1)
/* Equilibrium saturation index of AFt */
#define EquilibriumAFtSaturationIndex(r)  (exp(2*Gamma_AFt*V_AFt/(RT*r)))
#define dEquilibriumAFtSaturationIndex(r) (-2*Gamma_AFt*V_AFt/(RT*r*r)*EquilibriumAFtSaturationIndex(r))
#define InverseOfEquilibriumAFtSaturationIndex(b)  (2*Gamma_AFt*V_AFt/(RT*log(b)))
#define AFtSolidContent(n,s,dt)     MAX((n + dt*r_aft*(s - 1)),0.)



/* AFm Properties ((Ca4).(Al2).(SO4).12(OH).6(H2O)) */
/* Equilibrium constant (logK = -96 + (3*19)) */
#define LogK_AFm       (-39)


/* CSH2 Properties */
/* Molar volume of CSH2 crystal (m3/mole) */
#define V_CSH2     (75.e-6)      /* (75.e-3) */
/* Equilibrium constant (logK = -4.58 + (2*3)) */
#define LogK_CSH2       (1.42)


/* Physical constants (J/mol) */
#define RT     (2436.)



/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Fonctions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;



static double* ComputeVariables(Element_t*,double**,double**,double*,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double*,double,int) ;


static void    ComputeTransferCoefficients(Element_t*) ;
static double* ComputeVariableFluxes(Element_t*,int,int) ;

static double* ComputeFluxes(Element_t*,double*,int,int) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;


static double  Radius(double,double,double,Element_t*) ;




/* Parameters */
static double phi0 ;
static double n_ch0 ;
static double n_c3ah60 ;
static double ph ;
static double a_AFt ;
static double d_so4 ;



#include "DiffusionCoefficientOfMoleculeInWater.h"

void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (m2/s) */
  d_so4   = DiffusionCoefficientOfMoleculeInWater(SO4,TK) ;
}



#define NN                (2)
#define NbOfVariables    (NEQ+16)
static double Variables[NN][NbOfVariables] ;
static double dVariables[2*NbOfVariables] ;

#define I_C_SO4        (NEQ+0)

#define I_N_S          (NEQ+1)

#define I_N_CH         (NEQ+2)
#define I_N_AFt        (NEQ+3)
#define I_N_C3AH6      (NEQ+4)

#define I_Radius       (NEQ+8)
#define I_Radiusn      (NEQ+9)

#define I_PHI          (NEQ+10)
#define I_PHIn         (NEQ+11)

#define I_S_AFt        (NEQ+12)

#define I_S_C          (NEQ+13)

#define I_V_Cem0       (NEQ+15)

  
  

#define NbOfVariableFluxes    (1)
static double VariableFluxes[NbOfVariableFluxes] ;

#define I_W_S           (0)



int pm(const char *s)
{
  if(strcmp(s,"porosity") == 0)        return (0) ;
  else if(strcmp(s,"N_CH") == 0)       return (1) ;
  else if(strcmp(s,"N_C3AH6") == 0)    return (2) ;
  else if(strcmp(s,"pH") == 0)         return (3) ;
  else if(strcmp(s,"a_AFt") == 0)      return (4) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  phi0      = GetProperty("porosity") ;
  n_ch0     = GetProperty("N_CH") ;
  n_c3ah60  = GetProperty("N_C3AH6") ;
  ph        = GetProperty("pH") ;
  a_AFt     = GetProperty("a_AFt") ;
}




int SetModelProp(Model_t* model)
/** Set the model properties 
 *  Return 0 */
{
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_S,"sulfur") ;
  
  /** Names of the main (nodal) unknowns */
#if (U_SO4 == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_SO4,"logc_so4") ;
#else
  Model_CopyNameOfUnknown(model,U_C_SO4,"c_so4") ;
#endif

  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
  }
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 4 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE* ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- First equation      (sulfur)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- First unknown       (c_so4)\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"prop1 = 0.01   # Property 1\n") ;
  fprintf(ficd,"prop2 = 1.e-3  # Property 2\n") ;
  fprintf(ficd,"prop3 = 1.e6   # Property 3\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}




int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* f = Element_GetImplicitTerm(el) ;
  double* v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int i ;
  
  /* We skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  for(i = 0 ; i < nn ; i++) {
    double v_cem      = V_CH*n_ch0 + V_C3AH6*n_c3ah60 ;
    
    N_CH(i)          = n_ch0 ;
    N_C3AH6(i)       = n_c3ah60 ;
    N_AFt(i)         = 0 ;
    
    PHI(i)           = phi0 ;
    
    {
      PoreRadius(i)    = PoreEntryRadiusMax ;
    }
    
    V_Cem0(i)        = v_cem ;
  }
  
  
  /* Initialization */
  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,u,f,0,i) ;

    /* Back up */

    /* Molar contents */
    N_S(i)           = x[I_N_S] ;

    /* Others */
    N_CH(i)          = x[I_N_CH] ;
    N_C3AH6(i)       = x[I_N_C3AH6] ;
    N_AFt(i)         = x[I_N_AFt] ;
    
    PoreRadius(i)    = x[I_Radius] ;
    
    PHI(i)           = x[I_PHI] ;
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el) ;

  /* Flux */
  {
    double* w = ComputeVariableFluxes(el,0,1) ;

    W_S     = w[I_W_S  ] ;
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms.
 *  IMPORTANT: if needed use only with the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  /* If you need the implicit terms, use the previous ones */
  double* f = Element_GetPreviousImplicitTerm(el) ;
  /* If you need the nodal values, use the previous ones */
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  
  /* If needed ! */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Compute here ve */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      ComputeVariables(el,u,u,f,0,i) ;
    }
  }
  
  ComputeTransferCoefficients(el) ;
  
  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* f   = Element_GetImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;

  
  
  for(i = 0 ; i < nn ; i++) {
    /* molarities */
    double* x = ComputeVariables(el,u,u_n,f_n,dt,i) ;

    /* Back up */

    /* Molar contents */
    N_S(i)  = x[I_N_S] ;

    N_CH(i)          = x[I_N_CH] ;
    N_C3AH6(i)       = x[I_N_C3AH6] ;
    N_AFt(i)         = x[I_N_AFt] ;
    
    PoreRadius(i)    = x[I_Radius] ;
    
    PHI(i)           = x[I_PHI] ;


    {
      int test = 0 ;
      if(x[I_PHI]     < 0.) test = -1 ;
      
      if(test < 0) {
        double xx = Element_GetNodeCoordinate(el,i)[0] ;
        
        printf("x         = %e\n",xx) ;
        printf("phi       = %e\n",x[I_PHI]) ;
        printf("c_so4     = %e\n",x[I_C_SO4]) ;
        printf("n_ch      = %e\n",x[I_N_CH]) ;
        return(-1) ;
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Flux */
  {
    double* w = ComputeVariableFluxes(el,0,1) ;

    W_S     = w[I_W_S  ] ;
  }

  return(0) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  int    i ;

  /*
    Initialization 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  {
    FVM_t* fvm = FVM_GetInstance(el) ;
    double c[4*NEQ*NEQ] ;
    
    TangentCoefficients(el,dt,c) ;
    
    {
      double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
      for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
    }
  }
  

#if (U_SO4 == NOLOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_SO4)     /= Ln10*C_SO4(0) ;
      K(i,U_C_SO4+NEQ) /= Ln10*C_SO4(1) ;
    }
  }
#endif
  

  
  return(0) ;
#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double* area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /*
    Mass balance of elements S
  */
  R(0,E_S)  -= volume[0]*(N_S(0)  - N_Sn(0))  + dt*surf*W_S ;
  R(1,E_S)  -= volume[1]*(N_S(1)  - N_Sn(1))  - dt*surf*W_S ;
  
  
  return(0) ;
#undef R
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  int NbOfOutputs  = 10 ;
  int    i ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;


  /* Initialization */
  for(i = 0 ; i < NbOfOutputs ; i++) {
    Result_SetValuesToZero(r + i) ;
  }
  
  
  {
    FVM_t* fvm = FVM_GetInstance(el) ;
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    double* x = ComputeVariables(el,u,u,f,0,j) ;

    i = 0 ;
    
    {
      double c_so4 = pow(10,x[I_C_SO4]) ;
      Result_Store(r + i++,&c_so4       ,"Sulfate concentration",1) ;
    }
    Result_Store(r + i++,(x + I_PHI    ),"Porosity",1) ;
    Result_Store(r + i++,(x + I_S_AFt  ),"Saturation index of AFt",1) ;
    Result_Store(r + i++,(x + I_S_C    ),"Saturation degree of crystal",1) ;
    Result_Store(r + i++,(x + I_N_S    ),"Sulfur mole content",1) ;
    Result_Store(r + i++,(x + I_N_CH   ),"Portlandite mole content",1) ;
    Result_Store(r + i++,(x + I_N_AFt  ),"Ettringite mole content",1) ;
    Result_Store(r + i++,(x + I_N_C3AH6),"Hydrogarnet mole content",1) ;
    Result_Store(r + i++,(x + I_Radius ),"Pore entry radius",1) ;
    {
      double beta = exp(2*Gamma_AFt*V_AFt/(RT*x[I_Radius])) ;
      Result_Store(r + i++,&beta,"Equilibrium saturation index of AFt",1) ;
    }
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}



int TangentCoefficients(Element_t* el,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    dec = NEQ*NEQ ;
  double dxi[NEQ] ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < NEQ ; i++) {
    dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
  }

  
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double* x         = ComputeVariables(el,u,u_n,f_n,dt,i) ;
    int k ;

    #if (U_SO4 == NOLOG_U)
    dxi[U_C_SO4] =  1.e-2*ObVal_GetValue(obval + U_C_SO4)/(Ln10*x[U_C_SO4]) ;
    #endif
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double* dx    = ComputeVariableDerivatives(el,dt,x,dxk,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
        
        cii[E_S*NEQ    + k] = dx[I_N_S] ;
      }
      
      /* Transfer terms from node i to node j */
      {
        int j = 1 - i ;
        double* dw    = ComputeFluxes(el,dx,i,j) ;
        
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            double* cij = c + (i*nn + j)*NEQ*NEQ ;
        
            cij[E_S*NEQ    + k] = - dt*dw[I_W_S] ;
          }
        }
      }
    }
  }

  return(dec) ;
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double dt,int n)
{
  double* x = Variables[n] ;
  
  /* Primary Variables */
  x[U_C_SO4 ]   = LogC_SO4(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_Radiusn]  = PoreRadiusn(n) ;
  x[I_PHIn   ]  = PHIn(n) ;
  
  {
    double* v0  = Element_GetConstantTerm(el) ;
    
    x[I_V_Cem0 ]  = V_Cem0(n) ;
  }
    
  ComputeSecondaryVariables(el,dt,x) ;
  return(x) ;
}





void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  /* Cement chemistry */
  double logc_so4  = x[U_C_SO4] ;
  double logc_h    = -ph ;
  double logc_oh   = LogK_H2O - logc_h ;
  double logc_ca   = LogK_CH - 2*logc_oh ;
  double logs_csh2 = logc_ca + logc_so4 - LogK_CSH2 ;
  double logs_aft  = - LogK_AFt + LogK_C3AH6 + 3*LogK_CSH2 + 3*logs_csh2 ;
    
  /* Compute the crystal saturation as a function of s_aft */
  double r_n = x[I_Radiusn] ;
  double s_aft  = pow(10,logs_aft) ;
  double r   = Radius(r_n,s_aft,dt,el) ;
  double s_l = LiquidSaturationDegree(r) ;
  double s_c = 1 - s_l ;

  /* Porosity (explicit scheme) */
  double phi = x[I_PHIn] ;
  
  /* Solid contents */
  /* ... as components: CH, AFt, C3AH6 */
  double n_aft      = phi*s_c/V_AFt ;
  double n_ch       = n_ch0 - 3*n_aft ;
  double n_c3ah6    = n_c3ah60 - n_aft ;
  /* ... as elements: S */
  double n_s_s      = 3*n_aft ;
  /* ... as volumes: S */
  double v_cem      = V_CH*n_ch + V_C3AH6*n_c3ah6 ;
    
    
  /* Liquid contents */
  /* ... as elements: S */
  double c_so4  = pow(10,logc_so4) ;
  double n_s_l  = phi*s_l*c_so4 ;


  /* Updated porosity */
  double v_cem0     = x[I_V_Cem0] ;
  phi = phi0 + v_cem0 - v_cem ;
  

  /* Back up */
  
  /* Log ion concentrations in liquid phase */
  x[I_C_SO4     ] = logc_so4 ;
  
  /* Saturation indexes of dissolved solids */
  x[I_S_AFt     ] = s_aft ;

  /* Solid components */
  x[I_N_CH      ] = n_ch ;
  x[I_N_AFt     ] = n_aft ;
  x[I_N_C3AH6   ] = n_c3ah6 ;
  
  /* Saturation degree of crystal */
  x[I_S_C       ] = s_c ;
  
  /* Porosity */
  x[I_PHI       ] = phi ;
  
  /* Radius */
  x[I_Radius    ] = r ;
  
  
  
  /* Element contents */
  x[I_N_S       ] = n_s_l  + n_s_s ;
}




double* ComputeVariableDerivatives(Element_t* el,double dt,double* x,double dxi,int i)
{
  double* dx = dVariables ;
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



void ComputeTransferCoefficients(Element_t* el)
/* Explicit terms (va)  */
{
  double* va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialization */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* molarities */
    double* x = Variables[i] ;
    
    double c_so4      = pow(10,x[I_C_SO4]) ;

    /* porosity */
    double phi        = x[I_PHI] ;
    
    /* tortuosity coefficient */
    double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
    iff *= Ln10 ;
    
    KF_SO4        += c_so4*d_so4*iff ;
  }
  
  
  {
    FVM_t* fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    double dij   = dist[1] ;
    double lij = 0.5/dij ;
  
    /* Averaging */
    for(i = 0 ; i < NVE ; i++) va[i] *= lij ;
  }
}



double* ComputeVariableFluxes(Element_t* el,int i,int j)
{
  double* grdij = dVariables ;


  {
    {
      double* xi  = Variables[i] ;
      double* xj  = Variables[j] ;
      int k ;
      
      for(k = 0 ; k < NbOfVariables ; k++)  {
        grdij[k] = xj[k] - xi[k] ;
      }
    }
  }
  
  /* Fluxes */
  {
    {
      double* w = ComputeFluxes(el,grdij,i,j) ;
    
      return(w) ;
    }
  }
}





double* ComputeFluxes(Element_t* el,double* grd,int i,int j)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes ;

  /* Gradients */
  double grd_so4      = grd[I_C_SO4] ;
    

  /* Flux */
  double w_so4      = - KF_SO4 * (grd_so4)  ;
  double w_s_l      =   w_so4 ;

  w[I_W_S ]  = w_s_l ;
   
  return(w) ;
}



double Radius(double r_n,double s_aft,double dt,Element_t* el)
{
  double r_max = PoreEntryRadiusMax ;
  double beta_min  = EquilibriumAFtSaturationIndex(r_max) ;
  double r_inf = (s_aft > beta_min) ? InverseOfEquilibriumAFtSaturationIndex(s_aft) : r_max ;
  double r = r_n ;
  
  {
    double s_ln  = LiquidSaturationDegree(r_n) ;
    int iterations = 20 ;
    double tol = 1.e-6 ;
    int i ;
    
    if(r_n == r_inf) return(r_n) ;
    
    for(i = 0 ; i < iterations ; i++) {
      double s_l   = LiquidSaturationDegree(r) ;
      double beta  = EquilibriumAFtSaturationIndex(r) ;
      double eq    = s_l - s_ln + dt*a_AFt*(1 - beta/s_aft) ;
      double ds_l  = dLiquidSaturationDegree(r) ;
      double dbeta = dEquilibriumAFtSaturationIndex(r) ;
      double deq   = ds_l - dt*a_AFt*dbeta/s_aft ;
      double dr    = (fabs(deq) > 0) ? - eq/deq : 0. ;
      
      /* The solution r should be in the range between r_n and r_inf.
       * So let us assume that, at a given iteration, an estimation r
       * has been found between r_n and r_inf. Then we look for a new 
       * estimation r + dr in the range between r_0 and r_1 by using
       * an under-relaxation technique so that dr should be given by 
       *         dr = a*(r_1 - r) + (1 - a)*(r_0 - r)
       * with a being in the range between 0 and 1.
       * The under-relaxation technique is such that r_0 should be in 
       * the range between r_n and r and r_1 should be in the range 
       * between r_inf and r i.e
       *         r_0 = b0*r_n   + (1 - b0)*r
       *         r_1 = b1*r_inf + (1 - b1)*r
       * with b0 and b1 being in between 0 and 1. So we get
       *         dr = a*b1*(r_inf - r) + (1 - a)*b0*(r_n - r)
       * If b0=b1=b then
       *         dr = (a*(r_inf - r) + (1 - a)*(r_n - r))*b
       * The bigger b the larger the range where r is looked for, without
       * exceding the initial range defined by r_n and r_inf.
       * The value b=0.5 corresponds to half the initial range.
       */
      {
        /* The predicted a is computed from the predicted dr */
        double b = 0.5 ;
        double a = (dr/b - r_n + r)/(r_inf - r_n) ;
       
        /* If the predicted a is < 0 then the used value is set to 0 */
        if(a < 0) a = 0 ;
      
        /* if the predicted a is > 1 then the used value is set to 1 */
        if(a > 1) a = 1 ;
        
        {
          dr = 0.5*(a*(r_inf - r) + (1 - a)*(r_n - r)) ;
        }
      }
      
      r += dr ;
      if(fabs(dr/r_n) < tol) return(r) ;
    }
  }
  
  Message_Warning("Radius: not converged") ;
  Exception_Interrupt ;

  return(r) ;
}
