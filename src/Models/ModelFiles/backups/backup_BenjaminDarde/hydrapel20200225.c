#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "hydration of expansive pellet mixtures (2019)"
#define AUTHORS "Darde"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (1+dim)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (32)
#define NVE     (2)
#define NV0     (0)

/* Equation index */
#define E_water  (0+dim)
#define E_mec    (0)

/* Unknown index */
#define U_p_l2   (0+dim)
#define U_u      (0)

/* We define some names for implicit terms */
#define M_water         (vim   + 0)[0]
#define M_water_n       (vim_n + 0)[0]

#define M_L1          (vim   + 1)[0]
#define M_L1_n        (vim_n + 1)[0]

#define W_water           (vim   + 2)

#define SIG           (vim   + 5)
#define SIG_n         (vim_n + 5)

#define F_MASS        (vim   + 14)

#define EPS_P         (vim   + 17)
#define EPS_P_n       (vim_n + 17)

#define HARDV         (vim   + 26)[0]
#define HARDV_n       (vim_n + 26)[0]

#define CRIT          (vim   + 28)[0]

#define EPSV_1        (vim   + 29)[0]
#define EPSV_1n       (vim_n + 29)[0]

#define EPSV_2        (vim   + 30)[0]
#define EPSV_2n       (vim_n + 30)[0]

#define P_L_1         (vim   + 31)[0]
#define P_L_1n        (vim_n + 31)[0]


/* We define some names for explicit terms */
#define K_L           (vex + 0)[0]
#define D_V           (vex + 1)[0]

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

static void    ComputePhysicoChemicalProperties(double) ;

//static double pie(double,double,Curve_t*) ;
//static double dpiesdpl(double,double,Curve_t*) ;

#define ComputeFunctionGradients(...)  Plasticity_ComputeFunctionGradients(plasty,__VA_ARGS__)
#define ReturnMapping(...)             Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)         Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define UpdateElastoplasticTensor(...) Plasticity_UpdateElastoplasticTensor(plasty,__VA_ARGS__)



/* Physical properties
 * -------------------*/
#define TEMPERATURE      (293.)      /* Temperature (K) */


/* Water properties
 * ----------------*/
/* Molar mass */
#define M_H2O     (18.e-3)
/* Molar volume of liquid water */
#define V_H2O     (18.e-6)
/* Mass density */
//#define MassDensityOfWaterVapor(p_v)   (p_v*M_H2O/RT)
#define MassDensityOfWaterVapor(p_l)   (VaporPressure(p_l)*M_H2O/RT)
#define dMassDensityOfWaterVapor(p_l)  (dVaporPressure(p_l)*M_H2O/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(p_l)          (exp(V_H2O/RT*(p_l - p_l0)))
#define VaporPressure(p_l)             (p_v0*RelativeHumidity(p_l))
#define LiquidPressure(hr)             (p_l0 + RT/V_H2O*(log(hr)))
#define dRelativeHumidity(p_l)         (V_H2O/RT*exp(V_H2O/RT*(p_l - p_l0)))
#define dVaporPressure(p_l)            (p_v0*dRelativeHumidity(p_l))


/* Material properties 
 * -------------------*/
//#define SATURATION_CURVE        (Element_GetCurve(el))
//#define SaturationDegree(pc)    (Curve_ComputeValue(SATURATION_CURVE,pc))
//#define dSaturationDegree(pc)   (Curve_ComputeDerivative(SATURATION_CURVE,pc))

#define RELATIVEPERM_CURVE                (Element_GetCurve(el) + 1)
#define RelativePermeabilityToLiquid(pc)  1.   //(Curve_ComputeValue(RELATIVEPERM_CURVE,pc))
#define TortuosityToGas(f)            (tortuosity)
#define IntrinsicPermeability(rhod)   (k_int_M)

//#define CAPIHARDENING_CURVE     (Element_GetCurve(el) + 2)
//#define CapillaryHardening(pc)  (Curve_ComputeValue(CAPIHARDENING_CURVE,pc))

#define LOADINGCOLLAPSEFACTOR_CURVE     (Element_GetCurve(el) + 0)
#define LoadingCollapseFactor(s)  (Curve_ComputeValue(LOADINGCOLLAPSEFACTOR_CURVE,s))

//#define EquivalentPressure(pl,pg)   (pie(pl,pg,SATURATION_CURVE))
//#define dEquivalentPressure(pl,pg)  (dpiesdpl(pl,pg,SATURATION_CURVE))





/* Parameters */
static double  gravity ;
static double  rho_s ;
static double* sig0 ;
static double  rho_l0 ;
static double  p_g = 0 ;
static double  k_int_m ;
static double  k_int_M ;
static double  mu_l ;
static double  kappa_m ;
static double  kappa_M ;
static double  kappa_s ;
static double  poisson ;
static double  alpha ;
static double  beta ;
//static double  lambda_M0 ;
static double  p_c ;
static double  p_co0 ;
static double  e0 ;
static double  p_vap_sat ;
static double  B_a ;
static double  M_R_T ;
static double  B_b ;
static double  A_m ;
static double  B_m ;
static double  p_l0 ;
static double  p_limit ;
static double  k_s ;
static double  tortuosity ;


/* MODIF */
static double  initial_pel_dry_dens ;
static double  p_poudre ;
static double  p_pellet ;
static double  initial_equivalent_radius ;
static double  beta_hydro_0 ;
static double  beta_hydro_a ;
static double  rho_solid ;

static bool    isGran = true ;



static Elasticity_t* elasty ;
static Plasticity_t* plasty ;


static double RT ;
static double D_av0 ;
static double  p_l0 ;
static double  p_v0 ;




#include "DiffusionCoefficientOfMoleculeInAir.h"
#include "WaterVaporPressure.h"
#include "AtmosphericPressure.h"
#include "WaterViscosity.h"
#include "AirViscosity.h"
#include "PhysicalConstant.h"

void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Air (m2/s) */
  D_av0   = DiffusionCoefficientOfMoleculeInAir(H2O,TK) ;
  
  /* Water vapor pressure */
  p_v0    = WaterVaporPressure(TK) ;
  
  /* Reference pressures */
  //p_l0    = AtmosphericPressure ;
  p_l0    = p_g ;
  
  /* Mass densities */
  rho_l0 = 1000. ;
  
  /* Viscosities */
  mu_l    = WaterViscosity(TK) ;
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
}



/* Variable indexes */
enum {
  I_U     = 0,
  I_P_L1  = I_U     + 3,
  I_P_L2,
  I_EPS,
  I_SIG   = I_EPS   + 9,
  I_EPS_P = I_SIG   + 9,
  I_Fmass = I_EPS_P + 9,
  I_M_water = I_Fmass + 3,
  I_M_L1,
  I_S_L,
  I_EPSV_1,
  I_EPSV_2,
  I_W_water,
  I_HARDV = I_W_water   + 3,
  I_CRIT,
  I_RHO_L,
  I_PHI_M,
  I_K_water,
  I_GRD_P_L2,
  I_GRD_RHO_V = I_GRD_P_L2 + 3,
  I_Last = I_GRD_RHO_V + 3
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
  } else if(!strcmp(s,"Poisson")) { 
    return (2) ;
  } else if(!strcmp(s,"rho_l"))      { 
    return (3) ;
  } else if(!strcmp(s,"k_int_m"))      { 
    return (4) ;
  } else if(!strcmp(s,"mu_l"))       { 
    return (5) ;
  } else if(!strcmp(s,"initial_liquid_pressure"))       {
    return (6) ;
  } else if(!strcmp(s,"initial_stress"))       {
    return(7) ;
  } else if(!strncmp(s,"initial_stress_",15))   {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(7 + 3*i + j) ;
    
    /* Cam-clay */
  } else if(!strcmp(s,"kappa_m")) {
    return(16) ;
  } else if(!strcmp(s,"kappa_M")) {
    return(17) ;
  } else if(!strcmp(s,"kappa_s"))  {
    return(18) ;
  } else if(!strcmp(s,"initial_pre-consolidation_pressure")) {
    return(19) ;
  } else if(!strcmp(s,"initial_void_ratio")) {
    return(20) ;
  } else if(!strcmp(s,"p_c")) {
    return(21) ;
  } else if(!strcmp(s,"limit_liquid_pressure")) {
    return(22) ;
  } else if(!strcmp(s,"p_vap_sat")) {
    return(23) ;
  } else if(!strcmp(s,"B_a")) {
    return(24) ;
  } else if(!strcmp(s,"M_R_T")) {
    return(25) ;
  } else if(!strcmp(s,"B_b")) {
    return(26) ;
  } else if(!strcmp(s,"A_m")) {
    return(27) ;
  } else if(!strcmp(s,"B_m")) {
    return(28) ;
  } else if(!strcmp(s,"slope_of_critical_state_line"))  {
    return(29) ;
  } else if(!strcmp(s,"k_s"))  {
    return(30) ;
  } else if(!strcmp(s,"alpha"))  {
    return(31) ;
  } else if(!strcmp(s,"beta"))  {
    return(32) ;
  } else if(!strcmp(s,"lambda_M0"))  {
    return(33) ;

	
  } else if(!strcmp(s,"initial_pel_dry_dens"))  {
    return(34) ;
  } else if(!strcmp(s,"p_poudre"))  {
    return(35) ;
  } else if(!strcmp(s,"p_pellet"))  {
    return(36) ;	
  } else if(!strcmp(s,"initial_equivalent_radius"))  {
    return(37) ;	
  } else if(!strcmp(s,"beta_hydro_0"))  {
    return(38) ;	
  } else if(!strcmp(s,"beta_hydro_a"))  {
    return(39) ;
  } else if(!strcmp(s,"rho_solid"))  {
    return(40) ;
  } else if(!strcmp(s,"tortuosity"))  {
    return(41) ;
  } else if(!strcmp(s,"k_int_M"))      { 
    return (42) ;
	
	} else return(-1) ;
}

/* MODIFS :
static double  initial_pel_dry_dens ;
static double  p_poudre ;
static double  p_pellet ;
static double  initial_equivalent_radius ;
static double  beta_hydro_0 ;
static double  beta_hydro_a ;
*/



void GetProperties(Element_t* el)
{
  gravity     = Element_GetPropertyValue(el,"gravity") ;
  rho_s       = Element_GetPropertyValue(el,"rho_s") ;
  k_int_m       = Element_GetPropertyValue(el,"k_int_m") ;
  k_int_M       = Element_GetPropertyValue(el,"k_int_M") ;
  mu_l        = Element_GetPropertyValue(el,"mu_l") ;
  rho_l0      = Element_GetPropertyValue(el,"rho_l") ;
  sig0        = &Element_GetPropertyValue(el,"initial_stress") ;
  p_co0       = Element_GetPropertyValue(el,"initial_pre-consolidation_pressure") ;
  e0          = Element_GetPropertyValue(el,"initial_void_ratio") ;
  kappa_m     = Element_GetPropertyValue(el,"kappa_m") ;
  kappa_M     = Element_GetPropertyValue(el,"kappa_M") ;
  kappa_s     = Element_GetPropertyValue(el,"kappa_s") ;
  poisson     = Element_GetPropertyValue(el,"Poisson") ;
  alpha       = Element_GetPropertyValue(el,"alpha") ;
  beta        = Element_GetPropertyValue(el,"beta") ;
//  lambda_M0 = Element_GetPropertyValue(el,"lambda_M0") ;
  p_c         = Element_GetPropertyValue(el,"p_c") ;
  p_vap_sat         = Element_GetPropertyValue(el,"p_vap_sat") ;
  B_a         = Element_GetPropertyValue(el,"B_a") ;
  M_R_T         = Element_GetPropertyValue(el,"M_R_T") ;
  B_b         = Element_GetPropertyValue(el,"B_b") ;
  A_m         = Element_GetPropertyValue(el,"A_m") ;
  B_m         = Element_GetPropertyValue(el,"B_m") ;
  p_l0        = Element_GetPropertyValue(el,"initial_liquid_pressure") ;
  p_limit     = Element_GetPropertyValue(el,"limit_liquid_pressure") ;
  k_s         = Element_GetPropertyValue(el,"k_s") ;
  initial_pel_dry_dens = Element_GetPropertyValue(el,"initial_pel_dry_dens") ;
  p_poudre    = Element_GetPropertyValue(el,"p_poudre") ;
  p_pellet    = Element_GetPropertyValue(el,"p_pellet") ;
  initial_equivalent_radius = Element_GetPropertyValue(el,"initial_equivalent_radius") ;
  beta_hydro_0= Element_GetPropertyValue(el,"beta_hydro_0") ;
  beta_hydro_a= Element_GetPropertyValue(el,"beta_hydro_a") ;
  rho_solid= Element_GetPropertyValue(el,"rho_solid") ;
  tortuosity = Element_GetPropertyValue(el,"tortuosity") ;

  plasty  = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  elasty  = Plasticity_GetElasticity(plasty) ;
}




static double calcul_m(double,double,double) ;
static double calcul_delta_epsilon_m(double,double,double,double) ;
static double calcul_delta_epsilon_pel(double,double,double,double,double) ;
static double calcul_m_lim(void) ;
static double calcul_dFdm(double,double,double) ;
static double calcul_dmdp(double,double,double) ;
static double calcul_deps1dp(double,double,double) ;
static double calcul_deps1ds(double,double,double) ;
static double calcul_dmds(double,double,double) ;
static double calcul_kappa_eq(double,double,double,double,double,double) ;
static double calcul_delta_epsilon_sigm(double,double,double,double,double,double) ;
static double calcul_delta_epsilon_s1(double,double,double,double,double) ;
static double calcul_delta_epsilon_s2(double,double,double,double,double) ;
static double calcul_de1dp(double,double,double,double) ;
static double calcul_de2dp(double,double,double,double) ;
static double calcul_deMdp(double) ;
static double calcul_delta_epsilon_pel_cont(double,double,double,double,double,double) ;
static double calcul_delta_epsilon_pou_cont(double,double,double,double,double,double) ;
static double bulkModulusGran(double, double, double) ;
static double bulkModulusCont(double) ;
static double bulkModulus(double*, double*, double, double, double, double) ;


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
  Model_CopyNameOfEquation(model,E_water,"water") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_p_l2,"p_l2") ;
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
/* MODIF int  NbOfProp = 34 ; */
  int  NbOfProp = 43 ;
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
        //double young    = Material_GetPropertyValue(mat,"young") ;
        //double poisson  = Material_GetPropertyValue(mat,"poisson") ;
        
        Elasticity_SetToIsotropy(elasty) ;
        //Elasticity_SetParameters(elasty,young,poisson) ;
        //Elasticity_ComputeStiffnessTensor(elasty) ;
      }
    }
    {
      /* Cam-Clay with offset*/
      {
        double lambda_M0 = Material_GetPropertyValue(mat,"lambda_M0") ;
        double M      = Material_GetPropertyValue(mat,"slope_of_critical_state_line") ;
        double pc0    = Material_GetPropertyValue(mat,"initial_pre-consolidation_pressure") ;
        
        e0     = Material_GetPropertyValue(mat,"initial_void_ratio") ;
        kappa_M  = Material_GetPropertyValue(mat,"kappa_M") ;
        
        Plasticity_SetToCamClayOffset(plasty) ;
        Plasticity_SetParameters(plasty,kappa_M,lambda_M0,M,pc0,e0) ; // verifier e0
      }
    }

#if 0
    {      
      Elasticity_PrintStiffnessTensor(elasty) ;
    }
#endif
  }

  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
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
  fprintf(ficd,"k_int_m = 1e-19     # intrinsic permeability\n") ;
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_water) {
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
  int dim = Element_GetDimensionOfSpace(el) ;
  FEM_t* fem = FEM_GetInstance(el) ;
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
      double p_l2 = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l2) ;
      double* vim  = vim0 + p*NVI ;
      int    i ;

      
      /* Initial stresses, hardening variables */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
      } else {
        for(i = 0 ; i < 9 ; i++) SIG[i] = sig0[i] ;
        HARDV = p_co0 ;
      }
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
      
      EPSV_1 = 0. ;
      EPSV_2 = 0. ;
      P_L_1  = p_l2 ;
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
      
      M_water = x[I_M_water] ;
    
      for(i = 0 ; i < 3 ; i++) W_water[i] = x[I_W_water + i] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT  = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
      
      EPSV_1 = x[I_EPSV_1] ;
      EPSV_2 = x[I_EPSV_2] ;
      P_L_1  = x[I_P_L1] ;
    }
    
    
    /* storage in vex */
    {
      /* pressures */
      double p_l2 = x[I_P_L2] ;
      double pc = p_g - p_l2 ;
    
      /* saturation degrees */
      double s_l = x[I_S_L] ;
      double s_g = 1 - s_l ;
    
      /* porosity */
      double phi_M = x[I_PHI_M] ;
    
      /* tortuosity */
      double tau_g   = TortuosityToGas(phi_M) ;
      double* vex  = vex0 + p*NVE ;
      double rho_l = x[I_RHO_L] ;
      double rho_d = 0 ;
      double kint  = IntrinsicPermeability(rho_d) ;
      double k_h = rho_l*kint/mu_l*RelativePermeabilityToLiquid(pc) ;
      double d_v = phi_M * s_g * tau_g * D_av0 ;
    
      K_L = k_h ;
      D_V = d_v ;
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
    double p_l2 = x[I_P_L2] ;
    double pc = p_g - p_l2 ;
    
    /* saturation degrees */
    double s_l = x[I_S_L] ;
    double s_g = 1 - s_l ;
    
    /* porosity */
    double phi_M = x[I_PHI_M] ;
    
    /* tortuosity */
    double tau_g   = TortuosityToGas(phi_M) ;
    
    /* permeability */
    double rho_d = 0 ;
    double kint = IntrinsicPermeability(rho_d) ;
    double k_h = rho_l*kint/mu_l*RelativePermeabilityToLiquid(pc) ;
    double d_v = phi_M * s_g * tau_g * D_av0 ;
    
    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      
      K_L = k_h ;
      D_V = d_v ;
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
      
      M_water = x[I_M_water] ;
    
      for(i = 0 ; i < 3 ; i++) W_water[i] = x[I_W_water + i] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;

      EPSV_1 = x[I_EPSV_1] ;
      EPSV_2 = x[I_EPSV_2] ; 
      P_L_1  = x[I_P_L1] ;
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
    int dec = ComputeTangentCoefficients(fem,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    /* The matrix kp is stored as (u for displacement, s1,s2 for pressure)
     * | Kuu  Kup2  Kup1  |
     * | Kp2u Kp2p2 Kp2p1 |
     * | Kp1u Kp1p2 Kp1p1 |
     * i.e. the displacements u are in the positions 0 to dim-1 and
     * the pressure p is in the position dim.
     * So we need to store the matrix by accounting for the right indexes.
     */
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
        K(E_water + i*NEQ,U_p_l2 + j*NEQ) += dt*kc[i*nn + j] ;
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
  
  
  /* 2. Conservation of the mass of water */
  
  /* 2.1 Accumulation Terms */
  {
    double* vim = vim_1 ;
    double* vim_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_water - M_water_n ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_water) -= ra[i] ;
    }
  }
  
  /* 2.2 Transport Terms */
  {
    double* vim = vim_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_water,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_water) -= -dt*rf[i] ;
  }
  
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 8 ;
  double* vex  = Element_GetExplicitTerm(el) ;
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
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
    /* Variables */
    //double* x = ComputeVariables(el,u,u,vim,t,0,p) ;
    /* Pressures */
    double p_l1 = 0 ;
    double p_l2 = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l2) ;
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
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      double* def =  FEM_ComputeLinearStrainTensor(fem,u,intfct,i,U_u) ;
      int j ;
      
      p_l1 += P_L_1/np ;
      
      for(j = 0 ; j < 3 ; j++) w_l[j] += W_water[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps[j] += def[j]/np ;
      
      hardv += HARDV/np ;
      
      k_h += K_L/np ;
    }
    
    tre = eps[0] + eps[4] + eps[8] ;
    e   = (1 + e0) * tre ;
      
    i = 0 ;
    Result_Store(r + i++,&p_l1    ,"Pore pressure in pellets",1) ;
    Result_Store(r + i++,&p_l2    ,"Pore pressure in powder",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_l      ,"Fluid mass flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,&e       ,"Void ratio variation",1) ;
    Result_Store(r + i++,eps_p    ,"Plastic strains",9) ;
    Result_Store(r + i++,&hardv   ,"Hardening variable",1) ;
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
  double*  vim0   = Element_GetCurrentImplicitTerm(el) ;
  double*  vim0_n = Element_GetPreviousImplicitTerm(el) ;
//  double*  vex0  = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double dui[NEQ] ;
  int    dec = 100 ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  {
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }


  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < np ; p++) {
      double* vim   = vim0   + p*NVI ;
      double* vim_n = vim0_n + p*NVI ;
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vim0_n,t,dt,p) ;
      double* c0 = c + p*dec ;


      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
    

      /* Mechanics */
      {
    
        /* Tangent stiffness matrix */
        {
          double* c1 = c0 ;
        
          {
            double p_l1 = x[I_P_L1] ;
            double p_l2 = x[I_P_L2] ;

            double* sig   = x + I_SIG ; // contraintes

            double* eps   =  x   + I_EPS ;
            double eps_1  = x[I_EPSV_1] ; // deformation pellets
            double eps_2  = x[I_EPSV_2] ; // deformation poudre

            double bulk = bulkModulus(sig,eps,eps_1,eps_2,p_l1,p_l2);

            double young = 3 * bulk * (1 - 2*poisson) ;
          
            Elasticity_SetParameters(elasty,young,poisson) ;
          }

          Elasticity_ComputeStiffnessTensor(elasty,c1) ;
      
          {
            double crit = CRIT ;
            
            /* Criterion */
            if(crit >= 0.) {
              double p_l2 = x[I_P_L2] ;
              double s2 = p_g - p_l2 ;
              double p_con = HARDV_n ;
              double lcf   = LoadingCollapseFactor(s2) ;
              double logp_co = log(p_c) + log(p_con/p_c) * lcf ;
              double p_co    = exp(logp_co) ;
              double p_s     = k_s * s2 ;
              double hardv[2] = {p_co,p_s} ;
              double* sig = SIG ;
              double crit1 = ComputeFunctionGradients(sig,hardv) ;
              double fcg   = UpdateElastoplasticTensor(c1) ;
          
              if(fcg < 0) return(-1) ;
            }
          }
        }

      
        /* Coupling matrix (p_l2) */
        {
          double  dp_l2 = dui[U_p_l2] ;
          double* dxdp_l2 = ComputeVariableDerivatives(el,t,dt,x,dp_l2,I_P_L2) ;
          double* dsigdp_l2 = dxdp_l2 + I_SIG ;
          double  dtrsigdp_l2 = dsigdp_l2[0] + dsigdp_l2[4] + dsigdp_l2[8] ;
          double  biot2 = - dtrsigdp_l2 / 3 ;
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 3 ; i++) B1(i,i) = biot2 ;

        }
      }
    
    
      /* Conservation of the mass of water */
      {
    
        /* Coupling matrix */
        {
          double* c1 = c0 + 81 + 9 ;
          double deps = 1.e-6 ;
          //double* dxdeps = ComputeVariableDerivatives(el,t,dt,x,deps,I_EPS) ;
          //double dm_wdeps = dxdeps[I_EPS] ;
          double dm_wdeps = 0 ;
          int i ;

          for(i = 0 ; i < 3 ; i++) B1(i,i) = dm_wdeps ;
        }
      
      
        /* Storage matrix */
        {
          double  dp_l2 = dui[U_p_l2] ;
          double* dxdp_l2 = ComputeVariableDerivatives(el,t,dt,x,dp_l2,I_P_L2) ;
          double  dmdp_l2 = dxdp_l2[I_M_water] ;
          double* c1 = c0 + 81 + 9 + 9 ;
        
          c1[0] = dmdp_l2 ;
        }
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
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
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
      double pl = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l2) ;
      
      /* Permeability tensor */
      c1[0] = K_L + D_V * dMassDensityOfWaterVapor(pl) ;
      c1[4] = K_L + D_V * dMassDensityOfWaterVapor(pl) ;
      c1[8] = K_L + D_V * dMassDensityOfWaterVapor(pl) ;
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
    
    /* Pressures */
    x[I_P_L2] = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l2) ;
    
    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_p_l2) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P_L2 + i]  = grd[i] ;
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
    
    /* Pressures at previous time step */
    {
      double* vim_n = f_n + p*NVI ;
      
      x_n[I_P_L2] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l2) ;
      x_n[I_P_L1] = P_L_1n ;
    }
    
    /* Microscopic volumetric strains (pellets and grains) */
    {
      double* vim_n = f_n + p*NVI ;
      
      x_n[I_EPSV_1] = EPSV_1n ;
      x_n[I_EPSV_2] = EPSV_2n ;
    }
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;
      double pl = x_n[I_P_L2] ;
      
      x_n[I_K_water]  = K_L + D_V * dMassDensityOfWaterVapor(pl) ;
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
  /* Pressures */
  double  p_l1n  = x_n[I_P_L1] ;
  double  p_l2   = x[I_P_L2] ;
  double  p_l2n  = x_n[I_P_L2] ;
    

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
      
      /* Stresses at t_n */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
    
      /* Elastic trial stresses at t */
      {
        
        double tr_sig_n  = sig_n[0] + sig_n[4] + sig_n[8] ;
        double sigm_n    = tr_sig_n/3. ;
        double tr_eps_n  = eps_n[0] + eps_n[4] + eps_n[8] ;
        double eps_1_n   = x_n[I_EPSV_1] ;
        double eps_2_n   = x_n[I_EPSV_2] ;

        double initial_solid_fraction_1 = p_pellet*rho_s/initial_pel_dry_dens ; // compacités
        double initial_solid_fraction_2 = p_poudre*rho_s/initial_pel_dry_dens ;
        double solid_fraction_1         = initial_solid_fraction_1 * (1.+eps_1_n)/(1.+ tr_eps_n ) ;
        double solid_fraction_2         = initial_solid_fraction_2 * (1.+eps_2_n)/(1.+ tr_eps_n ) ;

// ---------------------------------------------       vérification des différentes conditions 
        if(p_l2 - p_g >= p_limit) {
// le mélange est "continu" car succion inf succion limite
             isGran    = false ;
        }

        if((solid_fraction_2 / (1.-solid_fraction_1)) >= solid_fraction_1) {
// le mélange est "continu" car compacite matrice = compacite assemblage de pellets
             isGran    = false ;
        }
// --------------------------------------------- 

  #if 1
// --------------------------------------------- calcul p_l1
        {
/* calcul temps caractéristique de diffusion */
          double deps1ds_n = calcul_deps1ds(p_l1n,sigm_n,solid_fraction_1) ; // > 0
          double C_diff_n = k_int_m / (mu_l * deps1ds_n) ; // > 0
          double tau_hydro_n = initial_equivalent_radius * initial_equivalent_radius * (pow((1+eps_1_n),2./3.)) / C_diff_n ; // temps caractéristique de diffusion au temps t_n

/* calcul beta, paramètre de transfert micro-macro */
          double beta_hydro_n = beta_hydro_0*pow((-p_l1n),beta_hydro_a) ;

/* calcul dmw/ds1 au temps t_n */
          double dmwds1_n = - rho_l0*initial_solid_fraction_1*deps1ds_n ; // s diminue, mw augmente

/* calcul succion micro */
          double C_transfert = dt / (dt - dmwds1_n * tau_hydro_n/beta_hydro_n) ;
          double p_l1 = p_l1n + C_transfert * (p_l2-p_l1n) ;
/*
          printf("p_l1 calcul. %e\n",p_l1) ;
          printf("C_transfert calcul. %e\n",C_transfert) ;
          printf("C_diff_n. %e\n",C_diff_n) ;
          printf("deps1ds_n. %e\n",deps1ds_n) ;
*/          
          x[I_P_L1] = p_l1 ;
          //x[I_P_L1] = p_l2 ;
        }
// --------------------------------------------- 
  #endif

        {
// ---------------------------------------------   calcul d_sigm par methode de Newton / sig_i+1 - sig_i = dsig_i = - f(sig_i) / ( df(sig_i)/dsig_i )
          double tr_deps = deps[0]+deps[4]+deps[8] ;
          double sigm = sigm_n ; // initialisation avant iterations
          double bulk ;
          double p_l1 = x[I_P_L1] ;
        
          if(isGran == true) {
// calcul granulaire

             double m_n = calcul_m(sigm_n,solid_fraction_1,p_l1n) ;
             double f_Mm = 1. ; // deformation variation volume pellets
             
             double tol = 1.e-8 ;

             int it = 0 ;
             do
             {
                 it = it + 1 ;

                 double m = calcul_m(sigm,solid_fraction_1,p_l1) ;

                 double Delta_F_i = calcul_delta_epsilon_m(m,m_n,B_a,B_b) ; // F_i - F_n ; // deformation variation de m
                 double Delta_G_i = calcul_delta_epsilon_pel(sigm,solid_fraction_1,p_l1,sigm_n,p_l1n) ; // G_i - G_n ; // deformation variation volume pellets

                 double f_test = tr_deps - Delta_F_i - f_Mm*Delta_G_i ;
                     
                 double dFdm        = calcul_dFdm(m,B_a,B_b); // " d_rond epsilon_V / d_rond m "
                 double dmdp        = calcul_dmdp(p_l1,sigm,solid_fraction_1) ; // " d_rond m / d_rond p "
                 double deps1dp     = calcul_deps1dp(p_l1,sigm,solid_fraction_1) ; // "d_rond epsilon_V1 / d_rond p "

                 double d_f_test = - dFdm * dmdp - f_Mm * deps1dp ; // "d f_test / dp "

                 double tmp_dsigm = (d_f_test != 0.) ? - f_test / d_f_test : 0 ;

                 sigm = sigm + tmp_dsigm ;
                     
                 if (fabs(f_test) <= tol) break ;

             } while ( it < 50 ) ;
             
             if(it >= 50) {
                arret("No convergence") ;
              }
             
             {
               double m       = calcul_m(sigm,solid_fraction_1,p_l1) ;
               double dFdm    = calcul_dFdm(m,B_a,B_b);
               double dmdp    = calcul_dmdp(p_l1,sigm,solid_fraction_1) ; // " d_rond m / d_rond p "
               double deps1dp = calcul_deps1dp(p_l1,sigm,solid_fraction_1) ; // "d_rond epsilon_V1 / d_rond p "

               bulk = bulkModulusGran(dFdm, dmdp, deps1dp);
             }

          } else {
// calcul continu
    
             double Delta_G_i = calcul_delta_epsilon_s1(sigm_n,p_l1,p_l1n,solid_fraction_1,solid_fraction_2) ; // deformation due aux variations de s1
             double Delta_H_i = calcul_delta_epsilon_s2(sigm_n,p_l2,p_l2n,solid_fraction_1,solid_fraction_2) ; // deformation due aux variations de s2

             double tol = 1.e-8 ;

             int it = 0 ;
             do
             {
                 it = it + 1 ;

                 double Delta_F_i = calcul_delta_epsilon_sigm(sigm,sigm_n,p_l1n,p_l2n,solid_fraction_1,solid_fraction_2) ;

                 double f_test = tr_deps - Delta_F_i - Delta_G_i - Delta_H_i ;
                 
                 double de1dp = calcul_de1dp(sigm,p_l1,solid_fraction_1,solid_fraction_2) ;

                 double de2dp = calcul_de2dp(sigm,p_l2,solid_fraction_1,solid_fraction_2) ;
                 double deMdp = calcul_deMdp(sigm) ;
                 
                 double d_f_test = - de1dp - de2dp - deMdp ;
                 
                 double tmp_dsigm = (d_f_test != 0.) ? - f_test / d_f_test : 0 ;
                 sigm = sigm + tmp_dsigm ;
                 
                 if (fabs(f_test) <= tol) break ;

             } while ( it < 50 ) ;
             
             if(it >= 50) {
                arret("No convergence") ;
              }
             
             {
               double propPel = p_pellet / (p_pellet + p_poudre) ;
               double kappa_eq = calcul_kappa_eq(propPel,sigm,p_l1,p_l2,solid_fraction_1,solid_fraction_2) ;
               // double initial_total_vr = rho_solid / rho_s - 1 ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange

               bulk = bulkModulusCont(kappa_eq);
             }
             
          }
// ---------------------------------------------   fin calcul d_sigm par methode de Newton



// ---------------------------------------------       calcul de sigma
          {
            double mu_elas = 3./2.*(1.-2.*poisson)/(1.+poisson)*bulk ; // G global
            double dsigm = sigm - sigm_n ;

            for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] + 2*mu_elas*deps[i] ;
      
            sig[0] = sig[0] + dsigm - 2./3.*mu_elas*tr_deps ;
            sig[4] = sig[4] + dsigm - 2./3.*mu_elas*tr_deps ;
            sig[8] = sig[8] + dsigm - 2./3.*mu_elas*tr_deps ; // sig_ij = 2G eps_ij + (  (K-2/3G) tr_eps - biot s  ) * I, avec K tr_eps - biot s = sigm
          }
// ---------------------------------------------
        }


      }
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Return mapping */
      {
        double s2    = p_g - p_l2 ;
        double p_con = x_n[I_HARDV] ;
        double lcf   = LoadingCollapseFactor(s2) ;
        double logp_co = log(p_c) + log(p_con/p_c) * lcf ;
        double p_co    = exp(logp_co) ;
        double p_s     = k_s * s2 ;
        double hardv[2] = {p_co,p_s} ;
        double crit  = ReturnMapping(sig,eps_p,hardv) ;
        
        p_co = p_c * exp(log(hardv[0]/p_c) / lcf) ;
        
        x[I_CRIT]  = crit ;
        x[I_HARDV] = p_co ;
      }
    }
  }
  
  /* Backup microscopic volumetric strains epsv1 and epsv2 */
  {
  
// ---------------------------------------------       sigma et epsilon 
    double* sig   = x + I_SIG ; // contraintes
    double sigm   = (sig[0] + sig[4] + sig[8])/3. ;
    double sigm_n = (sig_n[0] + sig_n[4] + sig_n[8])/3. ;
  
    double tr_eps = eps[0] + eps[4] + eps[8] ; // deformation volumique mélange
    double eps_1n = x_n[I_EPSV_1] ; // deformation pellets t_n
    double eps_2n = x_n[I_EPSV_2] ; // deformation poudre t_n
    double delta_eps_V1 ; // deformation volumique pellets
    double delta_eps_V2 ; // deformation volumique poudre
// --------------------------------------------- 


// ---------------------------------------------       fractions solides ("compacités") 
    double initial_solid_fraction_1 = p_pellet*rho_s/initial_pel_dry_dens ; // compacités
    double initial_solid_fraction_2 = p_poudre*rho_s/initial_pel_dry_dens ;
    double solid_fraction_1         = initial_solid_fraction_1 * (1.+eps_1n)/(1.+ tr_eps ) ;
    double solid_fraction_2         = initial_solid_fraction_2 * (1.+eps_2n)/(1.+ tr_eps ) ;
    double p_l1 = x[I_P_L1] ;
// ---------------------------------------------
    
 
// ---------------------------------------------       contraintes effectives et d epsilon
    if(isGran == true) {
         double sigm2 = 0. ;
         double sigm2_n = 0. ;

         delta_eps_V1 = calcul_delta_epsilon_pel(sigm,solid_fraction_1,p_l1,sigm_n,p_l1n) ;
         delta_eps_V2 = calcul_delta_epsilon_pel(sigm2,solid_fraction_2,p_l2,sigm2_n,p_l2n) ; // meme expression mais avec sig = 0 pour la poudre dans la partie "granulaire"
    } else {
         delta_eps_V1 = calcul_delta_epsilon_pel_cont(sigm,sigm_n,p_l1,p_l1n,solid_fraction_1,solid_fraction_2)  ;
         delta_eps_V2 = calcul_delta_epsilon_pou_cont(sigm,sigm_n,p_l2,p_l2n,solid_fraction_1,solid_fraction_2) ;
    }
// --------------------------------------------- 


// ---------------------------------------------       increments déformations
    x[I_EPSV_1] = eps_1n + delta_eps_V1 ; // delta epsilon V1
    x[I_EPSV_2] = eps_2n + delta_eps_V2 ; // delta epsilon V2

// --------------------------------------------- 

  }
  
  
  /* Backup mass flow */
  #if 1
  {
    /* Porosity */
    
    double tre  = eps[0] + eps[4] + eps[8] ;
    double eps1 = x[I_EPSV_1] ;
    double eps2 = x[I_EPSV_2] ;
    
    double initial_solid_fraction_1 = p_pellet*rho_s/initial_pel_dry_dens ; // compacités
    double initial_solid_fraction_2 = p_poudre*rho_s/initial_pel_dry_dens ;
    double solid_fraction_1         = initial_solid_fraction_1 * (1.+eps1)/(1.+ tre ) ;
    double solid_fraction_2         = initial_solid_fraction_2 * (1.+eps2)/(1.+ tre ) ;
    
    double initial_e1          = rho_solid / initial_pel_dry_dens - 1. ; // indice des vides a partir de densite mesuree au labo
    double initial_poro_1      = (initial_e1/(1.+initial_e1)) ; // porosite ini pellet
    
    double phi_M      = 1. - solid_fraction_1 - solid_fraction_2 ; // porosite macro = VvM/V = (V - V1 - V2)/V = 1 - sol.frac.1 - sol.frac.2
    double phi_m1   = (initial_poro_1 + eps1) / (1.+eps1) ; // Vv1 / V1, multiplier par sol.frac.1 pour rapporter au volume global
    double phi_m2   = (initial_poro_1 + eps2) / (1.+eps2) ; // idem pour Vv2/V2; densites poudre et pel identiques donc porosite ini identique
    double phi_tot   = phi_M + phi_m1*solid_fraction_1 + phi_m2*solid_fraction_2 ; // Vv_M / V   +   Vv_m1 / V1 * V1 / V   +   Vv_m2 / V2 * V2 / V
    
    /* Fluid mass density */
    double rho_l = rho_l0 ;
    
    /* Fluid mass flow */
    {
//  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      /* Transfer coefficient */ 
      double k_w = x_n[I_K_water] ;

	  
	  
	  
// FICK, besoin de tortuosite, phi_eq = VM2/(Vm2+VM2) ok via solid fractions, 1-SliqM ok, Dv, Ma masse molaire air environ 0.029 kg/ mol 
// rho_g = Ma*p_a / RT + M_R_T * pvapM
// pg = pa + pvap donc pour notre cas, pa = - pvap OU 0.1 - pvap
// fick : flux = - tortuosite * phi_eq * (1 - Sw) * Dv * rho_g grad(rho_vap / rho_g)
	  
	  
	  
// DARCY : flux = - K_int_grain_de_poudre * krw / mu_l * grad(succion2)
// K_int = B * exp ( - A * rho_d_matrice )
// krw = Swmat ^ "3 à 10" sorti du chapeau
// Swmat = omega_v2 / (omega 2 + omegaM)  * SwM * omegaM/ (omega 2 + omegaM) = phi_m2*solFrac2 / (solFrac2+1-solFrac1) * SwM * (1 - phi1-phi2)/(----------)
	  
	  
	  
//  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      /* Pressure gradient */
      double* gpl = x + I_GRD_P_L2 ;
    
      /* Mass flow */
      double* w_water = x + I_W_water ;
      int i ;
    
      for(i = 0 ; i < 3 ; i++) w_water[i] = - k_w*gpl[i] ;
    }
    
    /* Total liquid mass content and body force */
    {
      double  sat_liq_M = 0. ; //saturation liquide dans les macropores. fixée à 0. ajouter eventuellement une wrc si besoin.
      double  sat_gas_M = 1. - sat_liq_M ; //saturation vapeur eau dans les macropores.

      double  rho_vap = MassDensityOfWaterVapor(p_l2) ;

      double  m_w_M  = rho_l * phi_M * sat_liq_M + rho_vap * phi_M * sat_gas_M ; // dans les macropores ;
      double  m_w_1  = rho_l * phi_m1 * solid_fraction_1 ; // dans les pellets, revient à rho_l * (initial_poro + eps1)/(1+eps1) * initial_sol_frac * (1+eps1) / (1+epsV) = rho_l * (n_0 + eps1) * phi_1_0 / (1+tr_eps)
      double  m_w_2  = rho_l * phi_m2 * solid_fraction_2 ; // dans les grains de poudre 

      double  m_w    = m_w_1 + m_w_2 + m_w_M ;

      double* f_mass = x + I_Fmass ;
      int i ;
    
      x[I_S_L]     = sat_liq_M ;
      x[I_M_water] = m_w ;
      x[I_RHO_L]   = rho_l ;
      x[I_PHI_M]   = phi_M ; // porosite macro ici
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_w)*gravity ;
    }
    
    /* Liquid mass content in pellets */
    {
      double m_w_1 = rho_l * phi_m1 * solid_fraction_1 ; // m_w = rho * n * Phi [kg pour un m3 de matériau]  =  rho  *  (eps1 + n0)/(1+eps1)  *  Phi  
      x[I_M_L1]    = m_w_1 ;
    }
  }
  #endif
}






#if 1
double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double*  x_n = Variable_n ;
  double* dx = dVariable ;
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







double calcul_m(double sigm,double phi_1,double p_l1)
{
/* 
parametre m
*/
  double sigm1_eff = sigm/phi_1 + p_l1 ; // sigm <0 compression ; s < 0
  double bulk1 = exp(alpha*(-sigm1_eff))/beta ;
  double y_mod_eff_1 = bulk1 * 3.*(1. - 2.*poisson) / (1. - poisson*poisson) ;

  double m = (sigm < 0) ? - sigm / y_mod_eff_1 : 0 ;
  
  m        = pow(m,2./3.) ;
  
  return(m) ;
}



double calcul_delta_epsilon_m(double m,double m_n,double f_epsilon_a,double f_epsilon_b)
{
/* 
m = (p / E/(1-nu*nu) )^(2/3)
f_epsilon = delta epsilon_V / delta m

epsilon_V (m) = bilinear ; = f_epsilon_a*m (m<m_lim) ; = f_epsilon_b*m (m>=m_lim)
*/
  double res ;
  double m_lim = calcul_m_lim() ;
  if(m>=m_lim) {
     res  = f_epsilon_a*m_lim + f_epsilon_b*(m-m_lim) ;
  } else {
     res  = f_epsilon_a*m ;
  }
  if(m_n>=m_lim) {
     res  = res - (f_epsilon_a*m_lim + f_epsilon_b*(m_n-m_lim)) ;
  } else {
     res  = res - f_epsilon_a*m_n ;
  }


  return(res) ;
}



double calcul_delta_epsilon_pel(double sigm,double phi_1,double p_l1,double sigm_n,double p_l1n)
{
 
//pellet volumetric strain

  double sigm1_eff   = sigm/phi_1 + p_l1 ; // sigm < 0 compression ; s < 0
  double sigm1_eff_n = sigm_n/phi_1 + p_l1n ; // sigm < 0 compression ; s < 0

  double res = beta/alpha * ( exp(-alpha*(-sigm1_eff)) - exp(-alpha*(-sigm1_eff_n)) ) ;

  return(res) ;
}



double calcul_m_lim(void)
{
//m_lim = function of initial solid fraction

  double initial_solid_fraction_1 = p_pellet*rho_s/initial_pel_dry_dens ; 
  double res = A_m * initial_solid_fraction_1 + B_m ;
  
  return(res) ;
}



double calcul_dFdm(double m,double f_epsilon_a,double f_epsilon_b)
{
/* 
d_rond epsilon_V / d_rond m
*/
  double res ;
  double m_lim = calcul_m_lim() ;
  if(m>=m_lim) {
     res  = f_epsilon_b ;
  } else {
     res  = f_epsilon_a ;
  }
  return(res) ;
}
  
  

double calcul_deps1dp(double p_l1,double sigm,double phi_1)
{
/* 
d_rond epsV1 / d_rond p
*/
  double sigm1_eff = sigm/phi_1 + p_l1 ;
  double res = beta * exp(- alpha * ( - sigm1_eff ) ) / phi_1 ;
  
  return(res) ;
}



double calcul_deps1ds(double p_l1,double sigm,double phi_1)
{
/* 
d_rond epsV1 / d_rond s
s diminue, augmentation de volume
*/
  double sigm1_eff = sigm/phi_1 + p_l1 ;
  double res = beta * exp(- alpha * ( - sigm1_eff ) ) ;
  
  return(res) ;
}



double calcul_dmdp(double p_l1,double sigm,double phi_1)
{
/* 
d_rond m / d_rond p
*/
  
  double m_s = exp(-2./3.*alpha*(-p_l1)) ;
  double m_p = pow((-sigm),(2./3.))*exp(-2./3.*alpha*(-sigm)/phi_1);
  
  double coeff_m = beta / 3. * (1. - poisson*poisson) / (1. - 2. * poisson) ;
  coeff_m = pow(coeff_m , 2./3.) ;
  
  double res = - coeff_m * m_s * ( ( 2./3. * m_p / (-sigm) ) + ( -2./3. * alpha / phi_1 ) * m_p ) ;
  
  return(res) ;
}



double calcul_dmds(double p_l1,double sigm,double phi_1)
{
/* 
d_rond m / d_rond s
*/
  
  double m_s = exp(-2./3.*alpha*(-p_l1)) ;
  double m_p = pow((-sigm),(2./3.))*exp(-2./3.*alpha*(-sigm)/phi_1);
  
  double coeff_m = beta / 3. * (1. - poisson*poisson) / (1. - 2. * poisson) ;
  coeff_m = pow(coeff_m , 2./3.) ;
  
  double res = - (-2./3. * alpha) * coeff_m * m_p * m_s ;
  
  return(res) ;
}
    



double calcul_kappa_eq(double propPel,double sigm,double p_l1,double p_l2,double solid_fraction_1,double solid_fraction_2)
{
  double sigm1   = sigm / (solid_fraction_1 + solid_fraction_2) ;
  double sigm2   = sigm1 ;
          
  double sigmeff1 = sigm1 + p_l1 ; // p1
  double sigmeff2 = sigm2 + p_l2 ; // p2

  double res = - (propPel / sigmeff1 + (1-propPel) / sigmeff2) * kappa_m / (solid_fraction_1+solid_fraction_2) - kappa_M / sigm ;

  return (res) ;
}



double calcul_delta_epsilon_sigm(double sigm,double sigm_n,double p_l1n,double p_l2n,double solid_fraction_1,double solid_fraction_2)
{
  double propPel          = p_pellet / (p_pellet + p_poudre) ;
  double propPou          = p_poudre / (p_pellet + p_poudre) ;
  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange
  
  double sigm1_n = sigm_n / (solid_fraction_1 + solid_fraction_2) ; // p1_n
  double sigm2_n = sigm1_n ; // p2_n
             
  double sigmeff1_n = sigm1_n + p_l1n ; // p'1_n
  double sigmeff2_n = sigm2_n + p_l2n ; // p'2_n
  
  double delta_epsilon = - propPel * kappa_m / (1.+e0) * log( (sigm+(solid_fraction_1 + solid_fraction_2)*p_l1n) / ((solid_fraction_1 + solid_fraction_2)*sigmeff1_n ) ); // contribution de la deformation des pellets
  delta_epsilon       += - propPou * kappa_m / (1.+e0) * log( (sigm+(solid_fraction_1 + solid_fraction_2)*p_l2n) / ((solid_fraction_1 + solid_fraction_2)*sigmeff2_n ) ); // contribution de la déformation des grains de poudre
  delta_epsilon       += -           kappa_M / (1.+e0) * log( sigm / sigm_n ); // deformation macrostructure

  return(delta_epsilon) ;
}



double calcul_delta_epsilon_pel_cont(double sigm,double sigm_n,double p_l1,double p_l1n,double solid_fraction_1,double solid_fraction_2)
{
  double initial_e1 = rho_solid / initial_pel_dry_dens - 1. ; // indices des vides a partir de densite mesuree au labo
  
  double sigm1_n = sigm_n / (solid_fraction_1 + solid_fraction_2) ; // p1_n
  double sigmeff1_n = sigm1_n + p_l1n ; // p'1_n
  
  double delta_epsilon = - kappa_m / (1.+initial_e1) * log( (sigm+(solid_fraction_1 + solid_fraction_2)*p_l1n) / ((solid_fraction_1 + solid_fraction_2)*sigmeff1_n ) ); // sigm
  delta_epsilon       += - kappa_m / (1.+initial_e1) * log( (p_l1+sigm1_n) / (sigmeff1_n) ); // succion
  
  return(delta_epsilon) ;
}


double calcul_delta_epsilon_pou_cont(double sigm,double sigm_n,double p_l2,double p_l2n,double solid_fraction_1,double solid_fraction_2)
{
  double initial_e2 = rho_solid / initial_pel_dry_dens - 1. ; // indices des vides a partir de densite mesuree au labo
  
  double sigm2_n = sigm_n / (solid_fraction_1 + solid_fraction_2) ; // p2_n
  double sigmeff2_n = sigm2_n + p_l2n ; // p'2_n
  
  double delta_epsilon = - kappa_m / (1.+initial_e2) * log( (sigm+(solid_fraction_1 + solid_fraction_2)*p_l2n) / ((solid_fraction_1 + solid_fraction_2)*sigmeff2_n ) ); // sigm
  delta_epsilon       += - kappa_m / (1.+initial_e2) * log( (p_l2+sigm2_n) / (sigmeff2_n) ); // succion
  
  return(delta_epsilon) ;
}



double calcul_delta_epsilon_s1(double sigm_n,double p_l1,double p_l1n,double solid_fraction_1,double solid_fraction_2)
{
  double propPel          = p_pellet / (p_pellet + p_poudre) ;
  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange
  
  double sigm1_n = sigm_n / (solid_fraction_1 + solid_fraction_2) ; // p1_n
             
  double sigmeff1_n = sigm1_n + p_l1n ; // p'1_n
  
  double delta_epsilon = - propPel * kappa_m / (1.+e0) * log( (p_l1+sigm1_n) / (sigmeff1_n) ); // contribution de la deformation des pellets

  return(delta_epsilon) ;
}



double calcul_delta_epsilon_s2(double sigm_n,double p_l2,double p_l2n,double solid_fraction_1,double solid_fraction_2)
{
  double propPou          = p_poudre / (p_pellet + p_poudre) ;
  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange
  
  double sigm2_n = sigm_n / (solid_fraction_1 + solid_fraction_2) ; // p2_n
             
  double sigmeff2_n = sigm2_n + p_l2n ; // p'2_n
  
  double delta_epsilon = - propPou * kappa_m / (1.+e0) * log( (p_l2+sigm2_n) / (sigmeff2_n) ); // contribution de la déformation des grains de poudre
  delta_epsilon       += -           kappa_s / (1.+e0) * log(p_l2/p_l2n); // deformation macrostructure

  return(delta_epsilon) ;
}



double calcul_de1dp(double sigm,double p_l1,double solid_fraction_1,double solid_fraction_2)
{
  double propPel          = p_pellet / (p_pellet + p_poudre) ;
  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange
  
  double sigm1 = sigm / (solid_fraction_1 + solid_fraction_2) ; // p1
             
  double sigmeff1 = sigm1 + p_l1 ; // p'1
  
  double res = - propPel * kappa_m / (1.+ e0) / (solid_fraction_1 + solid_fraction_2) / sigmeff1 ;

  return(res) ;
}



double calcul_de2dp(double sigm,double p_l2,double solid_fraction_1,double solid_fraction_2)
{
  double propPou          = p_poudre / (p_pellet + p_poudre) ;
  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange
  
  double sigm2 = sigm / (solid_fraction_1 + solid_fraction_2) ; // p2
             
  double sigmeff2 = sigm2 + p_l2 ; // p'2

  double res = - propPou * kappa_m / (1.+ e0) / (solid_fraction_1 + solid_fraction_2) / sigmeff2 ;

  return(res) ;
}



double calcul_deMdp(double sigm)
{
//  // double initial_total_vr = rho_solid / rho_s -1. ; // indice des vides total : micro et macro, à partir de la masse vol. de la bentonite et de la masse vol seche du mélange

  double res = - kappa_M / (1.+ e0) / sigm ;

  return(res) ;
}



double bulkModulusGran(double dFdm, double dmdp, double deps1dp)
{
  double f_Mm = 1. ; // deformation variation volume pellets
  double res = dFdm * dmdp + f_Mm * deps1dp;
  res        = 1./res ;

  return(res) ;
}

double bulkModulusCont(double kappa_eq)
{

  double res = (1. + e0) / kappa_eq ;

  return(res) ;
}



double bulkModulus(double* sig, double* eps, double eps_1, double eps_2, double p_l1, double p_l2)
{
            double bulk = 0 ;

            double sigm   = (sig[0] + sig[4] + sig[8])/3. ;
            double tr_eps = eps[0] + eps[4] + eps[8] ; // deformation volumique mélange
            double initial_solid_fraction_1 = p_pellet*rho_s/initial_pel_dry_dens ; // compacités
            double initial_solid_fraction_2 = p_poudre*rho_s/initial_pel_dry_dens ;
            double solid_fraction_1         = initial_solid_fraction_1 * (1.+eps_1)/(1.+ tr_eps ) ;
            double solid_fraction_2         = initial_solid_fraction_2 * (1.+eps_2)/(1.+ tr_eps ) ;

            if(isGran == true) {
               double m       = calcul_m(sigm,solid_fraction_1,p_l1) ;
               double dFdm    = calcul_dFdm(m,B_a,B_b);
               double dmdp    = calcul_dmdp(p_l1,sigm,solid_fraction_1) ; // " d_rond m / d_rond p "
               double deps1dp = calcul_deps1dp(p_l1,sigm,solid_fraction_1) ; // "d_rond epsilon_V1 / d_rond p "
               bulk           = bulkModulusGran(dFdm, dmdp, deps1dp) ;
            } else {
               double propPel  = p_pellet / (p_pellet + p_poudre) ;
               double kappa_eq = calcul_kappa_eq(propPel,sigm,p_l1,p_l2,solid_fraction_1,solid_fraction_2) ;
               bulk            = bulkModulusCont(kappa_eq) ;
            }
   return (bulk) ;
}







