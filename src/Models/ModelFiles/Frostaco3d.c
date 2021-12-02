#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the finite element method */
#include "FEM.h"
#include "Elasticity.h"

#define TITLE "Frost actions in (3D) concrete"
#define AUTHORS "Tahiri-Dangla"

#include "PredefinedMethods.h"


/* Nb of equations of the model */
#define NEQ   (dim+3)

/* Nb of terms per point */
#define NVI   (21)    /*  nb of implicit terms per point */
#define NVE   (3)     /*  nb of explicit terms per point */
#define NV0   (0)     /*  nb of constant terms per point */


/* Indices of equations */
#define E_Mech   (0)
#define E_Mass   (dim)
#define E_Salt   (dim+1)
#define E_The    (dim+2)

/* Indices of unknowns (generic indices) */
#define U_Mech   E_Mech
#define U_Mass   E_Mass
#define U_Salt   E_Salt
#define U_The    E_The


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */
 
/* Mechanics */
#define U_dis    U_Mech

/* Mass */
//#define U_p_max  U_Mass
#define U_p_l    U_Mass

/* Salt */
#define U_c_s    U_Salt

/* Thermic */
#define U_tem    U_The




/* We define some names for implicit terms (vi must be used as pointer below) */
#define M_TOT         (vi)[0]
#define M_TOTn        (vi_n)[0]
#define W_TOT         (vi + 1) /* this is a 3D vector */

#define M_SALT        (vi + 4)[0]
#define M_SALTn       (vi_n + 4)[0]
#define W_SALT        (vi + 5)

#define S_TOT         (vi + 8)[0]
#define S_TOTn        (vi_n + 8)[0]
#define W_THE         (vi + 9)

#define SIG           (vi + 12) /* this a 3D tensor */



/* We define some names for explicit terms (ve must be used as pointer below) */
#define KD_LIQ        (ve)[0]
#define KF_SALT       (ve + 1)[0]
#define KTH           (ve + 2)[0]


/* We define some names for constant terms (v0 must be used as pointer below) */



/* Shorthands of some units */
#include "InternationalSystemOfUnits.h"

#define m        (InternationalSystemOfUnits_OneMeter)
#define m3       (m*m*m)
#define dm       (0.1*m)
#define cm       (0.01*m)
#define dm3      (dm*dm*dm)
#define cm3      (cm*cm*cm)
#define Pa       (InternationalSystemOfUnits_OnePascal)
#define MPa      (1.e6*Pa)
#define GPa      (1.e9*Pa)
#define Joule    (Pa*m3)
#define kg       (InternationalSystemOfUnits_OneKilogram)
#define mol      (InternationalSystemOfUnits_OneMole)
#define Kelvin   (InternationalSystemOfUnits_OneKelvin)
#define Watt     (InternationalSystemOfUnits_OneWatt)



/* Physical constants
 * ------------------ */
#include "PhysicalConstant.h"
/* Perfect gas constant (J/mol/K) */
static double R_g ;




/* Water properties
 * ---------------- */
 /* Molar mass (kg/mol) */
static double M_H2O ;
/* Molar volume of liquid water (m3/mol) */
static double V_H2O ;
/* Mass density of water */
static double rho_h2o_l0 ;
/* Molar volume of ice (m3/mol) */
static double V_Ice ;
/* Mass density of ice */
static double rho_h2o_i0 ;
/* Viscosity (Pa.s) */
#include "WaterViscosity.h"
/* Logarithm of water activity in brine */
#define LogActivityOfWater(c_s,T)      activity(c_s,T)
//#define LogActivityOfWater(c_s,T)      activity_w_ideal(c_s,T)
//#define LogActivityOfWater(c_s,T)      0
/* Entropy of fusion of ice (J/mol/K) */
static double S_m ;
/* Melting temperature of ice at atmospheric pressure */
static double T_m ;
/* Atmospheric pressure (Pa) */
static double p_m ;
/* Bulk modulus of liquid water (Pa) */
static double K_l ;
/* Bulk modulus of ice (Pa) */
static double K_i ;
/* Specific heat of liquid water (J/kg/K) */
static double C_l ;
/* Specific heat of ice (J/kg/K) */
static double C_i ;
/* Thermal conductivity of liquid water (W/m/K) */
#define LAM_l       ((0.6) * Watt / m / Kelvin)
static double lam_l ;
/* Thermal conductivity of ice (W/m/K) */
#define LAM_i       ((2.2) * Watt / m / Kelvin)
static double lam_i ;
/* Thermal volumetric expansion coefficient of liquid water (1/K) */
/* #define ALPHA_l    (-298.e-6) */
#define ALPHA_l(theta)    ((-68.7e-6 + 13.877e-6/Kelvin*(theta))/Kelvin)  /* After Speedy (1987) */
/* Thermal volumetric expansion coefficient of ice (1/K) */
#define ALPHA_i           ((155.e-6) * (1./Kelvin))
static double alpha_i ;



/* Salt properties
 * --------------- */
/* Types of salts: SALT = CAH (A=anion , C=cation , H=H2O) */
#define NaCl   0
#define CaCl2  1

/* Choose here the salt you want! */
#define SALT   NaCl

#include "MolarMassOfMolecule.h"
#include "PartialMolarVolumeOfMoleculeInWater.h"
#include "DiffusionCoefficientOfMoleculeInWater.h"
#if SALT == CaCl2
  /* Valencies */
  #define Z_A      (-1)
  #define Z_C      (+2)
  /* Stoichiometries */
  #define NU_A     2
  #define NU_C     1
  /* Partial molar volumes */
  #define V_A      PartialMolarVolumeOfMoleculeInWater(Cl)
  #define V_C      PartialMolarVolumeOfMoleculeInWater(Ca)
  /* Molar masses */
  #define M_ca     MolarMassOfMolecule(Ca)
  #define M_cl     MolarMassOfMolecule(Cl)
  #define M_Salt   (M_ca + 2*M_cl)
  /* Molecular diffusion coefficients (m2/s) */
  #define D_A      DiffusionCoefficientOfMoleculeInWater(Cl,T_m)
  #define D_C      DiffusionCoefficientOfMoleculeInWater(Ca,T_m)
#elif SALT == NaCl
  /* Valencies */
  #define Z_A      (-1)
  #define Z_C      (+1)
  /* Stoichiometries */
  #define NU_A     1
  #define NU_C     1
  /* Partial molar volumes */
  #define V_A      PartialMolarVolumeOfMoleculeInWater(Cl)
  #define V_C      PartialMolarVolumeOfMoleculeInWater(Na)
  /* Molar masses */
  #define M_na     MolarMassOfMolecule(Na)
  #define M_cl     MolarMassOfMolecule(Cl)
  #define M_Salt   (M_na + M_cl)
  /* Molecular diffusion coefficients (m2/s) */
  #define D_A      DiffusionCoefficientOfMoleculeInWater(Cl,T_m)
  #define D_C      DiffusionCoefficientOfMoleculeInWater(Na,T_m)
#else
  #error "Salt not available"
#endif
/* Molecular diffusion coefficient of the salt */
#define D_CA     ((Z_C - Z_A)*D_C*D_A/(Z_C*D_C - Z_A*D_A))
static double D_salt ;
/* Relative partial specific volume of the salt (m3/kg) */
static double vr_salt ;



/* Material Properties
 * ------------------- */
#define SaturationDegree(pc)                Curve_ComputeValue(Element_GetCurve(el),pc)
#define RelativePermeabilityToLiquid(pc)    Curve_ComputeValue(Element_GetCurve(el) + 1,pc)
#define Tortuosity(phi,sl)                  tortuosity(phi,sl)

static double  tortuosity(double,double) ;
double tortuosity(double p,double s)
{
  double tau_l_sat = 0.296e-3*exp(9.95*p)/p ;
  
  if(s > 0.) {
    return(tau_l_sat/(s*(1 + 625*pow(1 - s,4)))) ;
  } else {
    return(0.) ;
  }
}



/* The parameters below are read in the input data file */
static double phi ;
static double k_int ;
static double lam_s ;
static double C_s ;
static double k_s ;
static double g_s ;
static double alpha_s ;
static double p0 ;
static double T0 ;
static double* cijkl ;
static Elasticity_t* elasty ;



/* They are stored in the order specified in pm */
static int  pm(const char* s) ;
int pm(const char* s)
{
       if(strcmp(s,"porosity") == 0) return (0) ;
  else if(strcmp(s,"k_int") == 0)    return (1) ;
  else if(strcmp(s,"C_s") == 0)      return (2) ;
  else if(strcmp(s,"lam_s") == 0)    return (3) ;
  else if(strcmp(s,"k_s") == 0)      return (4) ;
  else if(strcmp(s,"g_s") == 0)      return (5) ;
  else if(strcmp(s,"alpha_s") == 0)  return (6) ;
  else if(strcmp(s,"p0") == 0)       return (7) ;
  else if(strcmp(s,"T0") == 0)       return (8) ;
  else return (-1) ;
}


/* They are retrieved automatically by calling the following function */
static void    GetProperties(Element_t*) ;
void GetProperties(Element_t* el)
{
/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
  phi     = GetProperty("porosity") ;
  k_int   = GetProperty("k_int") ;
  C_s     = GetProperty("C_s") ;
  lam_s   = GetProperty("lam_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;
  p0      = GetProperty("p0") ;
  T0      = GetProperty("T0") ;
  
  elasty  = Element_FindMaterialData(el,Elasticity_t,"Elasticity") ;
  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
#undef GetProperty
}



/* Functions used below */
static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static int     ComputeTransferCoefficients(Element_t*,double,double*) ;

static int     ComputeTangentCoefficients(Element_t*,double,double,double*) ;

static void    ComputePhysicoChemicalProperties(void) ;

static double  activity(double,double) ;
static double  activity_w_ideal(double,double) ;
static double  lna_i(double,double,double,double,double,double) ;



void ComputePhysicoChemicalProperties(void)
{
/* Physical constants
 * ------------------ */
/* Perfect gas constant (J/mol/K) */
  R_g = PhysicalConstant(PerfectGasConstant) ;


/* Water properties
 * ---------------- */
 /* Molar mass (kg/mol) */
  M_H2O = (18.e-3) * kg / mol ;
/* Molar volume of liquid water (m3/mol) */
  V_H2O = (18.e-6) * m3 / mol ;
/* Mass density of water */
  rho_h2o_l0 = M_H2O/V_H2O ;
/* Molar volume of ice (m3/mol) */
  V_Ice = (19.63e-6) * m3 / mol ;
/* Mass density of ice */
  rho_h2o_i0 = M_H2O/V_Ice ;
/* Entropy of fusion of ice (J/mol/K) */
  S_m = (23.54) * Joule / mol / Kelvin ;
/* Melting temperature of ice at atmospheric pressure */
  T_m = (273.) * Kelvin ;
/* Atmospheric pressure (Pa) */
  p_m = (1.e5) * Pa ;
/* Bulk modulus of liquid water (Pa) */
  K_l = (1.8e9) * Pa ;
/* Bulk modulus of ice (Pa) */
  K_i = (7.8e9) * Pa ;
/* Specific heat of liquid water (J/kg/K) */
  C_l = (4180.) * Joule / kg / Kelvin ;
/* Specific heat of ice (J/kg/K) */
  C_i = (2000.) * Joule / kg / Kelvin ;
/* Thermal conductivity of liquid water (W/m/K) */
  lam_l = LAM_l ;
/* Thermal conductivity of ice (W/m/K) */
  lam_i = LAM_i ;
/* Thermal volumetric expansion coefficient of ice (1/K) */
  alpha_i = ALPHA_i ;


/* Salt properties
 * --------------- */
  D_salt = D_CA ;
/* Relative partial specific volume of the salt (m3/kg) */
  vr_salt = ((NU_A*V_A + NU_C*V_C)/M_Salt - V_H2O/M_H2O) ;
}





/* We define some indices for the local variables */
enum {
I_DIS = 0,
I_U_Mass =  I_DIS + 3,
I_M_TOT,
I_M_SALT,
I_S_TOT,
I_S_L,
I_SD_L,
I_P_L,
I_P_I,
I_RHO_L,
I_RHO_I,
I_MU_L,
I_MU_I,
I_C_SALT,
I_TEM,
I_EPS,
I_SIG        = I_EPS        + 9,
I_W_TOT      = I_SIG        + 9,
I_W_SALT     = I_W_TOT      + 3,
I_W_THE      = I_W_SALT     + 3,
I_GRD_U_Mass = I_W_THE      + 3,
I_GRD_C_SALT = I_GRD_U_Mass + 3,
I_GRD_TEM    = I_GRD_C_SALT + 3,
I_KD_LIQ     = I_GRD_TEM    + 3,
I_KF_SALT,
I_KTH,
I_Last
} ;



/* Locally defined intern variables  */
#define NP               IntFct_MaxNbOfIntPoints
#define NbOfVariables    (I_Last)
static double Variable[NbOfVariables] ;
static double Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;





int SetModelProp(Model_t* model)
/** Set the model properties, return 0.
 *  Warning:
 *  Never call InternationalSystemOfUnits_UseAsLength() or similar
 *  to modify the units because this will also affect other models.
 */
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Salt,"salt") ;
  Model_CopyNameOfEquation(model,E_The,"the") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
#if defined (U_p_max)
  Model_CopyNameOfUnknown(model,U_Mass,"p_max") ;
#elif defined (U_p_l)
  Model_CopyNameOfUnknown(model,U_Mass,"p_l") ;
#endif
  Model_CopyNameOfUnknown(model,U_Salt,"c_s") ;
  Model_CopyNameOfUnknown(model,U_The,"tem") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_Mech + i,name_unk[i]) ;
  }
  
  Model_GetComputePropertyIndex(model) = pm ;
  
  Model_GetSequentialIndexOfUnknown(model)[E_The] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_Mass] = 1 ;
  Model_GetSequentialIndexOfUnknown(model)[E_Salt] = 1 ;
  for(i = 0 ; i < dim ; i++) {
    Model_GetSequentialIndexOfUnknown(model)[E_Mech+i] = 1 ;
  }

  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 9 ;

  /* By default */
  {
    int i ;
    
    for(i = 0 ; i < NbOfProp ; i++) {
      Material_GetProperty(mat)[i] = 0. ;
    }
  }
  
  ComputePhysicoChemicalProperties() ;
  
  /* Pre-initialization */
  {
    Material_GetProperty(mat)[pm("p0")]  = p_m ;
    Material_GetProperty(mat)[pm("T0")]  = T_m ;
  }
  
  Material_ScanProperties(mat,datafile,pm) ;
      
  /* Elasticity */
  {
    elasty = Elasticity_Create() ;
      
    Material_AppendData(mat,1,elasty,Elasticity_t,"Elasticity") ;
  }
  
  /* The 4th rank elastic tensor */
  {
    elasty = Material_FindData(mat,Elasticity_t,"Elasticity") ;
    
    k_s  = Material_GetPropertyValue(mat,"k_s") ;
    g_s  = Material_GetPropertyValue(mat,"g_s") ;
    phi  = Material_GetPropertyValue(mat,"porosity") ;

    /* isotropic Hooke's law */
    {    
      /* Elastic moduli (Mori-Tanaka) */
      double k_0 = k_s ;
      double g_0 = g_s ;
      double a_0 = 3 * k_0/(4 * g_0) ;
      double b_0 = 6 * (k_0 + 2 * g_0) / (9 * k_0 + 8 * g_0) ;
      double K = (1 - phi) * k_s / (1 + phi * a_0 * k_s / k_0) ;
      double G = (1 - phi) * g_s / (1 + phi * b_0 * g_s / g_0) ;
      double young =  9 * G * K / (3 * K + G) ;
      double poisson = 0.5 * (1 - young / (3 * K)) ;
    
      Elasticity_SetToIsotropy(elasty) ;
      Elasticity_SetParameters(elasty,young,poisson) ;
      
      {
        double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
        Elasticity_ComputeStiffnessTensor(elasty,c) ;
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




int PrintModelProp(Model_t* model,FILE* ficd)
/** Print the model properties 
 *  Return the nb of equations */
{ 
  int dim = Model_GetDimension(model) ;
  
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  
  printf("This model consists in %d equations:\n",NEQ) ;
  printf("\t 1. Conservation of mass   (mass)\n") ;
  printf("\t 2. Masse balance of salt  (salt)\n") ;
  printf("\t 3. Entropy balance        (the)\n") ;
  printf("\t 4. Equilibrium equation   (meca1,meca_2,meca_3)\n") ;
  
  
  printf("\n") ;
  
  printf("The primary unknowns are:\n") ;
#if defined (U_p_max)
  printf("\t 1. The max of liquid and ice pressures   (p_max)\n") ;
#elif defined (U_p_l)
  printf("\t 1. The liquid pressure     (p_l)\n") ;
#endif
  printf("\t 2. The salt concentration  (c_s)\n") ;
  printf("\t 3. The temperature         (tem)\n") ;
  printf("\t 4. The displacement vector (u_1,u_2,u_3)\n") ;


  printf("\n") ;
  
  printf("Example of input data\n\n") ;

  printf("porosity = 0.2   # Porosity \n") ;
  printf("k_int = 5.e-21   # Intrinsic permeability (m2) \n") ;
  printf("C_s = 2e+06      # Volumetric heat of solid (J/m3/K) \n") ;
  printf("lam_s = 1.       # Thermal conductivity of solid (W/m/K) \n") ;
  printf("k_s = 3.18e10    # Bulk modulus of the solid matrix (Pa)\n") ;
  printf("g_s = 1.91e10    # Shear modulus of the solid matrix (Pa) \n") ;
  printf("alpha_s = 54.e-6 # Thermal expansion of the solid matrix (1/K) \n") ;
  printf("Curves = myfile  # Name of file: p_c  S_l  k_rl\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;

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
#define R(n,i)    (r[(n)*NEQ+(i)])
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
    
    /* Thermic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_The) {
      
      for(i = 0 ; i < nn ; i++) {
        double tem = Element_GetValueOfNodalUnknown(el,u,i,U_The) ;
        
        R(i,E_The) /= tem ;
      }
    }
  }
  
  return(0) ;
#undef R
}



int ComputeInitialState(Element_t* el,double t)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* vi0 = Element_GetImplicitTerm(el) ;
  double* ve0 = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;


  /*
    We load some input data
  */
  GetProperties(el) ;


  /* Loop on integration points */
  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,vi0,t,0,p) ;
    
      /* storage in vi */
      {
        double* vi  = vi0 + p*NVI ;
        int    i ;
      
        M_TOT   = x[I_M_TOT] ;
        M_SALT  = x[I_M_SALT] ;
        S_TOT   = x[I_S_TOT] ;
    
        for(i = 0 ; i < 3 ; i++) W_TOT[i]  = x[I_W_TOT + i] ;
        for(i = 0 ; i < 3 ; i++) W_SALT[i] = x[I_W_SALT + i] ;
        for(i = 0 ; i < 3 ; i++) W_THE[i]  = x[I_W_THE + i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
      }
    
    
      /* storage in ve */
      {
        double p_l = x[I_P_L] ;
        double tem = x[I_TEM] ;
        double p_i = x[I_P_I] ;
        double p_c = p_i - p_l ;

        double sd_l = x[I_SD_L] ;
        double sd_i = 1 - sd_l ;
    
        /* Mass densities */
        double rho_l     = x[I_RHO_L] ;

        /* Darcy */
        double mu_l = WaterViscosity(tem) ;
        double kr_l = RelativePermeabilityToLiquid(p_c) ;
        double kd_liq = rho_l*k_int/mu_l*kr_l ;
    
        /* Fick */
        double tau_l   = Tortuosity(phi,sd_l) ;
        double kf_salt = M_Salt*phi*sd_l*tau_l*D_salt ;
    
        /* Fourier */
        double lam_h2o = sd_l*lam_l + sd_i*lam_i ;
        double lam_hom = lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
        double kth = lam_hom/tem ;
        
        double* ve  = ve0 + p*NVE ;
        
        KD_LIQ  = kd_liq ;
        KF_SALT = kf_salt ;
        KTH     = kth ;
      }
    }
  }

  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the (current) explicit terms.
 *  IMPORTANT: if needed use only the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  double* ve0 = Element_GetExplicitTerm(el) ;
  /* If you need the implicit terms, use the previous ones */
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  /* If you need the nodal values, use the previous ones */
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* If needed ! */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  
  
  /* Compute here ve */
  /* Loop on integration points */
  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,vi_n,t,0,p) ;
    
    
      /* storage in ve */
      {
        double p_l = x[I_P_L] ;
        double tem = x[I_TEM] ;
        double p_i = x[I_P_I] ;
        double p_c = p_i - p_l ;

        double sd_l = x[I_SD_L] ;
        double sd_i = 1 - sd_l ;
    
        /* Mass densities */
        double rho_l     = x[I_RHO_L] ;

        /* Darcy */
        double mu_l = WaterViscosity(tem) ;
        double kr_l = RelativePermeabilityToLiquid(p_c) ;
        double kd_liq = rho_l*k_int/mu_l*kr_l ;
    
        /* Fick */
        double tau_l   = Tortuosity(phi,sd_l) ;
        double kf_salt = M_Salt*phi*sd_l*tau_l*D_salt ;
    
        /* Fourier */
        double lam_h2o = sd_l*lam_l + sd_i*lam_i ;
        double lam_hom = lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
        double kth = lam_hom/tem ;
        
        double* ve  = ve0 + p*NVE ;
        
        KD_LIQ  = kd_liq ;
        KF_SALT = kf_salt ;
        KTH     = kth ;
      }
    }
  }
  
  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the (current) implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi0  = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  
  /* Compute here vi (with the help of vi_n if needed) */
  /* Loop on integration points */
  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vi_n,t,dt,p) ;
    
      /* storage in vi */
      {
        double* vi  = vi0 + p*NVI ;
        int    i ;
      
        M_TOT   = x[I_M_TOT] ;
        M_SALT  = x[I_M_SALT] ;
        S_TOT   = x[I_S_TOT] ;
    
        for(i = 0 ; i < 3 ; i++) W_TOT[i]  = x[I_W_TOT + i] ;
        for(i = 0 ; i < 3 ; i++) W_SALT[i] = x[I_W_SALT + i] ;
        for(i = 0 ; i < 3 ; i++) W_THE[i]  = x[I_W_THE + i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
      }
    }
  }

  return(0) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;

  /* Initialization */
  {
    int    i ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;
  }


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    We load some input data
  */
  GetProperties(el) ;


  /*
  ** Poromechanics matrix
  */
  {
    int n = 81 + 3*9 + 3*(9 + 3) ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTangentCoefficients(el,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,3) ;
    /* The matrix kp is stored as 
     * (u: displacement, p: pressure, c: concentration, t: temperature)
     * |Kuu  Kup  Kuc  Kut|
     * |Kpu  Kpp  Kpc  Kpt|
     * |Kcu  Kcp  Kcc  Kct|
     * |Ktu  Ktp  Ktc  Ktt|
     * i.e. the displacements u are in the positions 0 to dim-1,
     * the pressure p is in the position dim, c in dim+1 and t in dim+2.
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
    int n = 9*9 ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTransferCoefficients(el,dt,c) ;
    double* kc_mass = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    double* kc_salt = FEM_ComputeConductionMatrix(fem,intfct,c+9*4,dec) ;
    double* kc_the  = FEM_ComputeConductionMatrix(fem,intfct,c+9*8,dec) ;
    int    i ;
  
    for(i = 0 ; i < nn ; i++) {
      int    j ;
      
      for(j = 0 ; j < nn ; j++) {
        K(E_Mass + i*NEQ,U_Mass + j*NEQ) += dt*kc_mass[i*nn + j] ;
        K(E_Salt + i*NEQ,U_Salt + j*NEQ) += dt*kc_salt[i*nn + j] ;
        K(E_The  + i*NEQ,U_The  + j*NEQ) += dt*kc_the[i*nn + j] ;
      }
    }
  }
  
  return(0) ;
#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
  double* vi_1 = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }

  if(Element_IsSubmanifold(el)) return(0) ;
  

  /* Compute here the residu R(n,i) */
  

  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  if(Element_EquationIsActive(el,E_Mech)) 
  {
    double* vi = vi_1 ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  #if 0
  if(Element_EquationIsActive(el,E_Mech)) 
  {
    double* vi = vi_1 ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Mech + dim - 1) -= -rbf[i] ;
    }
    
  }
  #endif
  
  
  
  /* 2. Conservation of total mass */
  
  /* 2.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_Mass)) 
  {
    double* vi = vi_1 ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI , vi_n += NVI) g1[i] = M_TOT - M_TOTn ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Mass) -= ra[i] ;
    }
  }
  
  /* 2.2 Transport Terms */
  if(Element_EquationIsActive(el,E_Mass)) 
  {
    double* vi = vi_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_TOT,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Mass) -= -dt*rf[i] ;
  }
  
  
  
  /* 3. Conservation of salt */
  
  /* 3.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_Salt)) 
  {
    double* vi = vi_1 ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI , vi_n += NVI) g1[i] = M_SALT - M_SALTn ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Salt) -= ra[i] ;
    }
  }
  
  /* 3.2 Transport Terms */
  if(Element_EquationIsActive(el,E_Salt)) 
  {
    double* vi = vi_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_SALT,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Salt) -= -dt*rf[i] ;
  }
  
  
  
  /* 4. Entropy balance */
  
  /* 3.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_The)) 
  {
    double* vi = vi_1 ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI , vi_n += NVI) g1[i] = S_TOT - S_TOTn ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_The) -= ra[i] ;
    }
  }
  
  /* 3.2 Transport Terms */
  if(Element_EquationIsActive(el,E_The)) 
  {
    double* vi = vi_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_THE,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_The) -= -dt*rf[i] ;
  }
  
  return(0) ;
#undef R
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  int NbOfOutputs  = 17 ;
  double* ve  = Element_GetExplicitTerm(el) ;
  double* vi  = Element_GetCurrentImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  //int dim = Element_GetDimensionOfSpace(el) ;
  //FEM_t* fem = FEM_GetInstance(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /* Initialization */
  {
    int i ;
    
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
    double* x = ComputeVariables(el,u,u,vi,t,0,p) ;
    /* strains */
    double* eps = x + I_EPS ;
    double  eps_v = eps[0] + eps[4] + eps[8] ;
    /* pressures */
    double dp_l = x[I_P_L] - p0 ;
    double dp_i = x[I_P_I] - p0 ;
    /* temperature */
    double dtem = x[I_TEM] - T0 ;
    /* fluxes */
    double w_tot[3]  = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double w_the[3]  = {0,0,0} ;
    /* Poroelastic coefficients: homogenization
     * Mori-Tanaka:     k_0 = k_s and g_0 = g_s
     * Self-consistent: k_0 = K   and g_0 = G (for the record) */
    double k_0 = k_s ;
    double g_0 = g_s ;
    double a_0 = 3 * k_0/(4 * g_0) ;
    double b_0 = 6 * (k_0 + 2 * g_0) / (9 * k_0 + 8 * g_0) ;
    double K = (1 - phi) * k_s / (1 + phi * a_0 * k_s / k_0) ;
    double G = (1 - phi) * g_s / (1 + phi * b_0 * g_s / g_0) ;
    double b    = 1 - K/k_s ;
    double sd_l = x[I_SD_L] ;
    double sd_i = 1 - sd_l ;
    double b_i  = b*sd_i ;
    double b_l  = b*sd_l ;
    double N    = (b - phi)/k_s ;
    double N_il = sd_i*sd_l*(N - 0.75*phi/g_s) ;
    double N_ii = sd_i*N - N_il ;
    double N_ll = sd_l*N - N_il ;
    double alpha_phi = (b - phi)*alpha_s ;
    double alpha_phi_l = sd_l*alpha_phi ;
    double alpha_phi_i = sd_i*alpha_phi ;
    
    double k_oedo = K + 4/3.*G ; /* Oedometric modulus */
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    
    int    i ;
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vi += NVI , ve += NVE) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_tot[j]  += W_TOT[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_salt[j] += W_SALT[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_the[j]  += W_THE[j]/np ;
      for(j = 0 ; j < 9 ; j++) sig[j]    += SIG[j]/np ;
    }
    
    i = 0  ;
    /* quantites exploitees */
    Result_Store(r + i++,x + I_P_L,"Liquid pressure",1) ;
    Result_Store(r + i++,x + I_C_SALT,"Salt concentration",1) ;
    Result_Store(r + i++,x + I_TEM,"Temperature",1) ;
    Result_Store(r + i++,x + I_SD_L,"Liquid saturation",1) ;
    Result_Store(r + i++,w_tot,"Liquid mass flux vector",3) ;
    Result_Store(r + i++,w_salt,"Salt mass flux vector",3) ;
    Result_Store(r + i++,x + I_P_I,"Ice pressure",1) ;
    Result_Store(r + i++,x + I_RHO_L,"Solution density",1) ;
    Result_Store(r + i++,x + I_RHO_I,"Ice density",1) ;
    Result_Store(r + i++,x + I_EPS,"Strain tensor",9) ;
    Result_Store(r + i++,sig,"Stress tensor",9) ;
    Result_Store(r + i++,x + I_DIS,"Displacement vector",3) ;
      
    {
      double Eps_T  = alpha_s*K*dtem/k_oedo ;
      
      Result_Store(r + i++,&Eps_T,"Strain_temperature",1) ;
    }
    
    {
      double Eps_P  = (b_l*dp_l + b_i*dp_i)/k_oedo ;
      
      Result_Store(r + i++,&Eps_P,"Strain_pressure",1) ;
    }
    
    {
      double phi_i  = b_i*eps_v + N_il*dp_l + N_ii*dp_i - alpha_phi_i*dtem ;
      
      Result_Store(r + i++,&phi_i,"Deformation of ice pores",1) ;
    }
    
    {
      double phi_l  = b_l*eps_v + N_ll*dp_l + N_il*dp_i - alpha_phi_l*dtem ;
    
      Result_Store(r + i++,&phi_l,"Deformation of liquid pores",1) ;
    }
    
    {
      double lnaw = LogActivityOfWater(x[I_C_SALT],x[I_TEM]) ;
      
      Result_Store(r + i++,&lnaw,"Lna_w",1) ;
    }
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(Element_t* el,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define T2(a,i,j)      ((a)[(i)*3+(j)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double dui[NEQ] ;
  int    dec = 81 + 3*9 + 3*(9 + 3) ;
  
  
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
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vi_n,t,dt,p) ;
      double* c0 = c + p*dec ;


      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
      
      
      /* The derivative of equations w.r.t unknowns */
      

      /* Derivatives w.r.t strains */
      {
        /* Tangent stiffness matrix */
        {
          double* c1 = c0 ;
          int i ;
          
          for(i = 0 ; i < 81 ; i++) {
            c1[i] = cijkl[i] ;
          }

          //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
        }
        
        /* Coupling matrices */
        {
          double deps = 1.e-6 ;
          double* dx = ComputeVariableDerivatives(el,t,dt,x,deps,I_EPS) ;
          double* c00 = c0 + 81 + 3*9 ;
          
          /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
           * and zero for the derivatives w.r.t others */
          {
            double* c1 = c00 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 12 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_M_SALT] ;
          }
          {
            double* c1 = c00 + 2*12 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_S_TOT] ;
          }
        }
      }
      
      /* Derivatives w.r.t U_Mass */
      {
        double  dp = dui[U_Mass] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_Mass) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdp = dx + I_SIG ;
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dx[I_M_SALT] ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dx[I_S_TOT] ;
          }
        }
      }
      
      /* Derivatives w.r.t U_Salt */
      {
        double  dc = dui[U_Salt] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dc,I_C_SALT) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdc = dx + I_SIG ;
          double* c1 = c0 + 81 + 9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdc[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 + 1 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dx[I_M_SALT] ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dx[I_S_TOT] ;
          }
        }
      }
      
      /* Derivatives w.r.t  U_The */
      {
        double  dtem = dui[U_The] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dtem,I_TEM) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdtem = dx + I_SIG ;
          double* c1 = c0 + 81 + 2*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdtem[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 + 2 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dx[I_M_SALT] ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dx[I_S_TOT] ;
          }
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



int ComputeTransferCoefficients(Element_t* el,double dt,double* c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  double* ve0 = Element_GetExplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 * 9 ;
  int    p ;
  

  for(p = 0 ; p < np ; p++) {
    int i ;
    double* c1 = c + p*dec ;
    
    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
    
    {
      double* ve  = ve0 + p*NVE ;
      
      /* Permeability tensor */
      {
        double* c2 = c1 ;
        
        c2[0] = KD_LIQ ;
        c2[4] = KD_LIQ ;
        c2[8] = KD_LIQ ;
      }
      
      /* Fick tensor */
      {
        double* c2 = c1 + 4 * 9 ;
        
        c2[0] = KF_SALT ;
        c2[4] = KF_SALT ;
        c2[8] = KF_SALT ;
      }
      
      /* Thermal conductivity tensor */
      {
        double* c2 = c1 + 8 * 9 ;
        
        c2[0] = KTH ;
        c2[4] = KTH ;
        c2[8] = KTH ;
      }
    }
  }

  return(dec) ;
}








double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
/** This locally defined function compute the intern variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables. */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ; 
  double*  x   = Variable ;
  double*  x_n = Variable_n ;
  /* Variables is a locally defined array of array */
  //Model_t*  model  = Element_GetModel(el) ;
  //double*   x      = Model_GetVariable(model,p) ;
  
    
  /* Load the primary variables in x */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_DIS + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_DIS + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Pressure */
    x[I_U_Mass] = FEM_ComputeUnknown(fem,u,intfct,p,U_Mass) ;
    
    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Mass) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_Mass + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    
    /* Salt concentration */
    x[I_C_SALT] = FEM_ComputeUnknown(fem,u,intfct,p,U_Salt) ;
    
    /* Salt concentration gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Salt) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_C_SALT + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    
    /* Temperature */
    x[I_TEM] = FEM_ComputeUnknown(fem,u,intfct,p,U_The) ;
    
    /* Temperature gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_The) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_TEM + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    /* Transfer coefficient */
    {
      double* ve0 = Element_GetExplicitTerm(el) ;
      double* ve  = ve0 + p*NVE ;
      
      x[I_KD_LIQ]  = KD_LIQ ;
      x[I_KF_SALT] = KF_SALT ;
      x[I_KTH]     = KTH ;
    }

  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,const double t,const double dt,double* x_n,double* x)
/** Compute the secondary variables from the primary ones. */
{
#if defined (U_p_max)
  int dim = Element_GetDimensionOfSpace(el) ;
#endif
  /* Retrieve the primary variables from x */
  /* Strains */
  double* eps   = x + I_EPS ;
  /* Salt concentration */
  double c_s    = x[I_C_SALT] ;
  /* Temperature */
  double tem    = x[I_TEM] ;
  /* Log of water activity */
  double lna = LogActivityOfWater(c_s,tem) ;
  
  /* Pressures */
#if defined (U_p_max)
  double  p_max = x[I_U_Mass] ;
  //     0 = V_Ice*(p_i   - p_m) - V_H2O*(p_l   - p_m) + S_m*(tem - T_m) - R_g*tem*lna
  double c = V_Ice*(p_max - p_m) - V_H2O*(p_max - p_m) + S_m*(tem - T_m) - R_g*tem*lna ;
  double p_min = (c > 0) ? p_max - c/V_Ice : p_max + c/V_H2O ;
  double p_i   = (c > 0) ? p_min : p_max ;
  double p_l   = (c > 0) ? p_max : p_min ;
#elif defined (U_p_l)
  double p_l   = x[I_U_Mass] ;
  double p_i   = p_m + (V_H2O*(p_l - p_m) - S_m*(tem - T_m) + R_g*tem*lna)/V_Ice ;
#endif
  
  double p_c = p_i - p_l ;
  

  /* Saturations */
  double sd_l = SaturationDegree(p_c) ;
  double sd_i = 1 - sd_l ;
  

  /* Mass densities */
  double alpha_l   = ALPHA_l(tem - T_m) ;
  double rho_i     = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
  double rho_salt  = M_Salt*c_s ;
  double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - vr_salt*rho_salt) ;
  
  /* Mass contents */
  double m_l     = phi*sd_l*rho_l ;
  double m_i     = phi*sd_i*rho_i ;
  double m_salt  = phi*sd_l*rho_salt ;
  
  /* Entropies */
  double s_l   = C_l*log(tem/T_m) - alpha_l/rho_h2o_l0*(p_l - p_m) ;
  double s_i   = C_i*log(tem/T_m) - S_m/M_H2O - alpha_i/rho_h2o_i0*(p_i - p_m) ;
  double s_sol = C_s*log(tem/T_m) ;
  double s_tot = s_sol + m_l*s_l + m_i*s_i ;
    
  
  /* Backup */
  
  x[I_S_TOT   ] = s_tot ;
  x[I_S_L     ] = s_l ;
    
  x[I_M_TOT   ] = m_l + m_i ;
  x[I_M_SALT  ] = m_salt ;
    
  x[I_RHO_L   ] = rho_l ;
  x[I_RHO_I   ] = rho_i ;
    
  x[I_SD_L    ] = sd_l ;

  x[I_P_L     ] = p_l ;
  x[I_P_I     ] = p_i ;
  
  /* Chemical potentials */
  {
    double mu_l = V_H2O*(p_l - p_m) + R_g*tem*lna ;
    double mu_i = V_Ice*(p_i - p_m) + S_m*(tem - T_m) ;
    
    x[I_MU_L    ] = mu_l ;
    x[I_MU_I    ] = mu_i ;
  }
  
  /* Stresses */
  {
    double K    = (1 - phi)*k_s*g_s/(0.75*phi*k_s + g_s) ;
    double b    = 1 - K/k_s ;
    double b_i  = b*sd_i ;
    double b_l  = b*sd_l ;
    
    double dp_l = p_l - p0 ;
    double dp_i = p_i - p0 ;
    double dtem = tem - T0 ;
    
    double* sig = x + I_SIG ;
    
    {
      int    i ;
    
      /* Initialize stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = 0 ;
      
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        for(j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*eps[j] ;
        }
      }
      #undef C
      
      sig[0] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
      sig[4] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
      sig[8] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
    }
  }
  
  /* Flux */
  {
    double* w_tot   = x + I_W_TOT ;
    double* w_salt  = x + I_W_SALT ;
    double* w_the   = x + I_W_THE ; 
    double  kd_liq  = x[I_KD_LIQ] ;
    double  kf_salt = x[I_KF_SALT] ;
    double  kth     = x[I_KTH] ;
    
    double* grd_c_s   = x + I_GRD_C_SALT ;
    double* grd_tem   = x + I_GRD_TEM ;
    int i ;
    
#if defined (U_p_max)
    double  grd_p_l[3] ;
    double* grd_p_max = x + I_GRD_U_Mass ;  
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    double dc_s = 1.e-2 * ObVal_GetValue(obval + U_Salt) ;
    double dtlnadc_s = (tem * LogActivityOfWater(c_s + dc_s,tem) - tem * lna) / dc_s ;
    double dtem = 1.e-2 * ObVal_GetValue(obval + U_The) ;
    double dtlnadtem = ((tem + dtem) * LogActivityOfWater(c_s,tem + dtem) - tem * lna) / dtem ;
    
    for(i = 0 ; i < 3 ; i++) {
      double grd_tlna = dtlnadc_s * grd_c_s[i] + dtlnadtem * grd_tem[i] ;
      double grd_c    = (V_Ice - V_H2O) * grd_p_max[i] + S_m * grd_tem[i] - R_g * grd_tlna ;
      double grd_p_min = (c > 0) ? grd_p_max[i] - grd_c/V_Ice : grd_p_max[i] + grd_c/V_H2O ;
      
      grd_p_l[i]   = (c > 0) ? grd_p_max[i] : grd_p_min ;
    }
#elif defined (U_p_l)
    double* grd_p_l = x + I_GRD_U_Mass ;  
#endif
    
    
    for(i = 0 ; i < 3 ; i++) {
      double w_l     = - kd_liq * grd_p_l[i] ;
      double j_salt  = - kf_salt * grd_c_s[i] ;
      double kc_salt =   M_Salt * c_s / rho_l ;
      double q       = - kth * grd_tem[i] ;
  
      w_tot[i]   = w_l ;
      w_salt[i]  = kc_salt * w_l + j_salt ;
      w_the[i]   = q + s_l * w_l ;
    }
  }
  
}


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








double activity(double c_s1,double tem1)
/* activity of water */
{
  double c_s = (c_s1 > 0) ? c_s1 : 0 ;
  double tem = (tem1 > 200) ? tem1 : 200 ;
  double T_0  = T_m ;
  /* References */
  double b0   = sqrt(M_H2O) ;
  double S0   = pow(M_H2O,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_na_nacl = 4.352/b0 ;
  double b_cl_nacl = 1.827/b0 ;
  double S_na_nacl = 26.448/S0 ;
  double S_cl_nacl = 19.245/S0 ;
  /* CaCl2 (d'apres Lin et Lee) */
  double b_ca_cacl2 = 3.908/b0 ;
  double b_cl_cacl2 = 2.085/b0 ;
  double S_ca_cacl2 = 18.321/S0 ;
  double S_cl_cacl2 = 10.745/S0 ;

  double epsi = 0.0007*(tem - T_0)*(tem - T_0) - 0.3918*(tem - T_0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*tem,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;
  double c_ani ;
  double c_cat ;
  double z_ani ;
  double z_cat ;

  switch(SALT) {
    case(NaCl) : {
      b_cat = b_na_nacl ;
      b_ani = b_cl_nacl ;
      S_cat = S_na_nacl ;
      S_ani = S_cl_nacl ;
      c_ani = c_s ;
      c_cat = c_s ;
      z_ani = -1 ;
      z_cat = 1 ;
      break ;
    }
    case(CaCl2) : {
      b_cat = b_ca_cacl2 ;
      b_ani = b_cl_cacl2 ;
      S_cat = S_ca_cacl2 ;
      S_ani = S_cl_cacl2 ;
      c_ani = 2*c_s ;
      c_cat = c_s ;
      z_ani = -1 ;
      z_cat = 2 ;
      break ;
    }
    default : {
      arret("non prevu") ;
    }
  }
  
  {
    /* concentrations */
    double c_h2o = (1 - c_s*(NU_A*V_A + NU_C*V_C))/V_H2O ;
    /* molalites * M_H2O */
    double m_ani  = c_ani/c_h2o ;
    double m_cat  = c_cat/c_h2o ;

    /* ionic strength */
    double I     =  0.5*(z_ani*z_ani*m_ani + z_cat*z_cat*m_cat);
  
    double II_ani   = lna_i(tem,I,z_ani,b_ani,S_ani,A) ;
    double II_cat   = lna_i(tem,I,z_cat,b_cat,S_cat,A) ;

    /* activity of water */
    double lna_h2o = m_ani*II_ani + m_cat*II_cat ;
    /* linearized activity of water */
    //double lna_h2o = - (m_ani + m_cat) ;

    return(lna_h2o) ;
  }
}


double lna_i(double T,double I,double z,double b,double S,double A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29 ;
  double a1 = alpha/(1+alpha) ;
  double II = sqrt(I) ;
  double lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/T ;
  double lna_i = -1 + lna*z*z ;
  
  return(lna_i) ;
}


double activity_w_ideal(double c_s,double tem)
/* Natural log of water activity in ideal mixture */
{
  double c_cat = NU_C*c_s ;
  double c_ani = NU_A*c_s ;
  //double c_w   = (1. - (NU_A*V_A + NU_C*V_C)*c_s)/V_H2O ;
  //double c_t   = c_w + c_cat + c_ani ;
  //double lna_w = log(c_w/c_t) ;
  /* linearized term */
  double lna_w = - (c_cat + c_ani)*V_H2O ;
  return(lna_w) ;
}
