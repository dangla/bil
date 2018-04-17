#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the numerical method */
#include "FVM.h"

#define TITLE   "Frost actions in concrete"
#define AUTHORS "Zeng,Fen-Chong,Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ      (3)
//#define NEQ      (4)


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Nb of terms */
#define NVI      (3*NN*NN + 2*NN)
#define NVE      (5*NN)
#define NV0      (0)

/* Indices of equations */
#define E_Mass   (0)
#define E_Salt   (1)
#define E_The    (2)
//#define E_ChemP  (3)

#if defined (E_ChemP) && (NEQ != 4)
  #error "Check the nb of equation"
#endif


/* Indices of unknowns */
//#define U_p_l    (0)
//#define U_p_i    (0)
#define U_p_max  (0)
#define U_c_s    (1)
#define U_tem    (2)

#if defined (U_p_l) && defined (U_p_i) && (NEQ == 3)
  #error "One pressure exactly must be defined"
#endif


/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Names of nodal unknowns */
#if defined (U_p_l)
  #define P_l(n)   (UNKNOWN(n,U_p_l))
  #define P_ln(n)  (UNKNOWNn(n,U_p_l))
#endif

#if defined (U_p_i)
  #define P_i(n)   (UNKNOWN(n,U_p_i))
  #define P_in(n)  (UNKNOWNn(n,U_p_i))
#endif

#if defined (U_p_max)
  #define P_max(n)   (UNKNOWN(n,U_p_max))
  #define P_maxn(n)  (UNKNOWNn(n,U_p_max))
#endif

#define C_s(n)   (UNKNOWN(n,U_c_s))
#define C_sn(n)  (UNKNOWNn(n,U_c_s))

#define TEM(n)   (UNKNOWN(n,U_tem))
#define TEM_n(n) (UNKNOWNn(n,U_tem))


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define M_TOT(i)    MassAndFlux(f,i,i)
#define M_TOTn(i)   MassAndFlux(f_n,i,i)
#define W_TOT(i,j)  MassAndFlux(f,i,j)

#define M_SALT(i)   MassAndFlux(f + NN*NN,i,i)
#define M_SALTn(i)  MassAndFlux(f_n + NN*NN,i,i)
#define W_SALT(i,j) MassAndFlux(f + NN*NN,i,j)

#define S_TOT(i)    MassAndFlux(f + 2*NN*NN,i,i)
#define S_TOTn(i)   MassAndFlux(f_n + 2*NN*NN,i,i)
#define W_THE(i,j)  MassAndFlux(f + 2*NN*NN,i,j)

#define Mu_Liq(i)   (f + 3*NN*NN)[i]
#define Mu_Ice(i)   (f + 3*NN*NN + NN)[i]




/* Names used for explicit terms */
#define KD_LIQ      (va)
#define KF_SALT     (va + NN)
#define KC_SALT     (va + 2*NN)
#define KTH         (va + 3*NN)
#define S_L         (va + 4*NN)



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
//#include "Log10ActivityOfWaterInBrine.h"
#define LogActivityOfWater(c_s,T)      activity(c_s,T_m)
//#define LogActivityOfWater(c_s,T)      activite_w_ideal(c_s,T)
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
/* Relative specific partial volume of the salt (m3/kg) */
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
static double phi,k_int,lam_s,C_s,k_s,g_s,alpha_s,p0,T0 ;

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
#undef GetProperty
}



/* Functions used below */
static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;

static int     TangentCoefficients(Element_t*,double,double,double*) ;

static void    ComputePhysicoChemicalProperties(void) ;

static double  activity(double,double) ;
static double  activite_w_ideal(double,double) ;
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
/* Relative specific partial volume of the salt (m3/kg) */
  vr_salt = ((NU_A*V_A + NU_C*V_C)/M_Salt - V_H2O/M_H2O) ;
}



#define NbOfVariables    (NEQ+11)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;

#define I_M_TOT      (NEQ+0)
#define I_M_SALT     (NEQ+1)
#define I_S_TOT      (NEQ+2)
#define I_S_L        (NEQ+3)
#define I_SD_L       (NEQ+4)
#define I_P_L        (NEQ+5)
#define I_P_I        (NEQ+6)
#define I_RHO_L      (NEQ+7)
#define I_RHO_I      (NEQ+8)
#define I_MU_L       (NEQ+9)
#define I_MU_I       (NEQ+10)
  
  

#define NbOfVariableFluxes    (3)
static double VariableFluxes[NbOfVariableFluxes] ;

#define I_W_TOT           (0)
#define I_W_SALT          (1)
#define I_W_THE           (2)




int SetModelProp(Model_t* model)
/** Set the model properties 
 *  Return 0 */
{
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Salt,"salt") ;
  Model_CopyNameOfEquation(model,E_The,"the") ;
#if defined (E_ChemP)
  Model_CopyNameOfEquation(model,E_ChemP,"chemp") ;
#endif
  
  /** Names of the main (nodal) unknowns */
#if defined (U_p_l)
  Model_CopyNameOfUnknown(model,U_p_l,"p_l") ;
#endif
#if defined (U_p_i)
  Model_CopyNameOfUnknown(model,U_p_i,"p_i") ;
#endif
#if defined (U_p_max)
  Model_CopyNameOfUnknown(model,U_p_max,"p_max") ;
#endif
  Model_CopyNameOfUnknown(model,U_c_s,"c_s") ;
  Model_CopyNameOfUnknown(model,U_tem,"tem") ;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 9 ;
  
  ComputePhysicoChemicalProperties() ;
  
  /* Pre-initialization */
  {
    Material_GetProperty(mat)[pm("p0")]  = p_m ;
    Material_GetProperty(mat)[pm("T0")]  = T_m ;
  }
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}




int PrintModelProp(Model_t* model,FILE* ficd)
/** Print the model properties 
 *  Return the nb of equations */
{ 
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  
  printf("This model consists in %d equations:\n",NEQ) ;
  printf("\t 1. Conservation of mass   (mass)\n") ;
  printf("\t 2. Masse balance of salt  (salt)\n") ;
  printf("\t 3. Entropy balance        (the)\n") ;
  
  
  printf("\n") ;
  
  printf("The primary unknowns are:\n") ;
#if defined (U_p_l)
  printf("\t 1. The liquid pressure    (p_l)\n") ;
#endif
#if defined (U_p_i)
  printf("\t 1. The ice pressure       (p_i)\n") ;
#endif
#if defined (U_p_max)
  printf("\t 1. The max of liquid and ice pressures   (p_max)\n") ;
#endif
  printf("\t 2. The salt concentration (c_s)\n") ;
  printf("\t 3. The temperature        (tem)\n") ;


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

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    i ;
  

  {
    double* r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
    
    /* Thermic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_The) {
      R(0,E_The) /= TEM_n(0) ;
    }
  }
  
  return(0) ;
#undef R
}



int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* f  = Element_GetImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  
  
  /* We can skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Input data */
  GetProperties(el) ;
  
  /* Masses and entropies */
  {
    int    i ;
  
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,f,0,0,i) ;

      M_TOT(i)   = x[I_M_TOT] ;
      M_SALT(i)  = x[I_M_SALT] ;
      S_TOT(i)   = x[I_S_TOT] ;
      Mu_Liq(i)  = x[I_MU_L] ;
      Mu_Ice(i)  = x[I_MU_I] ;
    }
  }

  {
    ComputeTransferCoefficients(el,u,f) ;
  }

  /* Flux */
  {
    int    i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;

        W_TOT(i,j)  =   w[I_W_TOT] ;
        W_SALT(i,j) =   w[I_W_SALT] ;
        W_THE(i,j)  =   w[I_W_THE] ;
        
        W_TOT(j,i)  = - w[I_W_TOT] ;
        W_SALT(j,i) = - w[I_W_SALT] ;
        W_THE(j,i)  = - w[I_W_THE] ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t* el,double t)
/* Explicit terms (va)  */
{
  double* f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Transfer coefficients
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}




int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  

  
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  /* Masses and entropies */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;

      M_TOT(i)   = x[I_M_TOT] ;
      M_SALT(i)  = x[I_M_SALT] ;
      S_TOT(i)   = x[I_S_TOT] ;
      Mu_Liq(i)  = x[I_MU_L] ;
      Mu_Ice(i)  = x[I_MU_I] ;
    }
  }

  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;

        W_TOT(i,j)  =   w[I_W_TOT] ;
        W_SALT(i,j) =   w[I_W_SALT] ;
        W_THE(i,j)  =   w[I_W_THE] ;
        
        W_TOT(j,i)  = - w[I_W_TOT] ;
        W_SALT(j,i) = - w[I_W_SALT] ;
        W_THE(j,i)  = - w[I_W_THE] ;
      }
    }
  }

  return(0) ;
} 





int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  
  /*
    Initialization 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  {
    double c[Element_MaxNbOfDOF*Element_MaxNbOfDOF] ;
    int dec = TangentCoefficients(el,t,dt,c) ;
    double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }


  return(0) ;

#undef K
}



int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Conservation of the total mass: (M - M_n) + dt * div(W) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = M_TOT(i) - M_TOTn(i) ;
        } else {
          g[i*nn + j] = dt * W_TOT(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_Mass) -= r1[i] ;
      }
    }
  }
  
  /*
    Conservation of the salt mass: (M_SALT - M_SALTn) + dt * div(W_SALT) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = M_SALT(i) - M_SALTn(i) ;
        } else {
          g[i*nn + j] = dt * W_SALT(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_Salt) -= r1[i] ;
      }
    }
  }
  
  /*
    Entropy balance: (S - S_n) + dt * div(Q/T + S_l*W_l) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = S_TOT(i) - S_TOTn(i) ;
        } else {
          g[i*nn + j] = dt * W_THE(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_The) -= r1[i] ;
      }
    }
  }
  
  /*
   * Equality of the chemical potentials
   */
  #if defined (E_ChemP)
  {
    for(i = 0 ; i < nn ; i++) {
      double r1 = Mu_Liq(i) - Mu_Ice(i) ;
      
      R(i,E_ChemP) -= r1 ;
    }
  }
  #endif
  
  return(0) ;

#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nso = 16 ;
  int    i ;
  

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Variables */
    double* x = ComputeVariables(el,u,u,f,t,0,j) ;

    /* pressure */
    double p_l  =  x[I_P_L] ;
    /* temperature */
    double tem  =  x[U_tem] ;
    
    /* ice pressure */
    double p_i  = x[I_P_I] ;
    
    /* liquid saturation */
    double sd_l  = x[I_SD_L] ;
    double sd_i  = 1 - sd_l ;
    
    /* original Qiang !!!
    double K = k_s*(1 - phi*(1 + 3*k_s/(4*g_s))) ;
    double G = g_s*(1 - 5*phi*(4*g_s + 3*k_s)/(8*g_s + 9*k_s)) ;
    double b_i  = phi*sd_i*(1 + 3*k_s/(4*g_s)) ;
    double b_l  = phi*sd_l*(1 + 3*k_s/(4*g_s)) ;
    double N_ii = 3*phi*sd_i/(4*g_s) ;
    double N_ll = 3*phi*sd_l/(4*g_s) ; 
    double alpha_phi_l = alpha_s*(b_l - phi*sd_l) ;
    double alpha_phi_i = alpha_s*(b_i - phi*sd_i) ;
    double Epsi   = (b_l*p_l + b_i*p_i + alpha_s*(tem - T_m))/K ;
    double phi_l  = b_l*Epsi + p_l*N_ll - alpha_phi_l*(tem - T_m) ;
    double phi_i  = b_i*Epsi + p_i*N_ii - alpha_phi_i*(tem - T_m) ;
    double phi_p  = phi + phi_l + phi_i ;
    
    double Eps    = (b_l*(p_l - p_m) + b_i*(p_i - p_m) + alpha_s*K*(tem - T_m))/K ;
    double Eps1   = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double Eps2   = alpha_s*K*(tem - T_m)/K ;
    double Rho    = (phi_p<1) ? (phi_p*sd_l*(p_l - p_m) + phi_p*sd_i*(p_i - p_m))/(1-phi_p) : 0 ;
    double Eps_T  = alpha_s*K*(tem - T_m)/K ;
    double Eps_P  = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double v = 0.2 ;
    double Eps_xx = Eps*(1+v)/(1-v);
    double Rho_yy = 2*G*Eps_xx ;
    double Rho_zz = Rho_yy ;
    */
    
    /* Poroelastic coefficients (Mori-Tanaka)*/
    double K = (1 - phi)*k_s*g_s/(0.75*phi*k_s + g_s) ;
    double G = (1 - phi)*g_s*(9*k_s + 8*g_s)/(3*k_s*(3 + 2*phi) + 4*g_s*(2 + 3*phi)) ;
    double b    = 1 - K/k_s ;
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
    
    double dp_l = p_l - p0 ;
    double dp_i = p_i - p0 ;
    double dtem = tem - T0 ;
    
    /* Boundary conditions Sig_xx = Eps_yy = Eps_zz = 0 */
    double Eps_T  = alpha_s*K*dtem/k_oedo ;
    double Eps_P  = (b_l*dp_l + b_i*dp_i)/k_oedo ;
    double Eps_xx = Eps_P + Eps_T ;
    double Sig_yy = - 2*G*Eps_xx ;
    //double Sig_zz = Sig_yy ;
      
    double phi_l  = b_l*Eps_xx + N_ll*dp_l + N_il*dp_i - alpha_phi_l*dtem ;
    double phi_i  = b_i*Eps_xx + N_il*dp_l + N_ii*dp_i - alpha_phi_i*dtem ;
    //double phi_p  = phi + phi_l + phi_i ;
    
    
    
    if(Element_GetCoordinateSystem(el) != CARTESIAN) {
      Message_FatalError("Impossible") ;
    }
    
    i = 0 ;
    /* quantites exploitees */
    Result_Store(r + i++,x + I_P_L,"Liquid pressure",1) ;
    Result_Store(r + i++,x + U_c_s,"Salt concentration",1) ;
    Result_Store(r + i++,x + U_tem,"Temperature",1) ;
    Result_Store(r + i++,x + I_SD_L,"Liquid saturation",1) ;
    {
      double* w  = &M_TOT(0) ;
      double* wn = FVM_ComputeTheNodalFluxVector(fvm,w) ;
      
      Result_Store(r + i++,wn + 3*j,"Liquid mass flux vector",3) ;
    }
    {
      double* w  = &M_SALT(0) ;
      double* wn = FVM_ComputeTheNodalFluxVector(fvm,w) ;
      
      Result_Store(r + i++,wn + 3*j,"Salt mass flux vector",3) ;
    }
    Result_Store(r + i++,x + I_P_I,"Ice pressure",1) ;
    Result_Store(r + i++,x + I_RHO_L,"Solution density",1) ;
    Result_Store(r + i++,x + I_RHO_I,"Ice density",1) ;
    Result_Store(r + i++,&Eps_xx,"Strain_xx",1) ;
    Result_Store(r + i++,&Sig_yy,"Stress_yy",1) ;
    Result_Store(r + i++,&Eps_T,"Strain_temperature",1) ;
    Result_Store(r + i++,&Eps_P,"Strain_pressure",1) ;
    //Result_Store(r + i++,&phi_p,"porosity",1) ;
    Result_Store(r + i++,&phi_i,"Deformation of ice pores",1) ;
    Result_Store(r + i++,&phi_l,"Deformation of liquid pores",1) ;
    {
      double lnaw = LogActivityOfWater(x[U_c_s],tem) ;
      Result_Store(r + i++,&lnaw,"Lna_w",1) ;
    }
  }

  return (nso) ;
}



void  ComputeTransferCoefficients(Element_t* el,double** u,double* f)
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;

  /*
    Input data
  */
  GetProperties(el) ;

  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;

  /* Transfer coefficients */
  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;
    
    double p_l = x[I_P_L] ;
    double tem = x[U_tem] ;  
    double c_s = x[U_c_s] ;
    double p_i = x[I_P_I] ;
    double p_c = p_i - p_l ;

    double sd_l = x[I_SD_L] ;
    double sd_i = 1 - sd_l ;
    
    /* Mass densities */
    double rho_l     = x[I_RHO_L] ;

    /* Darcy */
    double mu_l = WaterViscosity(tem) ;
    double kr_l = RelativePermeabilityToLiquid(p_c) ;
    double kd_l = rho_l*k_int/mu_l*kr_l ;
    
    /* Fick */
    double tau_l   = Tortuosity(phi,sd_l) ;
    double kf_salt = M_Salt*phi*sd_l*tau_l*D_salt ;
    
    /* Fourier */
    double lam_h2o = sd_l*lam_l + sd_i*lam_i ;
    double lam_hom = lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
    double kth = lam_hom/tem ;

    KD_LIQ[i]  = kd_l ;
    KF_SALT[i] = kf_salt ;
    KC_SALT[i] = M_Salt*c_s/rho_l ;
    KTH[i]     = kth ;
    S_L[i]     = x[I_S_L] ;
  }
  
}



double* ComputeVariableFluxes(Element_t* el,int i,int j)
{
  int nn = Element_GetNbOfNodes(el) ;
  double* grdij = dVariables ;
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij  = dist[nn*i + j] ;


  {
    double* xi  = Variables[i] ;
    double* xj  = Variables[j] ;
    int k ;
      
    for(k = 0 ; k < NbOfVariables ; k++)  {
      grdij[k] = (xj[k] - xi[k])/dij ;
    }
  }
  
  /* Fluxes */
  {
    double* w = ComputeFluxes(el,grdij,i,j) ;
    
    return(w) ;
  }
}





double* ComputeFluxes(Element_t* el,double* grd,int i,int j)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes ;

  /* Gradients */
  double grd_p_l      = grd[I_P_L] ;
  double grd_c_s      = grd[U_c_s] ;
  double grd_tem      = grd[U_tem] ;
  
  /* Transfer terms */
  double kd_l    = 0.5*(KD_LIQ[i]  + KD_LIQ[j]) ;
  double kf_salt = 0.5*(KF_SALT[i] + KF_SALT[j]) ;
  double kc_salt = 0.5*(KC_SALT[i] + KC_SALT[j]) ;
  double kth     = 0.5*(KTH[i]     + KTH[j]) ;
  double s_l     = 0.5*(S_L[i]     + S_L[j]) ;

  /* Fluxes */
  double w_l    = - kd_l * grd_p_l ;
  double j_salt = - kf_salt * grd_c_s ;
  double w_salt =   kc_salt * w_l + j_salt ;
  double q      = - kth * grd_tem ;
  double w_the  =   q + s_l * w_l ;
  

  w[I_W_TOT ]  = w_l ;
  w[I_W_SALT]  = w_salt ;
  w[I_W_THE]   = w_the ;
  
  return(w) ;
}



int TangentCoefficients(Element_t* el,double t,double dt,double* c)
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
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < NEQ ; i++) {
    dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
  }

  
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x   = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
    int k ;
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double* dx    = ComputeVariableDerivatives(el,t,dt,x,dxk,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
        
        cii[E_Mass*NEQ   + k] = dx[I_M_TOT ] ;
        cii[E_Salt*NEQ   + k] = dx[I_M_SALT] ;
        cii[E_The*NEQ    + k] = dx[I_S_TOT ] ;
      #if defined (E_ChemP)
        cii[E_ChemP*NEQ  + k] = dx[I_MU_L] - dx[I_MU_I] ;
      #endif
      }
      
      /* Transfer terms from node i to node j: d(w_ij)/d(u_i) */
      {
        int j ;
        
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            double* dw = ComputeFluxes(el,dx,i,j) ;
            double* cij = c + (i*nn + j)*NEQ*NEQ ;
            double dij  = dist[nn*i + j] ;
            double dtdij = dt/dij ;
        
            cij[E_Mass*NEQ   + k] = - dtdij * dw[I_W_TOT ] ;
            cij[E_Salt*NEQ   + k] = - dtdij * dw[I_W_SALT] ;
            cij[E_The*NEQ    + k] = - dtdij * dw[I_W_THE ] ;
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int n)
{
  double* x = Variables[n] ;
  
  /* Primary Variables */
#if defined (U_p_l)
  x[U_p_l] = P_l(n) ;
#endif
#if defined (U_p_i)
  x[U_p_i] = P_i(n) ;
#endif
#if defined (U_p_max)
  x[U_p_max] = P_max(n) ;
#endif
  x[U_tem] = TEM(n) ;
  x[U_c_s] = C_s(n) ;
  
  /* Needed variables to compute secondary components */
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  /* Primary variables */
#if defined (U_p_l)
  double p_l    = x[U_p_l] ;
#endif
#if defined (U_p_i)
  double p_i    = x[U_p_i] ;
#endif
#if defined (U_p_max)
  double p_max    = x[U_p_max] ;
#endif
  double c_s    = x[U_c_s] ;
  double tem    = x[U_tem] ;
  
  double lna = LogActivityOfWater(c_s,tem) ;
  
  /* Pressures */
#if !defined (U_p_i) && defined (U_p_l)
  double p_i = p_m + (V_H2O*(p_l - p_m) - S_m*(tem - T_m) + R_g*tem*lna)/V_Ice ;
#endif
#if !defined (U_p_l) && defined (U_p_i)
  double p_l = p_m + (V_Ice*(p_i - p_m) + S_m*(tem - T_m) - R_g*tem*lna)/V_H2O ;
#endif
#if defined (U_p_max)
  //     0 = V_Ice*(p_i   - p_m) - V_H2O*(p_l   - p_m) + S_m*(tem - T_m) - R_g*tem*lna
  double c = V_Ice*(p_max - p_m) - V_H2O*(p_max - p_m) + S_m*(tem - T_m) - R_g*tem*lna ;
  double p_min = (c > 0) ? p_max - c/V_Ice : p_max + c/V_H2O ;
  double p_i   = (c > 0) ? p_min : p_max ;
  double p_l   = (c > 0) ? p_max : p_min ;
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
  
  /* Chemical potentials */
  double mu_l = V_H2O*(p_l - p_m) + R_g*tem*lna ;
  double mu_i = V_Ice*(p_i - p_m) + S_m*(tem - T_m) ;
  
  /* Backup */
  x[I_M_TOT   ] = m_l + m_i ;
  x[I_M_SALT  ] = m_salt ;
  x[I_S_TOT   ] = s_tot ;
  x[I_S_L     ] = s_l ;
  x[I_SD_L    ] = sd_l ;
  x[I_P_L     ] = p_l ;
  x[I_P_I     ] = p_i ;
  x[I_RHO_L   ] = rho_l ;
  x[I_RHO_I   ] = rho_i ;
  x[I_MU_L    ] = mu_l ;
  x[I_MU_I    ] = mu_i ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double* dx = dVariables ;
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


double activite_w_ideal(double c_s,double tem)
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
