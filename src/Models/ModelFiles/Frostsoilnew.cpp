#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the finite element method */
#include "FEM.h"
#include "Elasticity.h"

#define TITLE "Frost actions in soil"
#define AUTHORS "Li-Dangla"

#include "PredefinedMethods.h"


/* Nb of equations of the model */
#define NEQ   (dim+3)


/* Indices of equations */
#define E_MECH   (3)
#define E_MASS   (0)
#define E_SALT   (1)
#define E_THER   (2)

/* Indices of unknowns (generic indices) */
#define U_MECH   E_MECH
#define U_MASS   E_MASS
#define U_SALT   E_SALT
#define U_THER   E_THER


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */
 
/* Mechanics */
#define U_DISP   U_MECH

/* Mass */
//#define U_P_MAX  U_MASS
#define U_P_L    U_MASS

/* Salt */
#define U_C_SALT    U_SALT

/* Thermic */
#define U_TEMP    U_THER




/* We define some names for implicit terms */
#define NVI   (55)    /*  nb of implicit terms per point */
template<typename T = double>
struct ImplicitValues_t {
  T Displacement[3];
  T Strain[9];
  T U_mass;
  T GradU_mass[3];
  T Concentration_salt;
  T GradConcentration_salt[3];
  T Temperature;
  T GradTemperature[3];
  T Mass_total;
  T MassFlow_total[3];
  T Mass_salt;
  T MassFlow_salt[3];
  T Entropy_total;
  T HeatFlow[3];
  T Stress[9];
  T BodyForce[3];
  T SaturationDegree_liquid;
  T Pressure_liquid;
  T Pressure_ice;
  T MassDensity_liquid;
  T MassDensity_ice;
  T ChemicalPotential_liquid;
  T ChemicalPotential_ice;
} ;



/* We define some names for explicit terms */
#define NVE   (5)     /*  nb of explicit terms per point */
template<typename T = double>
struct ExplicitValues_t {
  T Permeability_liquid;
  T Diffusion_salt;
  T ThermalConductivity;
  T MassFraction_salt;
  T Entropy_liquid;
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
#define NV0   (0)     /*  nb of constant terms per point */
template<typename T = double>
struct ConstantValues_t {};



#include <CustomValues.h>

using Values_t = CustomValues_t<ImplicitValues_t<>,ExplicitValues_t<>,ConstantValues_t<>> ;



/* Shorthands of some units */
#include "InternationalSystemOfUnits.h"

#define m1       (InternationalSystemOfUnits_OneMeter)
#define m3       (m1*m1*m1)
#define dm       (0.1*m1)
#define cm       (0.01*m1)
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
#define LAM_l       ((0.6) * Watt / m1 / Kelvin)
static double lam_l ;
/* Thermal conductivity of ice (W/m/K) */
#define LAM_i       ((2.2) * Watt / m1 / Kelvin)
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

#define ScaleFactor(EPSV,A,C) \
        (((EPSV) > 0) ? 1 + ((C) - 1)*(1 - exp(-(A)*(EPSV))) : 1)

#define aa  (1000)
//#define capillarypressurefraction  (10)
#define CapillaryPressureFactor(EPSV) \
        ScaleFactor(EPSV,aa,capillarypressurefraction)

#define YoungModulus(EPSV) \
        young*ScaleFactor(EPSV,aa,youngfraction)

#define Porosity(EPSV) \
        phi0*(thickness + Math_Max(EPSV,0))



/* The parameters below are read in the input data file */
static double gravity ;
static double rho_s ;
static double phi0 ;
static double k_int ;
static double lam_s ;
static double C_s ;
static double alpha_s ;
static double p0 ;
static double T0 ;
static double* cijkl ;
static Elasticity_t* elasty ;
static double young ;
static double poisson ;
static double biot ;
static double thickness ;
static double youngfraction ;
static double capillarypressurefraction ;

static double  effectivediffusioncoef(double,double) ;
double effectivediffusioncoef(double phi,double s)
{
  double tau_l_sat = 0.296e-3*exp(9.95*phi)/phi ;
  double tau = 0 ;
  
  if(s > 0.) {
    tau = tau_l_sat/(s*(1 + 625*pow(1 - s,4))) ;
  } else {
    tau = 0 ;
  }
  
  return(thickness*phi*s*tau*D_salt) ;
}

static double thermalconductivity(double,double) ;
double thermalconductivity(double phi,double sd_l)
{
  double sd_i = 1 - sd_l ;
  double lam_h2o = sd_l*lam_l + sd_i*lam_i ;
  double lam_hom = lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
  
  return(thickness*lam_hom) ;
}



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
  else if(strcmp(s,"Young") == 0)    return (9) ;
  else if(strcmp(s,"Poisson") == 0)  return (10) ;
  else if(strcmp(s,"Biot") == 0)     return (11) ;
  else if(strcmp(s,"thickness") == 0) return (12) ;
  else if(strcmp(s,"YoungFraction") == 0)    return (13) ;
  else if(strcmp(s,"gravity") == 0)  return (14) ;
  else if(strcmp(s,"rho_s") == 0)    return (15) ;
  else return (-1) ;
}


/* They are retrieved automatically by calling the following function */
static void    GetProperties(Element_t*) ;
void GetProperties(Element_t* el)
{
/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
  gravity = GetProperty("gravity") ;
  rho_s   = GetProperty("rho_s") ;
  phi0    = GetProperty("porosity") ;
  k_int   = GetProperty("k_int") ;
  C_s     = GetProperty("C_s") ;
  lam_s   = GetProperty("lam_s") ;
  alpha_s = GetProperty("alpha_s") ;
  p0      = GetProperty("p0") ;
  T0      = GetProperty("T0") ;
  young   = GetProperty("Young") ;
  poisson = GetProperty("Poisson") ;
  biot    = GetProperty("Biot") ;
  
  thickness = 1 ;
  youngfraction  = 1 ;
  capillarypressurefraction = 1 ;
  
  if(Element_HasZeroThickness(el)) {
    thickness = GetProperty("thickness") ;
    youngfraction  = GetProperty("YoungFraction") ;
    capillarypressurefraction = 10 ;

    young    /= thickness ;
    k_int    *= thickness ;
    C_s      *= thickness ;
    alpha_s  *= thickness ;
  }
  
  elasty  = Element_FindMaterialData(el,Elasticity_t,"Elasticity") ;
  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
#undef GetProperty
}



/* Functions used below */
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
  //S_m = 0 ;
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


#include <ConstitutiveIntegrator.h>
#include <LocalVariables.h>

#define USEFUNCTOR 0

#if USEFUNCTOR
struct SetInputs_t {
  Values_t& operator()(const LocalVariables_t<Values_t>&,const int&,Values_t&) ;
} ;

struct Integrate_t {
  Values_t&  operator()(const Element_t*,const double&,const double&,const Values_t&,Values_t&) ;
} ;

using CI_t = ConstitutiveIntegrator_t<Values_t,SetInputs_t,Integrate_t>;
#else
using SetInputs_t = LocalVariables_SetInputs_t<Values_t> ;
using Integrate_t = ConstitutiveIntegrator_Integrate_t<Values_t> ;
using CI_t = ConstitutiveIntegrator_t<Values_t,SetInputs_t*,Integrate_t*>;
#endif


static Integrate_t Integrate;
static SetInputs_t SetInputs;



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
  Model_CopyNameOfEquation(model,E_MASS,"mass") ;
  Model_CopyNameOfEquation(model,E_SALT,"salt") ;
  Model_CopyNameOfEquation(model,E_THER,"the") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
#if defined (U_P_MAX)
  Model_CopyNameOfUnknown(model,U_MASS,"p_max") ;
#elif defined (U_P_L)
  Model_CopyNameOfUnknown(model,U_MASS,"p_l") ;
#endif
  Model_CopyNameOfUnknown(model,U_SALT,"c_s") ;
  Model_CopyNameOfUnknown(model,U_THER,"tem") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_MECH + i,name_unk[i]) ;
  }
  
  Model_GetComputePropertyIndex(model) = pm ;
  
  Model_GetSequentialIndexOfUnknown(model)[E_THER] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_MASS] = 1 ;
  Model_GetSequentialIndexOfUnknown(model)[E_SALT] = 1 ;
  for(i = 0 ; i < dim ; i++) {
    Model_GetSequentialIndexOfUnknown(model)[E_MECH+i] = 1 ;
  }

  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 16 ;

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
    Material_GetProperty(mat)[pm("Young")] = -1 ;
    Material_GetProperty(mat)[pm("YoungFraction")] = 1 ;
    Material_GetProperty(mat)[pm("thickness")] = 1 ;
    Material_GetProperty(mat)[pm("p0")]  = p_m ;
    Material_GetProperty(mat)[pm("T0")]  = T_m ;
  }
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  if(Material_GetProperty(mat)[pm("Young")] < 0) {
    double k_s  = Material_GetPropertyValue(mat,"k_s") ;
    double g_s  = Material_GetPropertyValue(mat,"g_s") ;
    
    phi0  = Material_GetPropertyValue(mat,"porosity") ;
    
    if(k_s > 0 && g_s > 0) {
      /* Elastic moduli (Mori-Tanaka) */
      double k_0 = k_s ;
      double g_0 = g_s ;
      double a_0 = 0.75 * k_0 / g_0 ;
      double b_0 = 6 * (k_0 + 2 * g_0) / (9 * k_0 + 8 * g_0) ;
      double K = (1 - phi0) * k_s / (1 + phi0 * a_0 * k_s / k_0) ;
      double G = (1 - phi0) * g_s / (1 + phi0 * b_0 * g_s / g_0) ;
      
      young =  9 * G * K / (3 * K + G) ;
      poisson = 0.5 * (1 - young / (3 * K)) ;
      //poisson = (1.5 * K - G) / (3 * K + G) ;
      biot = 1 - K/k_s ;
      
      Material_GetProperty(mat)[pm("Young")] = young ;
      Material_GetProperty(mat)[pm("Poisson")] = poisson ;
      Material_GetProperty(mat)[pm("Biot")] = biot ;
    } else {
      arret("ReadMatProp") ;
    }
  }
      
  /* Elasticity */
  {
    elasty = Elasticity_Create() ;
      
    Material_AppendData(mat,1,elasty,Elasticity_t,"Elasticity") ;
  }
  
  /* The 4th rank elastic tensor */
  {
    elasty = Material_FindData(mat,Elasticity_t,"Elasticity") ;

    /* isotropic Hooke's law */
    { 
      young =  Material_GetProperty(mat)[pm("Young")] ;
      poisson = Material_GetProperty(mat)[pm("Poisson")] ;
    
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
#if defined (U_P_MAX)
  printf("\t 1. The max of liquid and ice pressures   (p_max)\n") ;
#elif defined (U_P_L)
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


  /* Continuity of unknowns across zero-thickness element */
  {
    if(Element_HasZeroThickness(el)) {
      #if defined (U_P_MAX)
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"p_max") ;
      #elif defined (U_P_L)
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"p_l") ;
      #endif
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"mass") ;
      
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"tem") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"the") ;
      
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"c_s") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"salt") ;
    }
  }
  
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_THER) {
      
      for(i = 0 ; i < nn ; i++) {
        double tem = Element_GetValueOfNodalUnknown(el,u,i,U_THER) ;
        
        R(i,E_THER) /= tem ;
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
    
    CI_t ci(el,t,0,u,u,vi0,SetInputs,Integrate) ;

    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      Values_t& val = ci.Integrate(p) ;
      
      /* storage in vi */
      ci.StoreImplicitTerms(p,vi0) ;

      /* storage in ve */
      ci.StoreExplicitTerms(p,ve0) ;
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
  if(NVE) {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    CI_t ci(el,t,0,u,u,vi_n,SetInputs,Integrate) ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      Values_t& val = ci.Integrate(p) ;
    
      /* storage in ve */
      ci.StoreExplicitTerms(p,ve0) ;
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
    
    CI_t ci(el,t,dt,u,u_n,vi_n,SetInputs,Integrate) ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      Values_t& val = ci.Integrate(p) ;
      
      /* storage in vi */
      ci.StoreImplicitTerms(p,vi0) ;
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
  int dim = Element_GetDimensionOfSpace(el) ;
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
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,3,E_MECH) ;
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
    double* k_mass    = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    double* k_salt_pl = FEM_ComputeConductionMatrix(fem,intfct,c+9*3,dec) ;
    double* k_salt    = FEM_ComputeConductionMatrix(fem,intfct,c+9*4,dec) ;
    double* k_the_pl  = FEM_ComputeConductionMatrix(fem,intfct,c+9*6,dec) ;
    double* k_the     = FEM_ComputeConductionMatrix(fem,intfct,c+9*8,dec) ;
    int    i ;
  
    for(i = 0 ; i < nn ; i++) {
      int    j ;
      
      for(j = 0 ; j < nn ; j++) {
        K(E_MASS + i*NEQ,U_MASS + j*NEQ) += dt*k_mass[i*nn + j] ;
        K(E_SALT + i*NEQ,U_MASS + j*NEQ) += dt*k_salt_pl[i*nn + j] ;
        K(E_SALT + i*NEQ,U_SALT + j*NEQ) += dt*k_salt[i*nn + j] ;
        K(E_THER  + i*NEQ,U_MASS + j*NEQ) += dt*k_the_pl[i*nn + j] ;
        K(E_THER  + i*NEQ,U_THER  + j*NEQ) += dt*k_the[i*nn + j] ;
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
  double* vi1   = Element_GetCurrentImplicitTerm(el) ;
  double* vi1_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  using TI = CustomValues_TypeOfImplicitValues(Values_t);
  TI& val = *((TI*) vi1) ;
  TI& val_n = *((TI*) vi1_n) ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }

  if(Element_IsSubmanifold(el)) return(0) ;
  

  /* Compute here the residu R(n,i) */
  

  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  if(Element_EquationIsActive(el,E_MECH)) 
  {
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,val.Stress,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_MECH + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  #if 1
  if(Element_EquationIsActive(el,E_MECH)) 
  {
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,val.BodyForce + dim - 1,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_MECH + dim - 1) -= -rbf[i] ;
    }
    
  }
  #endif
  
  
  
  /* 2. Conservation of total mass */
  
  /* 2.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_MASS)) 
  {
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++) {
      g1[i] = (&val.Mass_total)[i*NVI] - (&val_n.Mass_total)[i*NVI] ;
    }
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_MASS) -= ra[i] ;
    }
  }
  
  /* 2.2 Transport Terms */
  if(Element_EquationIsActive(el,E_MASS)) 
  {
    double* rf = FEM_ComputeFluxResidu(fem,intfct,val.MassFlow_total,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_MASS) -= -dt*rf[i] ;
  }
  
  
  
  /* 3. Conservation of salt */
  
  /* 3.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_SALT)) 
  {
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++) {
      g1[i] = (&val.Mass_salt)[i*NVI] - (&val_n.Mass_salt)[i*NVI] ;
    }
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_SALT) -= ra[i] ;
    }
  }
  
  /* 3.2 Transport Terms */
  if(Element_EquationIsActive(el,E_SALT)) 
  {
    double* rf = FEM_ComputeFluxResidu(fem,intfct,val.MassFlow_salt,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_SALT) -= -dt*rf[i] ;
  }
  
  
  
  /* 4. Entropy balance */
  
  /* 3.1 Accumulation Terms */
  if(Element_EquationIsActive(el,E_THER)) 
  {
    double* vi = vi1 ;
    double* vi_n = vi1_n ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++) {
      g1[i] = (&val.Entropy_total)[i*NVI] - (&val_n.Entropy_total)[i*NVI] ;
    }
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_THER) -= ra[i] ;
    }
  }
  
  /* 3.2 Transport Terms */
  if(Element_EquationIsActive(el,E_THER)) 
  {
    double* vi = vi1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,val.HeatFlow,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_THER) -= -dt*rf[i] ;
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
  double** u_n  = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  //int dim = Element_GetDimensionOfSpace(el) ;
  FEM_t* fem = FEM_GetInstance(el) ;

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
  
  
  CI_t ci(el,t,0,u,u,vi,SetInputs,Integrate) ;
  
  
  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Variables */
    Values_t& val = ci.Integrate(p) ;
    /* strains */
    double* eps = val.Strain ;
    double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_MECH) ;
    double  epsv = eps[0] + eps[4] + eps[8] ;
    double  epsv_n = eps_n[0] + eps_n[4] + eps_n[8] ;
    /* pressures */
    double p_l  = val.Pressure_liquid ;
    double p_i  = val.Pressure_ice ;
    double dp_l = val.Pressure_liquid - p0 ;
    double dp_i = val.Pressure_ice - p0 ;
    /* temperature */
    double dtem = val.Temperature - T0 ;
    /* fluxes */
    double w_tot[3]  = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double w_the[3]  = {0,0,0} ;
    double young1 = YoungModulus(epsv_n) ;
    double K = young1 / (3 - 6*poisson) ;
    double G = 0.5 * young1 / (1 + poisson) ;
    double b = biot ;
    double sd_l = val.SaturationDegree_liquid ;
    double sd_i = 1 - sd_l ;
    double b_i  = b*sd_i ;
    double b_l  = b*sd_l ;
    double N_il = (K > 0) ? - (b - phi0) / K * b * sd_i * sd_l  : 0 ;
    double N_ii = (K > 0) ?   (b - phi0) / K * sd_i * (1 - b_i) : 0 ;
    double N_ll = (K > 0) ?   (b - phi0) / K * sd_l * (1 - b_l) : 0 ;
    double alpha_phi = (b - phi0)*alpha_s ;
    double alpha_phi_l = sd_l*alpha_phi ;
    double alpha_phi_i = sd_i*alpha_phi ;
    
    double k_oedo = (K + 4/3.*G) ; /* k_oedo Oedometric modulus */
    double Kk_oedo = (1 + poisson)/(3 - 3*poisson) ; /* K/k_oedo */
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    Values_t* val1 = (Values_t*) vi ;
    
    int    i ;
    
    /* Averaging */
    for(i = 0 ; i < np ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_tot[j]  += val1->MassFlow_total[i*NVI+j]/np ;
      for(j = 0 ; j < 3 ; j++) w_salt[j] += val1->MassFlow_salt[i*NVI+j]/np ;
      for(j = 0 ; j < 3 ; j++) w_the[j]  += val1->HeatFlow[i*NVI+j]/np ;
      for(j = 0 ; j < 9 ; j++) sig[j]    += val1->Stress[i*NVI+j]/np ;
    }
    
    i = 0  ;
    /* quantites exploitees */
    Result_Store(r + i++,&val.Pressure_liquid,"Liquid pressure",1) ;
    Result_Store(r + i++,&val.Concentration_salt,"Salt concentration",1) ;
    Result_Store(r + i++,&val.Temperature,"Temperature",1) ;
    Result_Store(r + i++,&val.SaturationDegree_liquid,"Liquid saturation",1) ;
    Result_Store(r + i++,w_tot,"Liquid mass flux vector",3) ;
    Result_Store(r + i++,w_salt,"Salt mass flux vector",3) ;
    Result_Store(r + i++,&val.Pressure_ice,"Ice pressure",1) ;
    Result_Store(r + i++,&val.MassDensity_liquid,"Solution density",1) ;
    Result_Store(r + i++,&val.MassDensity_ice,"Ice density",1) ;
    Result_Store(r + i++,val.Strain,"Strain tensor",9) ;
    Result_Store(r + i++,sig,"Stress tensor",9) ;
    Result_Store(r + i++,val.Displacement,"Displacement vector",3) ;
      
    {
      double Eps_T  = alpha_s*Kk_oedo*dtem ;
      
      Result_Store(r + i++,&Eps_T,"Strain_temperature",1) ;
    }
    
    {
      double Eps_P  = (k_oedo > 0) ? (b_l*dp_l + b_i*dp_i)/k_oedo : 0 ;
      
      Result_Store(r + i++,&Eps_P,"Strain_pressure",1) ;
    }
    
    {
      double phi_i  = b_i*epsv + N_il*dp_l + N_ii*dp_i - alpha_phi_i*dtem ;
      
      Result_Store(r + i++,&phi_i,"Deformation of ice pores",1) ;
    }
    
    {
      double phi_l  = b_l*epsv + N_ll*dp_l + N_il*dp_i - alpha_phi_l*dtem ;
    
      Result_Store(r + i++,&phi_l,"Deformation of liquid pores",1) ;
    }
    
    {
      double lnaw = LogActivityOfWater(val.Concentration_salt,val.Temperature) ;
      
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
    
    CI_t ci(el,t,dt,u,u_n,vi_n,SetInputs,Integrate) ;
    
    for(p = 0 ; p < np ; p++) {
      /* Variables */
      Values_t& val = ci.Integrate(p) ;
      double* c0 = c + p*dec ;
      #define Index(V)  CustomValues_Index(&val,V,double)


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
         
          /* Already updated in ComputeVariables */
          for(i = 0 ; i < 81 ; i++) {
            c1[i] = cijkl[i] ;
          }

          //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
        }
        
        /* Coupling matrices */
        {
          double deps = 1.e-6 ;
          int i_eps = Index(Strain[0]);
          Values_t& dval = ci.Differentiate(deps,i_eps) ;
          double* c00 = c0 + 81 + 3*9 ;
          
          /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
           * and zero for the derivatives w.r.t others */
          {
            double* c1 = c00 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass_total ;
          }
          {
            double* c1 = c00 + 12 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass_salt ;
          }
          {
            double* c1 = c00 + 2*12 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dval.Entropy_total ;
          }
        }
      }
      
      /* Derivatives w.r.t U_MASS */
      {
        double  dp = dui[U_MASS] ;
        int i_u_mass = Index(U_mass);
        Values_t& dval = ci.Differentiate(dp,i_u_mass) ;
        
        /* Tangent Biot's coefficient */
        {
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dval.Mass_total ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dval.Mass_salt ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dval.Entropy_total ;
          }
        }
      }
      
      /* Derivatives w.r.t U_SALT */
      {
        double  dc = dui[U_SALT] ;
        int i_c_salt = Index(Concentration_salt);
        Values_t& dval = ci.Differentiate(dc,i_c_salt) ;
        
        /* Tangent Biot's coefficient */
        {
          double* c1 = c0 + 81 + 9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 + 1 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dval.Mass_total ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dval.Mass_salt ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dval.Entropy_total ;
          }
        }
      }
      
      /* Derivatives w.r.t  U_THER */
      {
        double  dtem = dui[U_THER] ;
        int i_tem = Index(Temperature);
        Values_t& dval = ci.Differentiate(dtem,i_tem) ;
        
        /* Tangent Biot's coefficient */
        {
          double* c1 = c0 + 81 + 2*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 3*9 + 9 + 2 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dval.Mass_total ;
          }
          {
            double* c1 = c00 + 12 ;
        
            c1[0] = dval.Mass_salt ;
          }
          {
            double* c1 = c00 + 2*12 ;
        
            c1[0] = dval.Entropy_total ;
          }
        }
      }
      #undef Index
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
  double* vi0 = Element_GetPreviousImplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 * 9 ;
  int    p ;
  using TE = CustomValues_TypeOfExplicitValues(Values_t);
  TE* val0 = (TE*) ve0 ;
  

  for(p = 0 ; p < np ; p++) {
    int i ;
    double* c1 = c + p*dec ;
    TE& val = val0[p] ;
    
    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
    
    {
      /* Transfer coefficients for the liquid mass flux */
      {
        double* c2 = c1 ;
        
        c2[0] = val.Permeability_liquid ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      
      /* Transfer coefficients for the salt mass flux */
      {
        double* c2 = c1 + 3 * 9 ;
        
        c2[0] = val.MassFraction_salt*val.Permeability_liquid ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      {
        double* c2 = c1 + 4 * 9 ;
        
        c2[0] = val.Diffusion_salt ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      
      /* Transfer coefficients for the entropy flux */
      {
        double* c2 = c1 + 6 * 9 ;
        
        c2[0] = val.Entropy_liquid*val.Permeability_liquid ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      {
        double* c2 = c1 + 8 * 9 ;
        
        c2[0] = val.ThermalConductivity ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
    }
  }

  return(dec) ;
}





#if USEFUNCTOR
Values_t& SetInputs_t::operator()(const LocalVariables_t<Values_t>& var,const int& p,Values_t& val)
#else
Values_t& SetInputs(const LocalVariables_t<Values_t>& var,const int& p,Values_t& val)
#endif
{
  /* Displacements and strains */
  var.DisplacementVectorAndStrainFEM(p,U_MECH,val.Displacement) ;
    
  /* Pressure and pressure gradient */
  var.ValueAndGradientFEM(p,U_MASS,&val.U_mass) ;
  
  /* Salt concentration and concentration gradient */
  var.ValueAndGradientFEM(p,U_SALT,&val.Concentration_salt) ;
  
  /* Temperature and temperature gradient */
  var.ValueAndGradientFEM(p,U_THER,&val.Temperature) ;
  
  return(val) ;
}



#if USEFUNCTOR
Values_t& Integrate_t::operator()(const Element_t* el,const double& t,const double& dt,const Values_t& val_n,Values_t& val)
#else
Values_t&  Integrate(const Element_t* el,const double& t,const double& dt,const Values_t& val_n,Values_t& val)
#endif
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Retrieve the primary variables */
  /* Strains */
  double* eps   = val.Strain ;
  double* eps_n = val_n.Strain ;
  double  epsv  = eps[0] + eps[4] + eps[8] ;
  double  epsv_n = eps_n[0] + eps_n[4] + eps_n[8] ;
  /* Salt concentration */
  double c_s    = val.Concentration_salt ;
  /* Temperature */
  double tem    = val.Temperature ;
  /* Log of water activity */
  double lna = LogActivityOfWater(c_s,tem) ;
  
  /* Pressures */
#if defined (U_P_MAX)
  double  p_max = val.U_mass ;
  //     0 = V_Ice*(p_i   - p_m) - V_H2O*(p_l   - p_m) + S_m*(tem - T_m) - R_g*tem*lna
  double c = V_Ice*(p_max - p_m) - V_H2O*(p_max - p_m) + S_m*(tem - T_m) - R_g*tem*lna ;
  double p_min = (c > 0) ? p_max - c/V_Ice : p_max + c/V_H2O ;
  double p_i   = (c > 0) ? p_min : p_max ;
  double p_l   = (c > 0) ? p_max : p_min ;
#elif defined (U_P_L)
  double p_l   = val.U_mass ;
  double p_i   = p_m + (V_H2O*(p_l - p_m) - S_m*(tem - T_m) + R_g*tem*lna)/V_Ice ;
#endif
  double p_c = p_i - p_l ;
  
  /* Porosity */
  double phi1 = Porosity(epsv) ;
  
  /* Saturations */
  double sd_l = SaturationDegree(p_c*CapillaryPressureFactor(epsv)) ;
  double sd_i = 1 - sd_l ;
  
  /* Backup */
    
  val.SaturationDegree_liquid = sd_l ;
  val.Pressure_liquid = p_l ;
  val.Pressure_ice = p_i ;
  

  /* Mass contents and entropies */
  {
    /* Mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_i     = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_Salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - vr_salt*rho_salt) ;
  
    /* Mass contents */
    double m_l     = phi1*sd_l*rho_l ;
    double m_i     = phi1*sd_i*rho_i ;
    double m_salt  = phi1*sd_l*rho_salt ;
  
    /* Entropies */
    double s_l   = C_l*log(tem/T_m) - alpha_l/rho_h2o_l0*(p_l - p_m) ;
    double s_i   = C_i*log(tem/T_m) - alpha_i/rho_h2o_i0*(p_i - p_m) - S_m/M_H2O ;
    double s_sol = C_s*log(tem/T_m) ;
    double s_tot = s_sol + m_l*s_l + m_i*s_i ;
  
    /*
    if(dphi > 0) {
      printf("dphi = %e\n",dphi) ;
    }
    */
    val.Mass_total = m_l + m_i ;
    val.Mass_salt  = m_salt ;
  
    val.Entropy_total = s_tot ;
    val.Entropy_liquid = s_l ;
    
    val.MassDensity_liquid = rho_l ;
    val.MassDensity_ice = rho_i ;
  }
  
  /* Chemical potentials */
  {
    double mu_l = V_H2O*(p_l - p_m) + R_g*tem*lna ;
    double mu_i = V_Ice*(p_i - p_m) + S_m*(tem - T_m) ;
    
    val.ChemicalPotential_liquid = mu_l ;
    val.ChemicalPotential_ice = mu_i ;
  }
  
  /* Stresses */
  {
    double young1 = YoungModulus(epsv_n) ;
    double K = young1 / (3 - 6*poisson) ;
    double b    = biot ;
    double b_i  = b*sd_i ;
    double b_l  = b*sd_l ;
    
    double dp_l = p_l - p0 ;
    double dp_i = p_i - p0 ;
    double dtem = tem - T0 ;
    
    double* sig = val.Stress ;
          
    Elasticity_SetParameters(elasty,young1,poisson) ;
    Elasticity_UpdateStiffnessTensor(elasty) ;
    
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
    }
    
    {
      sig[0] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
      sig[4] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
      sig[8] += - b_l * dp_l - b_i * dp_i - alpha_s * K * dtem ;
    }
  }
  
  /* Fluxes */
  {
    /* Transfer coefficients at the previous time */
    double  kd_liq  = val_n.Permeability_liquid ;
    double  kf_salt = val_n.Diffusion_salt ;
    double  kth     = val_n.ThermalConductivity ;
    double  mc_salt = val_n.MassFraction_salt ;
    double  s_l     = val_n.Entropy_liquid ;
    
    /* Fluxes */
    double* w_tot   = val.MassFlow_total ;
    double* w_salt  = val.MassFlow_salt ;
    double* w_the   = val.HeatFlow ; 
    
    /* Gradients */
    double* grd_c_s   = val.GradConcentration_salt ;
    double* grd_tem   = val.GradTemperature ;
    int i ;
    
#if defined (U_P_MAX)
    double  grd_p_l[3] ;
    double* grd_p_max = val.GradU_mass ;  
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    double dc_s = 1.e-2 * ObVal_GetValue(obval + U_SALT) ;
    double dtlnadc_s = (tem * LogActivityOfWater(c_s + dc_s,tem) - tem * lna) / dc_s ;
    double dtem = 1.e-2 * ObVal_GetValue(obval + U_THER) ;
    double dtlnadtem = ((tem + dtem) * LogActivityOfWater(c_s,tem + dtem) - tem * lna) / dtem ;
    
    for(i = 0 ; i < 3 ; i++) {
      double grd_tlna = dtlnadc_s * grd_c_s[i] + dtlnadtem * grd_tem[i] ;
      double grd_c    = (V_Ice - V_H2O) * grd_p_max[i] + S_m * grd_tem[i] - R_g * grd_tlna ;
      double grd_p_min = (c > 0) ? grd_p_max[i] - grd_c/V_Ice : grd_p_max[i] + grd_c/V_H2O ;
      
      grd_p_l[i]   = (c > 0) ? grd_p_max[i] : grd_p_min ;
    }
#elif defined (U_P_L)
    double* grd_p_l = val.GradU_mass ;   
#endif
    
    /*
     * w_tot     | K_tot_U   K_tot_C   K_tot_T  | GradU
     * w_salt  = | K_salt_U  K_salt_C  K_salt_T | GradC
     * w_the     | K_the_U   K_the_C   K_the_T  | GradT
     * 
     * If U = PL or U = Pmax(c > 0):
     * K_tot_U = kd_liq                   , K_tot_C = 0       , K_tot_T  = 0
     * 
     * If U = Pmax(c < 0):
     * K_tot_U = kd_liq*(V_Ice/V_H2O - 1) , K_tot_C = - R_g/V_H2O* dtlnadc_s , K_tot_T  = S_m/V_H2O - R_g/V_H2O* dtlnadtem
     * 
     * 
     * K_salt_U = mc_salt*kd_liq , K_salt_C = kf_salt , K_salt_T = 0
     * K_the_U  = s_l*kd_liq     , K_the_C  = 0       , K_the_T  = kth
     */
    {
      double rho_l = rho_h2o_l0 ;
      
      for(i = 0 ; i < 3 ; i++) {
        w_tot[i] = - kd_liq * grd_p_l[i] ;
      }
      w_tot[dim - 1] += kd_liq*rho_l*gravity ;
    }

    for(i = 0 ; i < 3 ; i++) {
      double w_l     = w_tot[i] ;
      double j_salt  = - kf_salt * grd_c_s[i] ;
      double q       = - kth * grd_tem[i] ;

      w_salt[i]  = mc_salt * w_l + j_salt ;
      w_the[i]   =     s_l * w_l + q ;
    }
  }
  
  /* Transfer coefficients at the current time */
  {
    double rho_l   = val.MassDensity_liquid ;
    double mu_l    = WaterViscosity(tem) ;
    double kr_l    = RelativePermeabilityToLiquid(p_c) ;
    double kd_liq  = rho_l*k_int/mu_l*kr_l ;
    double kf_salt = M_Salt*effectivediffusioncoef(phi0,sd_l) ;
    double lam_hom = thermalconductivity(phi0,sd_l) ;
    double kth     = lam_hom/tem ;
    double mc_salt = M_Salt * c_s / rho_l ;
      
    val.Permeability_liquid  = kd_liq ;
    val.Diffusion_salt       = kf_salt ;
    val.ThermalConductivity  = kth ;
    val.MassFraction_salt    = mc_salt ;
  }

  /* Body force */
  {
    double m_tot = val.Mass_total ;
    double* f_mass = val.BodyForce ;
    int i ;
      
    for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
    f_mass[dim - 1] = (rho_s + m_tot)*gravity ;
  }
  
  return(val) ;
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
