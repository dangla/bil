#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#ifdef HAVE_AUTODIFF
#define USE_AUTODIFF
#endif


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



#include "BaseName.h"
#include "CustomValues.h"
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

template<typename T>
using Values_t = CustomValues_t<T,ImplicitValues_t,ExplicitValues_t,ConstantValues_t,OtherValues_t> ;

using Values_d = Values_t<double> ;

#define Values_Index(V)  CustomValues_Index(Values_d,V,double)


#define MPM_t      BaseName(_MPM_t)




struct MPM_t: public MaterialPointModel_t<Values_t> {
  MaterialPointModel_SetInputs_t<Values_t> SetInputs;
  template<typename T>
  MaterialPointModel_Integrate_t<Values_t,T> Integrate;
  MaterialPointModel_Initialize_t<Values_t>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_t> SetTangentMatrix;
  MaterialPointModel_SetTransferMatrix_t<Values_t> SetTransferMatrix;
  MaterialPointModel_SetIndexes_t SetIndexes;
  MaterialPointModel_SetIncrements_t SetIncrements;
} ;







/* We define some names for implicit terms */
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
} ;



template<typename T = double>
struct OtherValues_t {
  T SaturationDegree_liquid;
  T Pressure_liquid;
  T Pressure_ice;
  T MassDensity_liquid;
  T MassDensity_ice;
  T ChemicalPotential_liquid;
  T ChemicalPotential_ice;
};



/* We define some names for explicit terms */
template<typename T = double>
struct ExplicitValues_t {
  T Permeability_liquid;
  T Diffusion_salt;
  T ThermalConductivity;
  T MassFraction_salt;
  T Entropy_liquid;
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
template<typename T = double>
struct ConstantValues_t {};



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
#define CapillaryPressureFactor(EPSV) \
        ScaleFactor(EPSV,aa,param.CapillaryPressureFraction)

#define YoungModulusFactor(EPSV) \
        ScaleFactor(EPSV,aa,param.YoungFraction)

#define PorosityFactor(EPSV) \
        (param.Thickness + Math_Max(EPSV,0))



/* The parameters below are read in the input data file */

#define Parameters_t    BaseName(_Parameters_t)

struct Parameters_t {
  double Porosity;
  double IntrinsicPermeability;
  double SolidVolumetricHeat;
  double SolidThermalConductivity;
  double SolidMatrixBulkModulus;
  double SolidMatrixShearModulus;
  double SolidMatrixThermalDilation;
  double ReferencePressure;
  double ReferenceTemperature;
  double YoungModulus;
  double PoissonRatio;
  double BiotCoefficient;
  double Thickness;
  double YoungFraction;
  double Gravity;
  double SolidSkeletonMassDensity;
  double CapillaryPressureFraction;
};

static Parameters_t param;
static double* cijkl ;
static Elasticity_t* elasty ;

template<typename T>
static T effectivediffusioncoef(double,T);
template static double  effectivediffusioncoef(double,double) ;
#ifdef USE_AUTODIFF
template static real  effectivediffusioncoef(double,real) ;
#endif
template<typename T>
T effectivediffusioncoef(double phi,T s)
{
  double h = param.Thickness;
  T tau_l_sat = 0.296e-3*exp(9.95*phi)/phi ;
  T tau = 0 ;
  
  if(s > 0.) {
    tau = tau_l_sat/(s*(1 + 625*pow(1 - s,4))) ;
  } else {
    tau = 0 ;
  }
  
  return(h*phi*s*tau*D_salt) ;
}


template<typename T>
static T thermalconductivity(double,T);
template static double thermalconductivity(double,double) ;
#ifdef USE_AUTODIFF
template static real thermalconductivity(double,real) ;
#endif
template<typename T>
T thermalconductivity(double phi,T sd_l)
{
  double lam_s = param.SolidThermalConductivity;
  double h = param.Thickness;
  T sd_i = 1 - sd_l ;
  T lam_h2o = sd_l*lam_l + sd_i*lam_i ;
  T lam_hom = lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
  
  return(h*lam_hom) ;
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
static void    GetProperties(Element_t*,double) ;
void GetProperties(Element_t* el,double t)
{
/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
  param.Gravity                    = GetProperty("gravity") ;
  param.SolidSkeletonMassDensity   = GetProperty("rho_s") ;
  param.Porosity                   = GetProperty("porosity") ;
  param.IntrinsicPermeability      = GetProperty("k_int") ;
  param.SolidVolumetricHeat        = GetProperty("C_s") ;
  param.SolidThermalConductivity   = GetProperty("lam_s") ;
  param.SolidMatrixThermalDilation = GetProperty("alpha_s") ;
  param.ReferencePressure          = GetProperty("p0") ;
  param.ReferenceTemperature       = GetProperty("T0") ;
  param.YoungModulus               = GetProperty("Young") ;
  param.PoissonRatio               = GetProperty("Poisson") ;
  param.BiotCoefficient            = GetProperty("Biot") ;
  
  param.Thickness = 1 ;
  param.YoungFraction  = 1 ;
  param.CapillaryPressureFraction = 1 ;
  
  if(Element_HasZeroThickness(el)) {
    param.Thickness                 = GetProperty("thickness") ;
    param.YoungFraction             = GetProperty("YoungFraction") ;
    param.CapillaryPressureFraction = 10 ;

    param.YoungModulus               /= param.Thickness ;
    param.IntrinsicPermeability      *= param.Thickness ;
    param.SolidVolumetricHeat        *= param.Thickness ;
    param.SolidMatrixThermalDilation *= param.Thickness ;
  }
  
  elasty  = Element_FindMaterialData(el,Elasticity_t,"Elasticity") ;
  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
#undef GetProperty
}



/* Functions used below */
static void    ComputePhysicoChemicalProperties(void) ;

template<typename T>
static T activity(T,T);
template static double  activity(double,double) ;
#ifdef USE_AUTODIFF
template static real  activity(real,real) ;
#endif

static double  activity_w_ideal(double,double) ;

template<typename T>
static T lna_i(T,T,double,double,double,T);
template static double  lna_i(double,double,double,double,double,double) ;
#ifdef USE_AUTODIFF
template static real  lna_i(real,real,double,double,double,real) ;
#endif



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
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
  
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
    double phi0  = Material_GetPropertyValue(mat,"porosity") ;
    
    if(k_s > 0 && g_s > 0) {
      /* Elastic moduli (Mori-Tanaka) */
      double k_0 = k_s ;
      double g_0 = g_s ;
      double a_0 = 0.75 * k_0 / g_0 ;
      double b_0 = 6 * (k_0 + 2 * g_0) / (9 * k_0 + 8 * g_0) ;
      double K = (1 - phi0) * k_s / (1 + phi0 * a_0 * k_s / k_0) ;
      double G = (1 - phi0) * g_s / (1 + phi0 * b_0 * g_s / g_0) ;
      double young =  9 * G * K / (3 * K + G) ;
      double poisson = 0.5 * (1 - young / (3 * K)) ;
      //double poisson = (1.5 * K - G) / (3 * K + G) ;
      double biot = 1 - K/k_s ;
      
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
      double young =  Material_GetProperty(mat)[pm("Young")] ;
      double poisson = Material_GetProperty(mat)[pm("Poisson")] ;
    
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
  
  MaterialPointModel_DefineNbOfInternalValues(MPM_t,el,NbOfIntPoints);

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
  int i = MaterialPointModel_ComputeInitialStateByFEM(MPM_t,el,t);
  
  return(i);
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the (current) explicit terms.
 *  IMPORTANT: if needed use only the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  int i = MaterialPointModel_ComputeExplicitTermsByFEM(MPM_t,el,t);
  
  return(i);
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the (current) implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  int i = MaterialPointModel_ComputeImplicitTermsByFEM(MPM_t,el,t,dt);
  
  return(i);
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
  int i = MaterialPointModel_ComputePoromechanicalMatrixByFEM(MPM_t,el,t,dt,k,E_MECH);
  
  return(i);
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
  /* Initialization */
  {
    int ndof = Element_GetNbOfDOF(el) ;
    
    for(int i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }
  /* 1. Mechanics */
  if(Element_EquationIsActive(el,E_MECH)) {
    MaterialPointModel_ComputeMechanicalEquilibriumResiduByFEM(MPM_t,el,t,dt,r,E_MECH,Stress,BodyForce);
  }
  /* 2. Conservation of total mass */
  if(Element_EquationIsActive(el,E_MASS)) {
    MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_MASS,Mass_total,MassFlow_total);
  }
  /* 3. Conservation of salt */
  if(Element_EquationIsActive(el,E_SALT)) {  
    MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_SALT,Mass_salt,MassFlow_salt);
  }
  /* 4. Entropy balance */
  if(Element_EquationIsActive(el,E_THER)) {  
    MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_THER,Entropy_total,HeatFlow);
  }
  
  return(0);
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  int NbOfOutputs  = 17 ;
  double* vi  = Element_GetCurrentImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n  = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /* Initialization */
  {
    for(int i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }
  
  /*
    Input data
  */
  Element_ComputeMaterialProperties(el,t) ;

  
  //ci.Set(el,t,0,u,vi,u,vi) ;
  
  
  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    //Values_d& val = *ci.ExtractInputs(p) ;
    /* displacements */
    double* disp = Element_ComputeDisplacementVector(el,u,intfct,p,U_MECH) ;
    /* strains */
    double* eps =  Element_ComputeLinearStrainTensor(el,u,intfct,p,U_MECH) ;
    double* eps_n =  Element_ComputeLinearStrainTensor(el,u_n,intfct,p,U_MECH) ;
    double  epsv = eps[0] + eps[4] + eps[8] ;
    double  epsv_n = eps_n[0] + eps_n[4] + eps_n[8] ;
    /* temperature */
    double tem = Element_ComputeUnknown(el,u,intfct,p,U_THER) ;
    double dtem = tem - param.ReferenceTemperature ;
    /* salt concentration */
    double c_salt = Element_ComputeUnknown(el,u,intfct,p,U_SALT) ;
    /* Log of water activity */
    double lna = LogActivityOfWater(c_salt,tem) ;
    /* pressures */
    double u_mass = Element_ComputeUnknown(el,u,intfct,p,U_MASS) ;
    #if defined (U_P_MAX)
    double  p_max = u_mass ;
    //     0 = V_Ice*(p_i   - p_m) - V_H2O*(p_l   - p_m) + S_m*(tem - T_m) - R_g*tem*lna
    double c = V_Ice*(p_max - p_m) - V_H2O*(p_max - p_m) + S_m*(tem - T_m) - R_g*tem*lna ;
    double p_min = (c > 0) ? p_max - c/V_Ice : p_max + c/V_H2O ;
    double p_i   = (c > 0) ? p_min : p_max ;
    double p_l   = (c > 0) ? p_max : p_min ;
    #elif defined (U_P_L)
    double p_l   = u_mass ;
    double p_i   = p_m + (V_H2O*(p_l - p_m) - S_m*(tem - T_m) + R_g*tem*lna)/V_Ice ;
    #endif
    double dp_l = p_l - param.ReferencePressure ;
    double dp_i = p_i - param.ReferencePressure ;
    double p_c  = p_i - p_l ;
    /* fluxes */
    double w_tot[3]  = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double w_the[3]  = {0,0,0} ;
    double young1 = param.YoungModulus*YoungModulusFactor(epsv_n) ;
    double K = young1 / (3 - 6*param.PoissonRatio) ;
    double G = 0.5 * young1 / (1 + param.PoissonRatio) ;
    double b = param.BiotCoefficient ;
    double sd_l = SaturationDegree(p_c*CapillaryPressureFactor(epsv)) ;
    double sd_i = 1 - sd_l ;
    double b_i  = b*sd_i ;
    double b_l  = b*sd_l ;
    double N_il = (K > 0) ? - (b - param.Porosity) / K * b * sd_i * sd_l  : 0 ;
    double N_ii = (K > 0) ?   (b - param.Porosity) / K * sd_i * (1 - b_i) : 0 ;
    double N_ll = (K > 0) ?   (b - param.Porosity) / K * sd_l * (1 - b_l) : 0 ;
    double alpha_phi = (b - param.Porosity)*param.SolidMatrixThermalDilation ;
    double alpha_phi_l = sd_l*alpha_phi ;
    double alpha_phi_i = sd_i*alpha_phi ;
    
    double k_oedo = (K + 4/3.*G) ; /* k_oedo Oedometric modulus */
    double Kk_oedo = (1 + param.PoissonRatio)/(3 - 3*param.PoissonRatio) ; /* K/k_oedo */
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    CustomValues_t<double,ImplicitValues_t>* val1 = (CustomValues_t<double,ImplicitValues_t>*) vi ;
    int i;
    
    /* Averaging */
    for(int i = 0 ; i < np ; i++) {
      for(int j = 0 ; j < 3 ; j++) w_tot[j]  += val1[i].MassFlow_total[j]/np ;
      for(int j = 0 ; j < 3 ; j++) w_salt[j] += val1[i].MassFlow_salt[j]/np ;
      for(int j = 0 ; j < 3 ; j++) w_the[j]  += val1[i].HeatFlow[j]/np ;
      for(int j = 0 ; j < 9 ; j++) sig[j]    += val1[i].Stress[j]/np ;
    }
    
    i = 0  ;
    /* quantites exploitees */
    Result_Store(r + i++,&p_l,"Liquid pressure",1) ;
    Result_Store(r + i++,&c_salt,"Salt concentration",1) ;
    Result_Store(r + i++,&tem,"Temperature",1) ;
    Result_Store(r + i++,&sd_l,"Liquid saturation",1) ;
    Result_Store(r + i++,w_tot,"Liquid mass flux vector",3) ;
    Result_Store(r + i++,w_salt,"Salt mass flux vector",3) ;
    Result_Store(r + i++,&p_i,"Ice pressure",1) ;
    /* Mass densities */
    {
      double alpha_l   = ALPHA_l(tem - T_m) ;
      double rho_salt  = M_Salt*c_salt ;
      double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - vr_salt*rho_salt) ;
      
      Result_Store(r + i++,&rho_l,"Solution density",1) ;
    }
    {
      double rho_i     = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
      
      Result_Store(r + i++,&rho_i,"Ice density",1) ;
    }
    Result_Store(r + i++,eps,"Strain tensor",9) ;
    Result_Store(r + i++,sig,"Stress tensor",9) ;
    Result_Store(r + i++,disp,"Displacement vector",3) ;
      
    {
      double Eps_T  = param.SolidMatrixThermalDilation*Kk_oedo*dtem ;
      
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
      double lnaw = LogActivityOfWater(c_salt,tem) ;
      
      Result_Store(r + i++,&lnaw,"Lna_w",1) ;
    }
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}





int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
{
  int   ncols = 9 + 3;
  int   dec = ncols*ncols ;
  double* c0 = c + p*dec ;

    /* Derivatives w.r.t strains */
    #if 0
    if(k >= 0 && k < 9) {
      /* Tangent stiffness matrix */
      {
        double* c1 = c0 ;
         
        /* Already updated in Integrate */
        #define C1(i,j)     ((c1)[(i)*9+(j)])
        #define Cijkl(i,j)  ((cijkl)[(i)*9+(j)])
        for(int i = 0 ; i < 9 ; i++) {
          C1(i,k) = Cijkl(i,k) ;
        }
        #undef C1
        #undef Cijkl

        //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
      }
        
      /* Coupling matrices */
      {
        double* c00 = c0 + 9*ncols ;
          
        /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
          * and zero for the derivatives w.r.t others */
        {
          double* c1 = c00 + ncols*U_MASS ;
          
          if(k == 0 || k == 4 || k == 8) c1[k] = dval.Mass_total ;
        }
        {
          double* c1 = c00 + ncols*U_SALT ;

          if(k == 0 || k == 4 || k == 8) c1[k] = dval.Mass_salt ;
        }
        {
          double* c1 = c00 + ncols*U_THER ;

          if(k == 0 || k == 4 || k == 8) c1[k] = dval.Entropy_total ;
        }
      }
      return(dec) ;
    }
    #else
    if(k == 0) {
      /* Tangent stiffness matrix */
      {
        double* c1 = c0 ;
         
        /* Already updated in Integrate */
        for(int i = 0 ; i < 81 ; i++) {
          c1[i] = cijkl[i] ;
        }

        //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
      }
        
      /* Coupling matrices */
      {
        double* c00 = c0 + 9*ncols ;
          
        /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
          * and zero for the derivatives w.r.t others */
        #define B1(i,j)  ((c1)[(i)*3+(j)])
        {
          double* c1 = c00 + ncols*U_MASS ;
          
          for(int i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass_total ;
        }
        {
          double* c1 = c00 + ncols*U_SALT ;

          for(int i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass_salt ;
        }
        {
          double* c1 = c00 + ncols*U_THER ;

          for(int i = 0 ; i < 3 ; i++) B1(i,i) = dval.Entropy_total ;
        }
        #undef B1
      }
      return(dec) ;
    }
    #endif
      
    /* Derivatives w.r.t U_MASS/U_SALT/U_THER */
    if(k >= 9) {
      {
        /* Tangent Biot's coefficient */
        {
          double* c1 = c0 + 9*k ;

          for(int i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
        }
      }
        
      /* General storage matrix */
      {
        double* c00 = c0 + 9*ncols + k ;
          
        {
          double* c1 = c00 + ncols*U_MASS ;
        
          c1[0] = dval.Mass_total ;
        }
        {
          double* c1 = c00 + ncols*U_SALT ;
        
          c1[0] = dval.Mass_salt ;
        }
        {
          double* c1 = c00 + ncols*U_THER ;
        
          c1[0] = dval.Entropy_total ;
        }
      }
      return(dec) ;
    }

  return(dec) ;
}

  


void MPM_t::SetIndexes(Element_t* el,int* ind) {
    ind[0] = Values_Index(Strain[0]);
    
    for(int i = 1 ; i < 9 ; i++) {
      ind[i] = ind[0] + i;
    }
    
    ind[9]  = Values_Index(U_mass);  
    ind[10] = Values_Index(Concentration_salt);
    ind[11] = Values_Index(Temperature);
}




void MPM_t::SetIncrements(Element_t* el,double* dui) {  
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
        
    for(int i = 0 ; i < 9 ; i++) {
      dui[i] =  1.e-6 ;
    }
    
    dui[9]  = 1.e-2*ObVal_GetValue(obval + U_MASS) ;
    dui[10] = 1.e-2*ObVal_GetValue(obval + U_SALT) ;
    dui[11] = 1.e-2*ObVal_GetValue(obval + U_THER) ;
}




int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
{
  int    dec = 9 * 9 ;

  {
    double* c0 = c + p*dec ;
    
    {
      /* Transfer coefficients for the liquid mass flux */
      {
        double* c1 = c0 ;
        
        c1[0] = dt*val.Permeability_liquid ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      
      /* Transfer coefficients for the salt mass flux */
      {
        double* c1 = c0 + 3 * 9 ;
        
        c1[0] = dt*val.MassFraction_salt*val.Permeability_liquid ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      {
        double* c1 = c0 + 4 * 9 ;
        
        c1[0] = dt*val.Diffusion_salt ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      
      /* Transfer coefficients for the entropy flux */
      {
        double* c1 = c0 + 6 * 9 ;
        
        c1[0] = dt*val.Entropy_liquid*val.Permeability_liquid ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      {
        double* c1 = c0 + 8 * 9 ;
        
        c1[0] = dt*val.ThermalConductivity ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
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
    
  /* Pressure and pressure gradient */
  var.ValueAndGradientFEM(el,p,U_MASS,&val.U_mass) ;
  
  /* Salt concentration and concentration gradient */
  var.ValueAndGradientFEM(el,p,U_SALT,&val.Concentration_salt) ;
  
  /* Temperature and temperature gradient */
  var.ValueAndGradientFEM(el,p,U_THER,&val.Temperature) ;
  
  return(&val) ;
}



template <typename T>
Values_t<T>* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_t<T>& val)
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Retrieve the primary variables */
  /* Strains */
  T* eps   = val.Strain ;
  double const* eps_n = val_n.Strain ;
  T  epsv  = eps[0] + eps[4] + eps[8] ;
  double  epsv_n = eps_n[0] + eps_n[4] + eps_n[8] ;
  /* Salt concentration */
  T c_s    = val.Concentration_salt ;
  /* Temperature */
  T tem    = val.Temperature ;
  /* Log of water activity */
  T lna = LogActivityOfWater(c_s,tem) ;
  
  /* Pressures */
#if defined (U_P_MAX)
  T  p_max = val.U_mass ;
  //     0 = V_Ice*(p_i   - p_m) - V_H2O*(p_l   - p_m) + S_m*(tem - T_m) - R_g*tem*lna
  T c = V_Ice*(p_max - p_m) - V_H2O*(p_max - p_m) + S_m*(tem - T_m) - R_g*tem*lna ;
  T p_min = (c > 0) ? p_max - c/V_Ice : p_max + c/V_H2O ;
  T p_i   = (c > 0) ? p_min : p_max ;
  T p_l   = (c > 0) ? p_max : p_min ;
#elif defined (U_P_L)
  T p_l   = val.U_mass ;
  T p_i   = p_m + (V_H2O*(p_l - p_m) - S_m*(tem - T_m) + R_g*tem*lna)/V_Ice ;
#endif
  T p_c = p_i - p_l ;
  
  /* Porosity */
  T phi1 = param.Porosity*PorosityFactor(epsv) ;
  
  /* Saturations */
  T sd_l = SaturationDegree(p_c*CapillaryPressureFactor(epsv)) ;
  T sd_i = 1 - sd_l ;
  
  /* Backup */
    
  val.SaturationDegree_liquid = sd_l ;
  val.Pressure_liquid = p_l ;
  val.Pressure_ice = p_i ;
  

  /* Mass contents and entropies */
  {
    /* Mass densities */
    T alpha_l   = ALPHA_l(tem - T_m) ;
    T rho_i     = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    T rho_salt  = M_Salt*c_s ;
    T rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - vr_salt*rho_salt) ;
  
    /* Mass contents */
    T m_l     = phi1*sd_l*rho_l ;
    T m_i     = phi1*sd_i*rho_i ;
    T m_salt  = phi1*sd_l*rho_salt ;
  
    /* Entropies */
    T s_l   = C_l*log(tem/T_m) - alpha_l/rho_h2o_l0*(p_l - p_m) ;
    T s_i   = C_i*log(tem/T_m) - alpha_i/rho_h2o_i0*(p_i - p_m) - S_m/M_H2O ;
    T s_sol = param.SolidVolumetricHeat*log(tem/T_m) ;
    T s_tot = s_sol + m_l*s_l + m_i*s_i ;
  
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
    T mu_l = V_H2O*(p_l - p_m) + R_g*tem*lna ;
    T mu_i = V_Ice*(p_i - p_m) + S_m*(tem - T_m) ;
    
    val.ChemicalPotential_liquid = mu_l ;
    val.ChemicalPotential_ice = mu_i ;
  }
  
  /* Stresses */
  {
    double young1 = param.YoungModulus*YoungModulusFactor(epsv_n) ;
    double K = young1 / (3 - 6*param.PoissonRatio) ;
    double b    = param.BiotCoefficient ;
    T b_i  = b*sd_i ;
    T b_l  = b*sd_l ;
    
    T dp_l = p_l - param.ReferencePressure ;
    T dp_i = p_i - param.ReferencePressure ;
    T dtem = tem - param.ReferenceTemperature ;
    
    T* sig = val.Stress ;
          
    Elasticity_SetParameters(elasty,young1,param.PoissonRatio) ;
    Elasticity_UpdateStiffnessTensor(elasty) ;
    
    {    
      /* Initialize stresses */
      for(int i = 0 ; i < 9 ; i++) sig[i] = 0 ;
      
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(int i = 0 ; i < 9 ; i++) {
        for(int j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*eps[j] ;
        }
      }
      #undef C
    }
    
    {
      sig[0] += - b_l * dp_l - b_i * dp_i - param.SolidMatrixThermalDilation * K * dtem ;
      sig[4] += - b_l * dp_l - b_i * dp_i - param.SolidMatrixThermalDilation * K * dtem ;
      sig[8] += - b_l * dp_l - b_i * dp_i - param.SolidMatrixThermalDilation * K * dtem ;
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
    T* w_tot   = val.MassFlow_total ;
    T* w_salt  = val.MassFlow_salt ;
    T* w_the   = val.HeatFlow ; 
    
    /* Gradients */
    T* grd_c_s   = val.GradConcentration_salt ;
    T* grd_tem   = val.GradTemperature ;
    
#if defined (U_P_MAX)
    T  grd_p_l[3] ;
    T* grd_p_max = val.GradU_mass ;  
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    double dc_s = 1.e-2 * ObVal_GetValue(obval + U_SALT) ;
    T dtlnadc_s = (tem * LogActivityOfWater(c_s + dc_s,tem) - tem * lna) / dc_s ;
    double dtem = 1.e-2 * ObVal_GetValue(obval + U_THER) ;
    T dtlnadtem = ((tem + dtem) * LogActivityOfWater(c_s,tem + dtem) - tem * lna) / dtem ;
    
    for(int i = 0 ; i < 3 ; i++) {
      T grd_tlna = dtlnadc_s * grd_c_s[i] + dtlnadtem * grd_tem[i] ;
      T grd_c    = (V_Ice - V_H2O) * grd_p_max[i] + S_m * grd_tem[i] - R_g * grd_tlna ;
      T grd_p_min = (c > 0) ? grd_p_max[i] - grd_c/V_Ice : grd_p_max[i] + grd_c/V_H2O ;
      
      grd_p_l[i]   = (c > 0) ? grd_p_max[i] : grd_p_min ;
    }
#elif defined (U_P_L)
    T* grd_p_l = val.GradU_mass ;   
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
      
      for(int i = 0 ; i < 3 ; i++) {
        w_tot[i] = - kd_liq * grd_p_l[i] ;
      }
      w_tot[dim - 1] += kd_liq*rho_l*param.Gravity  ;
    }

    for(int i = 0 ; i < 3 ; i++) {
      T w_l     = w_tot[i] ;
      T j_salt  = - kf_salt * grd_c_s[i] ;
      T q       = - kth * grd_tem[i] ;

      w_salt[i]  = mc_salt * w_l + j_salt ;
      w_the[i]   =     s_l * w_l + q ;
    }
  }
  
  /* Transfer coefficients at the current time */
  {
    T rho_l   = val.MassDensity_liquid ;
    T mu_l    = WaterViscosity(tem) ;
    T kr_l    = RelativePermeabilityToLiquid(p_c) ;
    T kd_liq  = rho_l*param.IntrinsicPermeability/mu_l*kr_l ;
    T kf_salt = M_Salt*effectivediffusioncoef(param.Porosity,sd_l) ;
    T lam_hom = thermalconductivity(param.Porosity,sd_l) ;
    T kth     = lam_hom/tem ;
    T mc_salt = M_Salt * c_s / rho_l ;
      
    val.Permeability_liquid  = kd_liq ;
    val.Diffusion_salt       = kf_salt ;
    val.ThermalConductivity  = kth ;
    val.MassFraction_salt    = mc_salt ;
  }

  /* Body force */
  {
    T m_tot = val.Mass_total ;
    T* f_mass = val.BodyForce ;
      
    for(int i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
    f_mass[dim - 1] = (param.SolidSkeletonMassDensity + m_tot)*param.Gravity  ;
  }
  
  return(&val) ;
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  return(NULL);
}



template<typename T>
T activity(T c_s1,T tem1)
/* activity of water */
{
  T c_s = (c_s1 > 0) ? c_s1 : 0 ;
  T tem = (tem1 > 200) ? tem1 : 200 ;
  T T_0  = T_m ;
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

  T epsi = 0.0007*(tem - T_0)*(tem - T_0) - 0.3918*(tem - T_0) + 87.663 ;
  T A    = 1398779.816/pow(epsi*tem,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;
  T c_ani ;
  T c_cat ;
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
    T c_h2o = (1 - c_s*(NU_A*V_A + NU_C*V_C))/V_H2O ;
    /* molalites * M_H2O */
    T m_ani  = c_ani/c_h2o ;
    T m_cat  = c_cat/c_h2o ;

    /* ionic strength */
    T I     =  0.5*(z_ani*z_ani*m_ani + z_cat*z_cat*m_cat);
  
    T II_ani   = lna_i(tem,I,z_ani,b_ani,S_ani,A) ;
    T II_cat   = lna_i(tem,I,z_cat,b_cat,S_cat,A) ;

    /* activity of water */
    T lna_h2o = m_ani*II_ani + m_cat*II_cat ;
    /* linearized activity of water */
    //T lna_h2o = - (m_ani + m_cat) ;

    return(lna_h2o) ;
  }
}


template<typename T>
T lna_i(T tem,T I,double z,double b,double S,T A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29 ;
  double a1 = alpha/(1+alpha) ;
  T II = sqrt(I) ;
  T lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/tem ;
  T lna_i = -1 + lna*z*z ;
  
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
