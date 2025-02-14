#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "Barcelona Basic Model for unsaturated soils with gas (2025)"
#define AUTHORS "Eizaguirre-Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (2+dim)

/* Indices of equation */
#define E_MASS  (0)
#define E_AIR   (1)
#define E_MECH  (2)


/* Generic indices of nodal unknowns */
#define U_MASS   E_MASS
#define U_AIR    E_AIR
#define U_MECH   E_MECH

/* Method chosen at compiling time.
 * Each equation is associated with a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */
/* Indices of unknowns */
#define U_P_L   U_MASS
#define U_P_G   U_AIR
#define U_DIS   U_MECH


#include "BaseName.h"
#include "CustomValues.h"
#include "MaterialPointModel.h"


#define ImplicitValues_t BaseName(_ImplicitValues_t)
#define ExplicitValues_t BaseName(_ExplicitValues_t)
#define ConstantValues_t BaseName(_ConstantValues_t)
#define OtherValues_t    BaseName(_OtherValues_t)




template<typename T>
struct ImplicitValues_t;

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



template<typename T = double>
struct ImplicitValues_t {
  T Displacement[3];
  T Strain[9];
  T Pressure_liquid;
  T GradPressure_liquid[3];
  T Pressure_gas;
  T GradPressure_gas[3];
  T Mass_total;
  T MassFlow_total[3];
  T Mass_air;
  T MassFlow_air[3];
  T MassFlow_gas[3];
  T Stress[9];
  T BodyForce[3];
  T PlasticStrain[9];
  T HardeningVariable;
  T YieldCriterion;
  T PlasticMultiplier;
  T Porosity;
};



template<typename T = double>
struct OtherValues_t {
  T SaturationDegree_liquid;
  T MassDensity_liquid;
  T MassDensity_gas;
  T MassDensity_air;
};



/* We define some names for explicit terms */
template<typename T = double>
struct ExplicitValues_t {
  T TransportMass_liquid;
  T TransportMass_gas;
  T TransportAir_liquid;
  T TransportAir_gas;
  T MassFraction_air;
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
template<typename T = double>
struct ConstantValues_t {};




/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*,double) ;
static void   ComputePhysicoChemicalProperties(double) ;
static double saturationdegree(double,double,Curve_t*) ;



#define ComputeTangentStiffnessTensor(...)  Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...)                  Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)              Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define UpdateElastoplasticTensor(...)      Plasticity_UpdateElastoplasticTensor(plasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...)     Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)



/* Units
 * ----- */
#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define meter (InternationalSystemOfUnits_OneMeter)
#define dm    (0.1*meter)
#define cm    (0.01*meter)
#define dm2   (dm*dm)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)
#define Pascal   (InternationalSystemOfUnits_OnePascal)
#define MPa   (1.e6*Pascal)
#define GPa   (1.e3*MPa)
#define mol   InternationalSystemOfUnits_OneMole
#define sec   InternationalSystemOfUnits_OneSecond
#define kg    InternationalSystemOfUnits_OneKilogram
#define gr    (0.001*kg)


#define TEMPERATURE  (293)


/* Material properties */
//#define SaturationDegree(pc)    (Curve_ComputeValue(saturationcurve,pc))
#define SaturationDegree(pc)    (saturationdegree(pc,p_c3,saturationcurve))

#define RelativePermeabilityToLiquid(pc)  (Curve_ComputeValue(relativepermliqcurve,pc))
#define RelativePermeabilityToGas(pc)     (Curve_ComputeValue(relativepermgascurve,pc))

#define TortuosityToGas(f,sg)     (0.1) //((sg > 0) ? pow(f,aa)*pow(sg,bb) : 0)
//#define TortuosityToGas(f,sg)   ((sg > 0) ? pow(f*sg,4./3) : 0) // After HYDRUS
#define aa                        (0.33)  /* 1/3 Millington, Thiery 1.74 */
#define bb                        (2.33)   /* 7/3 Millington, Thiery 3.2 */

#define KappaSuction(p)    ((kappasuctioncurve) ? Curve_ComputeValue(kappasuctioncurve,p) : kappa_s)
#define Kappa(s)           ((kappacurve) ? Curve_ComputeValue(kappacurve,s) : kappa)



#include "MolarMassOfMolecule.h"

/* Water property
 * -------------- */
 /* Molar mass */
#define M_H2O          MolarMassOfMolecule(H2O)
/* Molar volume of liquid water */
#define V_H2O          (18 * cm3)
/* Mass density */
#define MassDensityOfWaterVapor(pv)    (M_H2O*(pv)/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(pc)           (exp(-V_H2O/RT*(pc)))
#define VaporPressure(pc)              (p_v0*RelativeHumidity(pc))
#define dRelativeHumidity(pc)          (-V_H2O/RT*RelativeHumidity(pc))
#define dVaporPressure(pc)             (p_v0*dRelativeHumidity(pc))


/* Gas property
 * -------------- */
 /* Molar mass of the inert gas */
#define M_AIR          (28.97*gr)



/* Parameters */
static double  e0 ;
static double  phi0 ;
static double  hardv0 ;
static Elasticity_t* elasty ;
static Plasticity_t* plasty ;
static Curve_t* saturationcurve ;
static Curve_t* relativepermliqcurve ;
static Curve_t* relativepermgascurve ;
static Curve_t* kappasuctioncurve ;
static Curve_t* kappacurve ;
static double  kappa ;
static double  kappa_s ;
static double  p_atm ;
static double  p_v0 ;
static double  RT ;
static double  p_c3 ;

/* The parameters below are read in the input data file */

#define Parameters_t    BaseName(_Parameters_t)

struct Parameters_t {
  double Gravity;
  double MassDensity_solidskeleton;
  double InitialStress[9];
  double MassDensity_liquid;
  double IntrinsicPermeability_liquid;
  double IntrinsicPermeability_gas;
  double Viscosity_liquid;
  double Viscosity_gas;
  double PoissonRatio;
  double InitialVoidRatio;
  double InitialPorosity;
  double InitialHardeningVariable;
  double DiffusionCoefficient_watervapor;
  double DiffusionCoefficient_dissolvedair;
  double SlopeOfElasticLoadingLine;
  double SlopeOfVirginConsolidationLine;
  double SlopeOfCriticalStateLine;
  double InitialPreconsolidationPressure;
  double SlopeOfElasticWettingDryingLine;
  double HenryConstant;
  double ShearModulus;
  double SuctionCohesionCoefficient;
  double ReferenceConsolidationPressure;
  double MinimumSuction;
};


#include "PhysicalConstant.h"
//#include "AtmosphericPressure.h"
//#include "WaterViscosity.h"
//#include "AirViscosity.h"
#include "WaterVaporPressure.h"
//#include "DiffusionCoefficientOfMoleculeInAir.h"

void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Air (dm2/s) */
  
  /* Viscosities */
  //mu_l    = WaterViscosity(TK) ;
  //mu_g    = AirViscosity(TK) ;
  
  /* Water vapor pressure */
  p_v0    = WaterVaporPressure(TK) ;
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
  p_atm = 101325.*Pascal ;
  
  /* Liquid mass density */
  //rho_l0 = 1 * kg/dm3 ;
}





int pm(const char *s)
{
#define Parameters_Index(V)  CustomValues_Index(Parameters_t,V,double)
         if(!strcmp(s,"gravity"))    { 
    return (Parameters_Index(Gravity)) ;
  } else if(!strcmp(s,"rho_s"))      { 
    return (Parameters_Index(MassDensity_solidskeleton)) ;
  } else if(!strcmp(s,"shear_modulus")) { 
    return (Parameters_Index(ShearModulus)) ;
  } else if(!strcmp(s,"rho_l"))      { 
    return (Parameters_Index(MassDensity_liquid)) ;
  } else if(!strcmp(s,"kl_int"))      { 
    return (Parameters_Index(IntrinsicPermeability_liquid)) ;
  } else if(!strcmp(s,"mu_l"))       { 
    return (Parameters_Index(Viscosity_liquid)) ;
  } else if(!strcmp(s,"initial_stress"))       {
    return(Parameters_Index(InitialStress[0])) ;
  } else if(!strncmp(s,"initial_stress_",15))   {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(Parameters_Index(InitialStress[3*i + j])) ;
    
    /* BBM */
  } else if(!strcmp(s,"slope_of_swelling_line")) {
    return(Parameters_Index(SlopeOfElasticLoadingLine)) ;
  } else if(!strcmp(s,"slope_of_virgin_consolidation_line")) {
    return(Parameters_Index(SlopeOfVirginConsolidationLine)) ;
  } else if(!strcmp(s,"slope_of_critical_state_line"))  {
    return(Parameters_Index(SlopeOfCriticalStateLine)) ;
  } else if(!strcmp(s,"initial_pre-consolidation_pressure")) {
    return(Parameters_Index(InitialPreconsolidationPressure)) ;
  } else if(!strcmp(s,"initial_porosity")) {
    return(Parameters_Index(InitialPorosity)) ;
  } else if(!strcmp(s,"kappa_s")) {
    return(Parameters_Index(SlopeOfElasticWettingDryingLine)) ;
  } else if(!strcmp(s,"suction_cohesion_coefficient")) {
    return(Parameters_Index(SuctionCohesionCoefficient)) ;
  } else if(!strcmp(s,"reference_consolidation_pressure")) {
    return(Parameters_Index(ReferenceConsolidationPressure)) ;
  } else if(!strcmp(s,"kg_int"))      { 
    return (Parameters_Index(IntrinsicPermeability_gas)) ;
  } else if(!strcmp(s,"mu_g"))       { 
    return (Parameters_Index(Viscosity_gas)) ;
  } else if(!strcmp(s,"vapor_diffusion_coefficient"))       { 
    return (Parameters_Index(DiffusionCoefficient_watervapor)) ;
  } else if(!strcmp(s,"poisson"))       { 
    return (Parameters_Index(PoissonRatio)) ;
  } else if(!strcmp(s,"minimum_suction"))  {
    return (Parameters_Index(MinimumSuction)) ;
  } else if(!strcmp(s,"henry_law_constant"))  {
    return (Parameters_Index(HenryConstant)) ;
  } else if(!strcmp(s,"dissolved_air_diffusion_coefficient"))       { 
    return (Parameters_Index(DiffusionCoefficient_dissolvedair)) ;
  } else return(-1) ;
#undef Parameters_Index
}


void GetProperties(Element_t* el,double t)
{
/* To retrieve the material properties */
  Parameters_t& param = ((Parameters_t*) Element_GetProperty(el))[0] ;
  
  phi0    = param.InitialPorosity;
  e0      = phi0/(1 - phi0) ;
  kappa_s = param.SlopeOfElasticWettingDryingLine;
  p_c3    = param.MinimumSuction;
  
  plasty  = (Plasticity_t*) Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  elasty  = Plasticity_GetElasticity(plasty) ;
  
  hardv0  = Plasticity_GetHardeningVariable(plasty)[0] ;
  
  saturationcurve = Element_FindCurve(el,"sl") ;
  relativepermliqcurve = Element_FindCurve(el,"kl") ;
  relativepermgascurve = Element_FindCurve(el,"kg") ;
  kappasuctioncurve = Element_FindCurve(el,"kappa_s") ;
  kappacurve = Element_FindCurve(el,"kappa") ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
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
  Model_CopyNameOfEquation(model,E_MASS,"mass") ;
  Model_CopyNameOfEquation(model,E_AIR,"air") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  Model_CopyNameOfUnknown(model,U_P_G,"p_g") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_DIS + i,name_unk[i]) ;
  }
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
    
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = ((int) sizeof(Parameters_t)/sizeof(double)) ;
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
      /* Barcelona Basic model */
      {
        double lambda = Material_GetPropertyValue(mat,"slope_of_virgin_consolidation_line") ;
        double M      = Material_GetPropertyValue(mat,"slope_of_critical_state_line") ;
        double pc0    = Material_GetPropertyValue(mat,"initial_pre-consolidation_pressure") ;
        double coh    = Material_GetPropertyValue(mat,"suction_cohesion_coefficient") ;
        double p_ref  = Material_GetPropertyValue(mat,"reference_consolidation_pressure") ;
        Curve_t* lc   = Material_FindCurve(mat,"lc") ;
        
        kappa   = Material_GetPropertyValue(mat,"slope_of_swelling_line") ;
        phi0    = Material_GetPropertyValue(mat,"initial_porosity") ;
        e0      = phi0/(1 - phi0) ;
        
        Plasticity_SetTo(plasty,BBM) ;
        Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0,coh,p_ref,lc) ;
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
\t 1. The mass balance equation for the total mass (mass)\n\
\t 2. The mass balance equation for the air mass (air)\n\
\t 3. The equilibrium equation (meca_1,meca_2,meca_3)\n") ;

  printf("\n\
The primary unknowns are\n\
\t 1. The liquid pressure (p_l)\n\
\t 2. The gas pressure (p_g)\n\
\t 3. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # mass density of solid skeleton\n") ;
  fprintf(ficd,"young = 5.8e+09   # Young's modulus\n") ;
  fprintf(ficd,"poisson = 0.3     # Poisson's ratio\n") ;
  fprintf(ficd,"porosity = 0.15   # porosity\n") ;
  fprintf(ficd,"rho_l = 1000      # mass density of fluid\n") ;
  fprintf(ficd,"p_g = 0           # gas pressure\n") ;
  fprintf(ficd,"kl_int = 1e-19     # intrinsic permeability\n") ;
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
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;
  
  MaterialPointModel_DefineNbOfInternalValues(MPM_t,el,NbOfIntPoints);
  
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_MASS) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
      
    /* other */
    } else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
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
  {
    MaterialPointModel_ComputeMechanicalEquilibriumResiduByFEM(MPM_t,el,t,dt,r,E_MECH,Stress,BodyForce);
  }
  /* 2. Conservation of total mass */
  {
    MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_MASS,Mass_total,MassFlow_total);
  }
  /* 3. Conservation of air */
  {  
    MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_AIR,Mass_air,MassFlow_air);
  }
  
  return(0);
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 12 ;
  double* vex0  = Element_GetExplicitTerm(el) ;
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Element_GetDimensionOfSpace(el) ;

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
  GetProperties(el,t) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Pressure */
    double p_l = Element_ComputeUnknown(el,u,intfct,p,U_P_L) ;
    double p_g = Element_ComputeUnknown(el,u,intfct,p,U_P_G) ;
    double pc  = p_g - p_l ;
    /* saturation */
    double sl = SaturationDegree(pc) ;
    /* Displacement */
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps[9] = {0,0,0,0,0,0,0,0,0} ;
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double tre,e ;
    double w_tot[3] = {0,0,0} ;
    double w_air[3] = {0,0,0} ;
    double w_gas[3] = {0,0,0} ;
    double j_vap[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double dlambda = 0 ;
    CustomValues_t<double,ImplicitValues_t>* val1 = (CustomValues_t<double,ImplicitValues_t>*) vim0 ;
    CustomValues_t<double,ExplicitValues_t>* val2 = (CustomValues_t<double,ExplicitValues_t>*) vex0 ;
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      dis[i] = Element_ComputeUnknown(el,u,intfct,p,U_DIS + i) ;
    }
    
    /* Averaging */
    for(int i = 0 ; i < np ; i++) {
      for(int j = 0 ; j < 3 ; j++) w_tot[j] += val1[i].MassFlow_total[j]/np ;

      for(int j = 0 ; j < 9 ; j++) sig[j] += val1[i].Stress[j]/np ;
      
      for(int j = 0 ; j < 9 ; j++) eps_p[j] += val1[i].PlasticStrain[j]/np ;
      
      for(int j = 0 ; j < 9 ; j++) eps[j] += val1[i].Strain[j]/np ;
      
      for(int j = 0 ; j < 3 ; j++) w_air[j] += val1[i].MassFlow_air[j]/np ;
      
      for(int j = 0 ; j < 3 ; j++) w_gas[j] += val1[i].MassFlow_gas[j]/np ;
    
      
      for(int j = 0 ; j < 3 ; j++) {
        j_vap[j] += val2[i].MassFraction_air*val1[i].MassFlow_gas[j] - val1[i].MassFlow_air[j] ;
      }
      
      hardv += val1[i].HardeningVariable/np ;
      dlambda += val1[i].PlasticMultiplier/np ;
    }
    
    tre = eps[0] + eps[4] + eps[8] ;
    e   = (1 + e0) * tre ;
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Liquid_pore_pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_tot    ,"Total_mass_flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,&sl      ,"Saturation_degree",1) ;
    Result_Store(r + i++,&e       ,"Void_ratio_variation",1) ;
    Result_Store(r + i++,eps_p    ,"Plastic_strains",9) ;
    Result_Store(r + i++,&hardv   ,"Hardening_variable",1) ;
    Result_Store(r + i++,&p_g     ,"Gas_pore_pressure",1) ;
    Result_Store(r + i++,w_gas    ,"Gas_mass_flow",3) ;
    Result_Store(r + i++,j_vap    ,"Diffusive_vapor_mass_flow",3) ;
    Result_Store(r + i++,&dlambda ,"Plastic_multiplier",1) ;
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }
  
  return(NbOfOutputs) ;
}



double saturationdegree(double pc,double pc3,Curve_t* curve)
/* Saturation degree: regularization around 1 */
{
  double* x = Curve_GetXRange(curve) ;
  double* y = Curve_GetYValue(curve) ;
  double pc1 = x[0] ;
  double sl ;
  
  if(pc >= pc3 || pc1 >= pc3) {
    sl = Curve_ComputeValue(curve,pc) ;
  } else {
    double sl1 = y[0] ;
    double sl3 = Curve_ComputeValue(curve,pc3) ;
    
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  
  return(sl) ;
}





int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
{
  int ncols = 9 + 2;
  int dec = ncols*ncols;
  double* c0 = c + p*dec;

    /* 1. Derivatives w.r.t. the strain tensor
     * --------------------------------------- */
    if(k == 0) {
      double p_l = val.Pressure_liquid;
      double p_g = val.Pressure_gas;
      double pc  = p_g - p_l ;
      double sl  = SaturationDegree(pc);
      double sg  = 1 - sl;
      
      /* 1.1 Tangent stiffness matrix */
      {
        double* c1 = c0 ;
        double crit = val.YieldCriterion ;
        
        {
          /* Parameters */
          Parameters_t& par = ((Parameters_t*) Element_GetProperty(el))[0] ;
          double poisson = par.PoissonRatio;
          double* vi_n = Element_GetPreviousImplicitTerm(el) ;
          CustomValues_t<double,ImplicitValues_t>& val_n = ((CustomValues_t<double,ImplicitValues_t>*) vi_n)[p] ;
          double* sig_n   = val_n.Stress ;
          double  p_ln    = val_n.Pressure_liquid; // FEM_ComputeUnknown(fem,u_n,intfct,p,U_P_L) ;
          double  p_gn    = val_n.Pressure_gas; // FEM_ComputeUnknown(fem,u_n,intfct,p,U_P_G) ;
          double  pc_n    = p_gn - p_ln ;
          double signet_n = (sig_n[0] + sig_n[4] + sig_n[8])/3. + p_gn ;
          double kappa1   = Kappa(pc_n) ;
          double bulk     = - signet_n*(1 + e0)/kappa1 ;
          //double lame     = bulk - 2*shear/3. ;
          //double poisson  = 0.5 * lame / (lame + mu) ;
          //double young    = 2 * mu * (1 + poisson) ;
          double young    = 3 * bulk * (1 - 2*poisson) ;
          
          Elasticity_SetParameters(elasty,young,poisson) ;
          Elasticity_UpdateStiffnessTensor(elasty) ;
        }

        Elasticity_CopyStiffnessTensor(elasty,c1) ;
      
        /* Criterion */
        if(crit >= 0.) {
          double logp_co  = val.HardeningVariable;// HARDV ;
          double hardv[2] = {logp_co,pc} ;
          double dlambda = val.PlasticMultiplier;
          double sig[9] ;
    
          for(int i = 0 ; i < 9 ; i++) sig[i] = val.Stress[i] ;
    
          /* Net stresses */
          sig[0] += p_g ;
          sig[4] += p_g ;
          sig[8] += p_g ;
            
          /* Continuum tangent stiffness matrix */
          //ComputeTangentStiffnessTensor(sig,hardv) ;
          /* Consistent tangent stiffness matrix */
          //ComputeTangentStiffnessTensor(sig,hardv,&DLAMBDA) ;
          ComputeTangentStiffnessTensor(sig,hardv,&dlambda) ;

          CopyTangentStiffnessTensor(c1) ;
        }
      }
      
      #define B1(i,j) c1[(i)*3+(j)]
      {
      /* 1.2 Coupling matrix for the total mass  */
      {
        double* c1 = c0 + 9*ncols + ncols*((E_MASS < E_MECH) ? E_MASS : E_MECH - E_MASS) ;
        double rho_l = val.MassDensity_liquid; ;
        double rho_g = val.MassDensity_gas; ;
        
        for(int i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*sl + rho_g*sg ;
      }
    
      /* 1.3 Coupling matrix for the air mass */
      {
        double* c1 = c0 + 9*ncols + ncols*((E_AIR < E_MECH) ? E_AIR : E_MECH - E_AIR) ;
        double rho_air = val.MassDensity_air; ;
        
        for(int i = 0 ; i < 3 ; i++) B1(i,i) = rho_air*sg ;
      }
      }
      #undef B1
    }
    
    
    /* 2. Derivatives w.r.t. U_P_L/U_P_G
     * --------------------------------- */
    if(k >= 9) {
      /* 2.1 Mechanical coupling matrix */
      {
        double* c1 = c0 + 9*k ;

        for(int i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
      }
      
      /* 2.2 General storage matrix */
      {
        double* c00 = c0 + 9*ncols + k ;
        
        {
          double* c1 = c00 + ncols*((E_MASS < E_MECH) ? E_MASS : E_MECH - E_MASS);
        
          c1[0] = dval.Mass_total;
        }
        {
          double* c1 = c00 + ncols*((E_AIR < E_MECH) ? E_AIR : E_MECH - E_AIR);
        
          c1[0] = dval.Mass_air;
        }
      }
    }

  return(dec) ;
}

  


void MPM_t::SetIndexes(Element_t* el,int* ind) {
  ind[0] = Values_Index(Strain[0]);
    
  for(int i = 1 ; i < 9 ; i++) {
    ind[i] = ind[0] + i;
  }
    
  ind[9]  = Values_Index(Pressure_liquid);  
  ind[10] = Values_Index(Pressure_gas);
}




void MPM_t::SetIncrements(Element_t* el,double* dui) {  
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
        
  for(int i = 0 ; i < 9 ; i++) {
    dui[i] =  1.e-6 ;
  }
    
  dui[9]  = 1.e-2*ObVal_GetValue(obval + U_P_L) ;
  dui[10] = 1.e-2*ObVal_GetValue(obval + U_P_G) ;
}




int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
{
  int    dec = 4 * 9 ;

  {
    double* c0 = c + p*dec ;
    
    {
      /* Transfer coefficients for the total mass flux */
      {
        double* c1 = c0 ;
        
        c1[0] = dt*val.TransportMass_liquid ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      {
        double* c1 = c0 + 9 ;
        
        c1[0] = dt*val.TransportMass_gas ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      
      /* Transfer coefficients for the air mass flux */
      {
        double* c1 = c0 + 2 * 9 ;
        
        c1[0] = dt*val.TransportAir_liquid ;
        c1[4] = c1[0] ;
        c1[8] = c1[0] ;
      }
      {
        double* c1 = c0 + 3 * 9 ;
        
        c1[0] = dt*val.TransportAir_gas ;
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
  var.DisplacementVectorAndStrainFEM(el,p,U_DIS,val.Displacement) ;
    
  /* Pressure and pressure gradient */
  var.ValueAndGradientFEM(el,p,U_P_L,&val.Pressure_liquid) ;
  
  /* Salt concentration and concentration gradient */
  var.ValueAndGradientFEM(el,p,U_P_G,&val.Pressure_gas) ;
  
  return(&val) ;
}



template <typename T>
Values_t<T>* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_t<T>& val)
/** Compute the secondary variables from the primary ones. */
{
  /* Parameters */
  Parameters_t& par = ((Parameters_t*) Element_GetProperty(el))[0] ;
  double gravity = par.Gravity;
  double rho_s   = par.MassDensity_solidskeleton;
  double kl_int  = par.IntrinsicPermeability_liquid;
  double kg_int  = par.IntrinsicPermeability_gas;
  double mu_l    = par.Viscosity_liquid;
  double mu_g    = par.Viscosity_gas;
  double rho_l0  = par.MassDensity_liquid;
  double* sig0   = par.InitialStress;
  double kappa   = par.SlopeOfElasticLoadingLine;
  //mu      = par.shear_modulus;
  double poisson = par.PoissonRatio;
  double phi0    = par.InitialPorosity;
  double e0      = phi0/(1 - phi0) ;
  double kappa_s = par.SlopeOfElasticWettingDryingLine;
  double d_vap   = par.DiffusionCoefficient_watervapor;
  double d_airliq= par.DiffusionCoefficient_dissolvedair;
  double p_c3    = par.MinimumSuction;
  double henry   = par.HenryConstant;
  /* Inputs 
   * ------*/
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  T* eps   =  val.Strain ;
  double* eps_n =  val_n.Strain ;
  /* Stresses */
  double* sig_n =  val_n.Stress ;
  /* Plastic strains */
  double* eps_pn = val_n.PlasticStrain ;
  /* Pressures */
  T  p_l   = val.Pressure_liquid ;
  double  p_ln  = val_n.Pressure_liquid ;
  T  p_g   = val.Pressure_gas ;
  double  p_gn  = val_n.Pressure_gas ;
  T  pc    = p_g - p_l ;
  double  pc_n  = p_gn - p_ln ;
    

  /* Outputs 
   * ------*/

  /* Backup stresses, plastic strains */
  {
    T* sig   = val.Stress ;
    T* eps_p = val.PlasticStrain ;
    
    {
      T deps[9] ;
      
      /* Incremental deformations */
      for(int i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses at t */
      {
        T trde      = deps[0] + deps[4] + deps[8] ;
        T dlns      = log((pc + p_atm)/(pc_n + p_atm)) ;
        double signet_n  = (sig_n[0] + sig_n[4] + sig_n[8])/3. + p_gn ;
        double kappa1    = Kappa(pc_n) ;
        double bulk      = - signet_n*(1 + e0)/kappa1 ;
        double kappa1_s1  = KappaSuction(-signet_n) ;
        T dsigm     = bulk*trde - signet_n*kappa1_s1/kappa1*dlns ;
        //double lame      = bulk - 2*mu/3. ;
        //double poisson   = 0.5 * lame / (lame + mu) ;
        //double young     = 2 * mu * (1 + poisson) ;
        double young     = 3*bulk*(1 - 2*poisson) ;
        double mu        = 0.5*young/(1 + poisson) ;
        
        #if 1
        if(young < 0) {
          printf("stress:\n") ;
          Math_PrintMatrix(sig_n,3) ;
          printf("signet = %g\n",signet_n) ;
          printf("pg_n = %g\n",p_gn) ;
          printf("pc_n = %g\n",pc_n) ;
          printf("bulk = %g\n",bulk) ;
          printf("young = %g\n",young) ;
          arret("") ;
        }
        #endif
          
        Elasticity_SetParameters(elasty,young,poisson) ;
        Elasticity_UpdateStiffnessTensor(elasty) ;
        
        for(int i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] + 2*mu*deps[i] ;
      
        sig[0] += dsigm - 2*mu*trde/3. ;
        sig[4] += dsigm - 2*mu*trde/3. ;
        sig[8] += dsigm - 2*mu*trde/3. ;
      }
      
      /* Elastic trial net stresses at t  */
      {
        sig[0] += p_g ;
        sig[4] += p_g ;
        sig[8] += p_g ;
      }
    
      /* Plastic strains */
      for(int i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Return mapping */
      {
        double logp_con = val_n.HardeningVariable ; /* log(pre-consolidation pressure) at 0 suction at the previous time step */
        T hardv[2] = {logp_con,pc} ;
        T* crit  = ReturnMapping(sig,eps_p,hardv) ;
        T* dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        val.YieldCriterion  = crit[0] ;
        val.HardeningVariable = hardv[0] ;
        val.PlasticMultiplier = dlambda[0] ;
      }
    
      /* Total stresses */
      sig[0] -= p_g ;
      sig[4] -= p_g ;
      sig[8] -= p_g ;
    }
  }
  
  
  /* Backup mass flows */
  {
    /* Porosity */
    T tre   = eps[0] + eps[4] + eps[8] ;
    T phi   = phi0 + tre ;
    
    /* Partial pressures */
    T p_vap = VaporPressure(pc) ;
    T p_air = p_g - p_vap ;
    
    /* Fluid mass densities */
    double rho_l   = rho_l0 ;
    T rho_vap = M_H2O/RT*p_vap ;
    T rho_air = M_AIR/RT*p_air ;
    T rho_g   = rho_vap + rho_air ;
    T rho_airliq = henry*p_air ;
    
    /* saturation */
    T  sl  = SaturationDegree(pc) ;
    T  sg  = 1 - sl ;
    
    /* Fluid mass flow */
    {
      /* Transfer coefficients at the previous time */
      double kd_mass_l = val_n.TransportMass_liquid ;
      double kd_mass_g = val_n.TransportMass_gas ;
      double k_air_l   = val_n.TransportAir_liquid ;
      double k_air_g   = val_n.TransportAir_gas ;
      //double mc_air    = x_n[I_MC_air] ;
    
      /* Pressure gradients */
      T* gpl = val.GradPressure_liquid ;
      T* gpg = val.GradPressure_gas ;
    
      /* Mass flows */
      T* w_tot = val.MassFlow_total ;
      T* w_air = val.MassFlow_air ;
      T* w_gas = val.MassFlow_gas ;
      
      if(p_air < 0) {
        Message_RuntimeError("Negative air pressure!") ;
      }
    
      for(int i = 0 ; i < 3 ; i++) {
        T w_l   = - kd_mass_l*gpl[i] ;
        T w_g   = - kd_mass_g*gpg[i] ;
        /* w_a      = w_airgas + w_airliq
         * w_airgas = mc_airgas*w_g + j_airgas
         * j_airgas = - rho_g * Fick_vap * Grad(mc_airgas)
         * w_airliq = rho_airliq/rho_l*w_l + j_airliq ~ 0 + j_airliq
         * j_airliq = - Fick_air * Grad(rho_airliq)
         */
        T w_a   = - k_air_g*gpg[i] - k_air_l*gpl[i] ;
        
        w_tot[i] = w_l + w_g ;
        w_air[i] = w_a ;
        w_gas[i] = w_g ;
      }
      w_tot[dim - 1] += (kd_mass_l*rho_l + kd_mass_g*rho_g)*gravity ;
      w_air[dim - 1] += kd_mass_g*rho_air*gravity ;
    }
    
    /* Transfer coefficients at the current time */
    {
      T Darcy_L  = kl_int/mu_l*RelativePermeabilityToLiquid(pc) ;
      T Darcy_G  = kg_int/mu_g*RelativePermeabilityToGas(pc) ;
      T Fick_vap = d_vap*TortuosityToGas(phi,sg) ;
      T Fick_air = d_airliq ;
      T mc_air   = rho_air/rho_g ;
      T mc_vap   = rho_vap/rho_g ;
      T dp_vap   = dVaporPressure(pc) ;
      
      /* w_l = - rho_l*Darcy_L*GradPL */
      val.TransportMass_liquid = rho_l*Darcy_L ;
      /* w_g = - rho_g*Darcy_G*GradPG */
      val.TransportMass_gas = rho_g*Darcy_G ;
      /* w_airgas = mc_airgas*w_g + j_airgas */
      val.TransportAir_gas = rho_air*Darcy_G ;
      /* j_airgas = - rho_g * Fick_vap * Grad(mc_airgas) */
      val.TransportAir_liquid = (mc_air*M_H2O+mc_vap*M_AIR)/RT*dp_vap*Fick_vap ;
      val.TransportAir_gas += - val.TransportAir_liquid + mc_vap*M_AIR/RT*Fick_vap;
      /* w_airliq = rho_airliq/rho_l*w_l + j_airliq = 0 + j_airliq */
      /* j_airliq = - Fick_air * Grad(rho_airliq) */
      val.TransportAir_liquid += dp_vap*Fick_air*henry ;
      val.TransportAir_gas += (1 - dp_vap)*Fick_air*henry ;
    }
    
    /* Mass contents, body force */
    {
      T  m_l   = rho_l*phi*sl ;
      T  m_g   = rho_g*phi*sg ;
      T  m_tot = m_l + m_g ;
      T  m_air = rho_air*phi*sg + rho_airliq*phi*sl ;
      T* f_mass = val.BodyForce ;
    
      val.Mass_total = m_tot ;
      val.Mass_air = m_air ;
      val.MassDensity_liquid = rho_l ;
      val.MassDensity_gas = rho_g ;
      val.MassDensity_air = rho_air ;
      val.Porosity  = phi ;
      val.MassFraction_air  = rho_air/rho_g ;
      
      for(int i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_l)*gravity ;
    }
  }
  
  return(&val) ;
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  /* Parameters */
  Parameters_t& par = ((Parameters_t*) Element_GetProperty(el))[0] ;
  double* sig0   = par.InitialStress;
  DataFile_t* datafile = Element_GetDataFile(el) ;
  
  /* Initial stresses, hardening variables */
  if(DataFile_ContextIsPartialInitialization(datafile)) {
  } else {
    for(int i = 0 ; i < 9 ; i++) val.Stress[i] = sig0[i] ;
    val.HardeningVariable = hardv0 ;
  }
      
  for(int i = 0 ; i < 9 ; i++) val.PlasticStrain[i]  = 0 ;
  val.PlasticMultiplier = 0 ;
    
  return(&val);
}
