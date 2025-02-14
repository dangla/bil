#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"
#include "FEM.h"

#ifdef HAVE_AUTODIFF
//#define USE_AUTODIFF
#endif

#define TITLE   "Richards Equation (3D)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ     (1)

/* Equation index */
#define E_LIQ   (0)

/* Unknown index (generic) */
#define U_LIQ   (0)

/* Unknown index */
#define U_P_L   U_LIQ



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
  #ifndef USE_AUTODIFF
  void SetIncrements(Element_t* el,double* dui) { 
    ObVal_t* obval = Element_GetObjectiveValue(el);

    dui[0]  = 1.e-2*ObVal_GetValue(obval + U_LIQ);
  }
  #endif
} ;





/* We define the implicit terms */
template<typename T = double>
struct ImplicitValues_t {
  T Pressure_liquid;
  T GradPressure_liquid[3];
  T Mass_liquid;
  T MassFlow_liquid[3];
} ;



/* We define the explicit terms */
template<typename T = double>
struct ExplicitValues_t {
  T Permeability_liquid;
} ;



/* We define the constant terms */
template<typename T = double>
struct ConstantValues_t {
};



/* We define other terms */
template<typename T = double>
struct OtherValues_t {
  T SaturationDegree_liquid;
};




/* Material properties */
#define SATURATION_CURVE       (Element_GetCurve(el))
#define RELATIVEPERM_CURVE     (Element_GetCurve(el) + 1)

#define SaturationDegree(pc)                Curve_ComputeValue(SATURATION_CURVE,pc)
#define RelativePermeabilityToLiquid(pc)    Curve_ComputeValue(RELATIVEPERM_CURVE,pc)




/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*,double) ;
static double* MacroGradient(Element_t*,double) ;


/* Parameters */
/* The parameters below are read in the input data file */

#define Parameters_t    BaseName(_Parameters_t)


struct Parameters_t {
  double Porosity;
  double IntrinsicPermeability;
  double ReferenceGasPressure;
  double Gravity;
  double LiquidMassDensity;
  double LiquidViscosity;
  double MacroGradient[3];
  double MacroFunctionIndex[3];
};

static Parameters_t  param;
static double macrogradient[3];


double* MacroGradient(Element_t* el,double t)
{
  Functions_t* fcts = Material_GetFunctions(Element_GetMaterial(el)) ;
  Function_t*  fct = Functions_GetFunction(fcts) ;
  int nf = Functions_GetNbOfFunctions(fcts) ;
  double  f[3] = {0,0,0} ;
    
  for(int i = 0 ; i < 3 ; i++) {
    double fctindex = param.MacroFunctionIndex[i] ;
    int idx = floor(fctindex + 0.5) ;
    
    if(0 < idx && idx < nf + 1) {
      Function_t* macrogradfct = fct + idx - 1 ;
      
      f[i] = Function_ComputeValue(macrogradfct,t) ;
    }
        
    macrogradient[i] = param.MacroGradient[i] * f[i] ;
  }
    
  return(macrogradient);
}



int pm(const char *s) {
#define Parameters_Index(V)  CustomValues_Index(Parameters_t,V,double)
  if(!strcmp(s,"Gravity")) {
    return (Parameters_Index(Gravity)) ;
  } else if(!strcmp(s,"Porosity")) {
    return (Parameters_Index(Porosity)) ;
  } else if(!strcmp(s,"LiquidMassDensity")) {
    return (Parameters_Index(LiquidMassDensity)) ;
  } else if(!strcmp(s,"IntrinsicPermeability")) {
    return (Parameters_Index(IntrinsicPermeability)) ;
  } else if(!strcmp(s,"LiquidViscosity")) {
    return (Parameters_Index(LiquidViscosity)) ;
  } else if(!strcmp(s,"ReferenceGasPressure")) {
    return (Parameters_Index(ReferenceGasPressure)) ;
  } else if(!strcmp(s,"MacroGradient")) {
    return(Parameters_Index(MacroGradient[0])) ;
  } else if(!strncmp(s,"MacroGradient_",14)) {
    int i = (strlen(s) > 14) ? s[14] - '1' : 0 ;
    
    if(i < 0 || i >= 3) i = 0 ;
    
    return(Parameters_Index(MacroGradient[i]));
  } else if(!strcmp(s,"MacroFunctionIndex")) {
    return(Parameters_Index(MacroFunctionIndex[0])) ;
  } else if(!strncmp(s,"MacroFunctionIndex_",19)) {
    int i = (strlen(s) > 19) ? s[19] - '1' : 0 ;
    
    if(i < 0 || i >= 3) i = 0 ;
    
    return(Parameters_Index(MacroFunctionIndex[i]));
  } else return(-1) ;
#undef Parameters_Index
}



void GetProperties(Element_t* el,double t)
{
/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
  param = ((Parameters_t*) Element_GetProperty(el))[0] ;
#undef GetProperty
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_LIQ, "liq") ;

  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 12 ;

  /* Par defaut tout a 0 */
  for(int i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;

  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The equation to be solved is:\n") ;
  printf("\t- Mass balance of liquid (liq)\n") ;
  
  printf("\n") ;
  printf("The primary unknown is:\n") ;
  printf("\t- Liquid pressure        (p_l)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"Gravity = -9.81  # La gravite\n") ;
  fprintf(ficd,"Porosity = 0.38       # La porosite\n") ;
  fprintf(ficd,"LiquidMassDensity = 1000     # La masse volumique du fluide\n") ;
  fprintf(ficd,"IntrinsicPermeability = 8.9e-12  # La permeabilite intrinseque\n") ;
  fprintf(ficd,"LiquidViscosity = 0.001     # La viscosite du fluide\n") ;
  fprintf(ficd,"ReferenceGasPressure = 1.e5       # La pression du gaz\n") ;
  fprintf(ficd,"Curves = billes  # Le nom du fichier p_c S_l k_rl\n") ;
  
  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  IntFct_t *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;
    
  MaterialPointModel_DefineNbOfInternalValues(MPM_t,el,NbOfIntPoints);
  
  /* Continuity of pressure across zero-thickness element */
  {
    if(Element_HasZeroThickness(el)) {
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"p_l") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"liq") ;
    }
  }
  
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  IntFct_t *fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  
  {
    double *r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    int i ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}




int ComputeInitialState(Element_t* el,double t)
{
  int i = MaterialPointModel_ComputeInitialStateByFEM(MPM_t,el,t);
  
  return(i);
}


int  ComputeExplicitTerms(Element_t* el,double t)
{
  int i = MaterialPointModel_ComputeExplicitTermsByFEM(MPM_t,el,t);
  
  return(i);
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  int i = MaterialPointModel_ComputeImplicitTermsByFEM(MPM_t,el,t,dt);
  
  return(i);
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
{
  int i = MaterialPointModel_ComputeMassConservationMatrixByFEM(MPM_t,el,t,dt,k);
  
  #if 0
  {
    int ndof = Element_GetNbOfDOF(el) ;
    printf("\n");
    printf("Element %d\n",Element_GetElementIndex(el));
    Math_PrintMatrix(k,ndof);
  }
  #endif
  
  return(i);
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
{
  /* Initialization */
  {
    int ndof = Element_GetNbOfDOF(el) ;
    
    for(int i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }
  
  /* Conservation of liquid mass */
  MaterialPointModel_ComputeMassConservationResiduByFEM(MPM_t,el,t,dt,r,E_LIQ,Mass_liquid,MassFlow_liquid);
  
  return(0);
}





int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 3 ;
  double* vim = Element_GetImplicitTerm(el) ;
  double* vex = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;

  //if(Element_IsSubmanifold(el)) return(0) ;

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
    double pg  = param.ReferenceGasPressure;
    /* pressures */
    double pl  = Element_ComputeUnknown(el,u,intfct,p,U_P_L) ;
    double pc  = pg - pl ;
    
    /* saturation */
    double sl  = SaturationDegree(pc) ;
    
    double w_l[3] = {0,0,0} ;
    int    i ;
    
    /* Averaging */
    for(int i = 0 ; i < np ; i++) {
      CustomValues_t<double,ImplicitValues_t>* val1 = (CustomValues_t<double,ImplicitValues_t>*) vim ;
      
      for(int j = 0 ; j < 3 ; j++) w_l[j]  += val1[i].MassFlow_liquid[j]/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&pl,"pressure",1) ;
    Result_Store(r + i++,w_l,"flow",3) ;
    Result_Store(r + i++,&sl,"saturation",1) ;
  }
  
  return(NbOfOutputs) ;
}






int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
{
  int   ncols = 1;
  int   dec = ncols*ncols ;
  double* c0 = c + p*dec ;
      
  /* Derivatives w.r.t U_LIQ */
  if(k == 0) {
    /* General storage matrix */
    double* c1 = c0 ;

    c1[0] = dval.Mass_liquid ;
  }

  return(dec) ;
}




int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
{
  int    dec = 9 ;

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
    }
  }

  return(dec) ;
}



void MPM_t::SetIndexes(Element_t* el,int* ind) {
  ind[0]  = Values_Index(Pressure_liquid); 
}




Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
{
  LocalVariables_t<Values_d> var(u,NULL);
    
  /* Pressure and pressure gradient */
  var.ValueAndGradientFEM(el,p,U_LIQ,&val.Pressure_liquid) ;
  
  if(Geometry_IsPeriodic(Element_GetGeometry(el))) {
    double* grad = MacroGradient(el,t);
    IntFct_t* intfct = Element_GetIntFct(el) ;
    double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* coor = Element_ComputeCoordinateVector(el,h) ;
      
    for(int i = 0 ; i < 3 ; i++) {
      val.GradPressure_liquid[i] += grad[i] ;
      val.Pressure_liquid        += grad[i]*coor[i] ;
    }
  }
  
  return(&val) ;
}



template <typename T>
Values_t<T>* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_t<T>& val)
{
  int dim = Element_GetDimensionOfSpace(el);
  /* Parameters */
  double rho_l   = param.LiquidMassDensity;
  double phi     = param.Porosity;
  double gravity = param.Gravity;
  double pg      = param.ReferenceGasPressure;
  double mu_l    = param.LiquidViscosity;
  double k_int   = param.IntrinsicPermeability;
  /* Pressures */
  T pl  = val.Pressure_liquid;
  T pc  = pg - pl;
    
  /* saturation */
  T sl  = SaturationDegree(pc);
    
  /* liquid mass content */
  T m_l = rho_l*phi*sl;
  
  
  /* Backup */

  val.SaturationDegree_liquid = sl;
  val.Mass_liquid = m_l;
  
  /* Fluxes */
  {
    /* Transfer coefficients at the previous time */
    double k_l = val_n.Permeability_liquid;
    
    /* Fluxes */
    T* w_l = val.MassFlow_liquid;
    
    /* Gradients */
    T* gpl = val.GradPressure_liquid;
    
    for(int i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i];
    w_l[dim-1] += k_l*rho_l*gravity;
  }
  
  /* Transfer coefficients at the current time */
  {
    T kr_l    = RelativePermeabilityToLiquid(pc);
    T k_l     = rho_l*k_int/mu_l*kr_l;
      
    val.Permeability_liquid  = k_l;
  }
  
  return(&val) ;
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  /* Parameters */
  double rho_l   = param.LiquidMassDensity;
  double pg      = param.ReferenceGasPressure;
  double mu_l    = param.LiquidViscosity;
  double k_int   = param.IntrinsicPermeability;
  /* Pressures */
  double pl   = val.Pressure_liquid;
  double pc   = pg - pl;
  double kr_l = RelativePermeabilityToLiquid(pc);
  double k_l  = rho_l*k_int/mu_l*kr_l;

  val.Permeability_liquid  = k_l;
  
  return(&val);
}
