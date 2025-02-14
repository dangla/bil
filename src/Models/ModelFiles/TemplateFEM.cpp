#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Plasticity interface */
#include "Plasticity.h"

/* Choose the finite element method */
#include "FEM.h"

#define TITLE "Short title of my model"
#define AUTHORS "Authors"

#include "PredefinedMethods.h"

/*
 * The numbers below are arbitrary and serve only as example
 */

/* Nb of equations of the model */
#define NEQ   (dim+2) /* Let's consider a mechanical eq. and 2 scalar eqs. */


/* Indices of equations */
#define E_MASS1    (0)
#define E_MASS2    (1)
#define E_MECH     (2)

/* Indices of unknowns (generic indices) */
#define U_MECH     E_MECH
#define U_MASS1    E_MASS1
#define U_MASS2    E_MASS2


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */

/* Mechanics */
#define U_DISP     U_MECH
/* Unknown for equation E_MASS1 */
#define U_P1       U_MASS1
/* Unknown for equation E_MASS2 */
#define U_P2       U_MASS2



#include "BaseName_.h"
#include "CustomValues.h"
#include "MaterialPointModel.h"


#define ImplicitValues_t BaseName_(ImplicitValues_t)
#define ExplicitValues_t BaseName_(ExplicitValues_t)
#define ConstantValues_t BaseName_(ConstantValues_t)
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
using Values_t = CustomValues_t<T,ImplicitValues_t,ExplicitValues_t,ConstantValues_t> ;

using Values_d = Values_t<double> ;

#define Values_Index(V)  CustomValues_Index(Values_d,V,double)


#define MPM_t      BaseName(_MPM_t)




struct MPM_t: public MaterialPointModel_t<Values_d> {
  MaterialPointModel_SetInputs_t<Values_d> SetInputs;
  MaterialPointModel_Integrate_t<Values_d> Integrate;
  MaterialPointModel_Initialize_t<Values_d>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_d> SetTangentMatrix;
  MaterialPointModel_SetTransferMatrix_t<Values_d> SetTransferMatrix;
  MaterialPointModel_SetIndexes_t SetIndexes;
  MaterialPointModel_SetIncrements_t SetIncrements;
} ;



using CI_t = ConstitutiveIntegrator_t<Values_d,MPM_t>;





struct ImplicitValues_t {
  double Displacement[3];
  double Strain[9];
  double Pressure1;
  double GradPressure1[3];
  double Pressure2;
  double GradPressure2[3];
  double Mass1;
  double MassFlow1[3];
  double Mass2;
  double MassFlow2[3];
  double Stress[9];
  double BodyForce[3];
} ;



struct ExplicitValues_t {
  double Permeability1;
  double Permeability2;
} ;



struct ConstantValues_t {
  double InitialStress[9];
};




/* Shorthands of some units */
#include "InternationalSystemOfUnits.h"

#define meter    (InternationalSystemOfUnits_OneMeter)
#define m3       (meter*meter*meter)
#define dm       (0.1*meter)
#define cm       (0.01*mmeter
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



/* Material Properties 
 * ------------------- */
#define PCURVE1      (curve1)
#define CURVE1(x)    (Curve_ComputeValue(PCURVE1,x))
#define PCURVE2      (curve2)
#define CURVE2(x)    (Curve_ComputeValue(PCURVE2,x))




#define ComputeTangentStiffnessTensor(...) \
        Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...) \
        Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...) \
        Elasticity_CopyStiffnessTensor(elasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...) \
        Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)




/* Parameters */
/* The parameters below are read in the input data file */

#define Parameters_t    BaseName(_Parameters_t)


struct Parameters_t {
  double Coef1;
  double Coef2;
  double Coef3;
};

static Parameters_t  param;


/* The parameters below are read in the input data file */
static double coef1 ;
static double coef2 ;
static double coef3 ;
static double* sig0 ;
static Curve_t* curve1 ;
static Curve_t* curve2 ;
static Elasticity_t* elasty ;



/* They are stored in the order specified in pm */
static int  pm(const char* s) ;
int pm(const char* s)
{
#define Parameters_Index(V)  CustomValues_Index(Parameters_t,V,double)
  if(strcmp(s,"prop1") == 0)        return (0) ;
  else if(strcmp(s,"prop2") == 0)   return (1) ;
  else if(strcmp(s,"prop3") == 0)   return (2) ;
  else if(!strcmp(s,"sig0"))        return (3) ;
  else if(!strncmp(s,"sig0_",5)) {
     int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
     int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;

     return(3 + 3*i + j) ;
  } else return(-1) ;
#undef Parameters_Index
}



/* They are retrieved automatically by calling the following function */
static void   GetProperties(Element_t*) ;
void GetProperties(Element_t* el)
{
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) 
  coef1  = GetProperty("prop1") ;
  coef2  = GetProperty("prop2") ;
  coef3  = GetProperty("prop3") ;
  sig0   = &GetProperty("sig0") ;
  
  curve1   = Element_FindCurve(el,"y1-axis") ;
  curve2   = Element_FindCurve(el,"y2-axis") ;
#undef GetProperty
}



int SetModelProp(Model_t* model)
/** Set the model properties, return 0.
 *  Warning:
 *  Never call InternationalSystemOfUnits_UseAsLength() or similar
 *  to modify the units because this will also affect other models.
 */
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"mech_1","mech_2","mech_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */ 
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn[i]) ;
  }
  Model_CopyNameOfEquation(model,E_MASS1,"name1") ;
  Model_CopyNameOfEquation(model,E_MASS2,"name2") ;
  
  /** Names of the main (nodal) unknowns */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_DISP + i,name_unk[i]) ;
  }

#if defined (U_P1)
  Model_CopyNameOfUnknown(model,U_P1,"p1") ;
#else
  #error "Ambiguous or undefined unknown"
#endif

#if defined (U_P2)
  Model_CopyNameOfUnknown(model,U_P2,"p2") ;
#else
  #error "Ambiguous or undefined unknown"
#endif
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 3 ;
  
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
  printf("\t- Mechanical Equilibrium           (mech)\n") ;
  printf("\t- First conservation equation      (name1)\n") ;
  printf("\t- Second conservation equation     (name2)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Displacements                    (u_1,u_2,u_3) \n") ;
  printf("\t- First unknown                    (unk1)\n") ;
  printf("\t- Second unknown                   (unk2)\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"prop1 = 0.01   # Property 1\n") ;
  fprintf(ficd,"prop2 = 1.e-3  # Property 2\n") ;
  fprintf(ficd,"prop3 = 1.e6   # Property 3\n") ;
  fprintf(ficd,"Curves = my_file  # File name: x-axis y1-axis y2-axis\n") ;  

  return(0) ;
}



int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;
  int const nvi = CustomValues_NbOfImplicitValues(Values_d);
  int const nve = CustomValues_NbOfExplicitValues(Values_d);
  int const nv0 = CustomValues_NbOfConstantValues(Values_d);

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = nvi*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = nve*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = nv0*NbOfIntPoints ;
  
  /* Skip the rest of code for basic development.
   * For advanced developments find below 
   * some examples of possible operations */
  
  /* Find (and compute) some new interpolation functions */
  #if 0
  {
    int dim = Element_GetDimension(el) ;
    int nn  = Element_GetNbOfNodes(el) ;
    IntFct_t* intfct = Element_GetIntFct(el) ;
    char* typ1 = IntFct_GetType(intfct) ;
    char* typ2 = IntFct_GetType(intfct+1) ;
    
    /* Replace "Type1" and "Type2" by existing types */
    if(strcmp(typ1,"Type1") || strcmp(typ2,"Type2")) {
      int i   = IntFcts_FindIntFct(intfcts,nn,dim,"Type1") ;
      int j   = IntFcts_FindIntFct(intfcts,nn,dim,"Type2") ;
      
      if(j != i+1) {
        i   = IntFcts_AddIntFct(intfcts,nn,dim,"Type1") ;
        j   = IntFcts_AddIntFct(intfcts,nn,dim,"Type2") ;
      }
      /* here j is equal to i + 1 ! */
      Element_GetIntFct(el) = IntFcts_GetIntFct(intfcts) + i ;
    }
    
    NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

    /** Re-define the length of tables */
    Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
    Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
    Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  }
  #endif
  
  
  #if 0 // This part is to be tested!
  {
    /* Remove dof on some nodes 
     * (e.g. at the middle of the edges of a triangle) */
    {
      int nn  = Element_GetNbOfNodes(el) ;
      int j = 3 ; /* dof = 3 to be suppressed */
      int i ;
    
      for(i = nn/2 ; i < nn ; i++) { /* edge nodes */
        Element_GetUnknownPosition(el)[i*NEQ + j]  = -1 ;
        Element_GetEquationPosition(el)[i*NEQ + j] = -1 ;
      }
    }
  }
  #endif
  
  /* Dealing with zero-thickness interface element.
   * (e.g. continuous pressure accross the interface element)
   */
  #if 0
  {
    if(Element_HasZeroThickness(el)) {
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"second_unk") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"second_eqn") ;
    }
  }
  #endif
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* load,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  unsigned int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;

  /* Initialization */
  {
    double zero = 0. ;
    int    i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = zero ;
  }

  {
    
    /* First conservation equation */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(load)) == E_MASS1) {
      double* r1 = FEM_ComputeSurfaceLoadResidu(fem,intfct,load,t,dt) ;
      int    i ;
      
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
    
    /* Second conservation equation */
    } else if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(load)) == E_MASS2) {
    
    /* Other if any */
    } else {
      double* r1 = FEM_ComputeSurfaceLoadResidu(fem,intfct,load,t,dt) ;
      int    i ;
      
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
  double* vi0 = Element_GetImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  CI_t ci(mpm) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

    
  ci.Set(el,t,0,u,vi0,u,vi0) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  
  return(ci.ComputeInitialStateByFEM());
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the (current) explicit terms.
 *  IMPORTANT: if needed use only the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  /* If you need the nodal values, use the previous ones */
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  CI_t ci(mpm) ;
  
  /* If needed ! */
  if(Element_IsSubmanifold(el)) return(0) ;
    
  ci.Set(el,t,0,u,vi_n,u,vi_n) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  
  return(ci.ComputeExplicitTermsByFEM());
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the (current) implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi    = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  CI_t ci(mpm) ;

  if(Element_IsSubmanifold(el)) return(0) ;
    
  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  return(ci.ComputeImplicitTermsByFEM()) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
  double*  vi   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;
  CI_t ci(mpm) ;

  /* Initialization */
  for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  

  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;


  /*
  ** Poromechanical matrix
  */
  {
    double* kp = ci.ComputePoromechanicalMatrixByFEM(E_MECH);

    for(int i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  
  return(0) ;
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
  double* vi1   = Element_GetCurrentImplicitTerm(el) ;
  double* vi1_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;
  CI_t ci(mpm) ;
  
  /* Initialization */
  {    
    for(int i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }

  if(Element_IsSubmanifold(el)) return(0) ;
      
  ci.Set(el,t,dt,u_n,vi1_n,u,vi1) ;
  

  /* Compute here the residu R(n,i) */
  

  /* 1. Mechanics */
  if(Element_EquationIsActive(el,E_MECH)) 
  {
    int istress = Values_Index(Stress[0]);
    int ibforce = Values_Index(BodyForce[0]);
    double* rw = ci.ComputeMechanicalEquilibiumResiduByFEM(istress,ibforce);
    
    for(int i = 0 ; i < nn ; i++) {      
      for(int j = 0 ; j < dim ; j++) R(i,E_MECH + j) -= rw[i*dim + j] ;
    }
  }
  
  
  
  /* 2. Conservation of mass 1 */
  if(Element_EquationIsActive(el,E_MASS1))
  {  
    int imass = Values_Index(Mass1);
    int iflow = Values_Index(MassFlow1[0]);
    double* ra =  ci.ComputeMassConservationResiduByFEM(imass,iflow);
    
    for(int i = 0 ; i < nn ; i++) R(i,E_MASS1) -= ra[i] ;
  }
  
  
  
  /* 3. Conservation of mass 2 */
  if(Element_EquationIsActive(el,E_SALT))
  {  
    int imass = Values_Index(Mass2);
    int iflow = Values_Index(MassFlow2[0]);
    double* ra =  ci.ComputeMassConservationResiduByFEM(imass,iflow);
    
    for(int i = 0 ; i < nn ; i++) R(i,E_MASS2) -= ra[i] ;
  }
  
  return(0) ;
#undef R
}




int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  int NbOfOutputs  = 5 ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
  double* vi = Element_GetImplicitTerm(el) ;
  double* ve = Element_GetExplicitTerm(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  CI_t ci(mpm) ;
  
  /* initialization */
  {
    int    i ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }

  
  ci.Set(el,t,0,u,vi,u,vi) ;
  
  {
    int dim = Element_GetDimensionOfSpace(el) ;
    FEM_t* fem = FEM_GetInstance(el) ;
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Pressure */
    double p1 = Element_ComputeUnknown(el,u,intfct,p,U_P1) ;
    /* Displacement */
    double* dis = Element_ComputeDisplacementVector(el,u,intfct,p,U_DISP) ;
    /* Other quantities */
    double scalar    = 1 ;
    double vector[3] = {1,2,3} ;
    double tensor[9] = {11,12,13,21,22,23,31,32,33} ;
    CustomValues_t<double,ImplicitValues_t>* vali = (CustomValues_t<double,ImplicitValues_t>*) vi ;
    CustomValues_t<double,ExplicitValues_t>* vale = (CustomValues_t<double,ExplicitValues_t>*) ve ;
    
    /* Average some quantities over the element from vi and ve */
    double x_av = 0;
    for(int i = 0 ; i < np ; i++) {
      x_av += vali[i].SomeOtherQuantity/np ;
    }
    
    i = 0  ;
    Result_Store(r + i++,&p1,"Pressure1",1) ;
    
    Result_Store(r + i++,&scalar,"NameOfView_x",1) ; /* scalar */
    Result_Store(r + i++,vector ,"NameOfView_v",3) ; /* vector */
    Result_Store(r + i++,tensor ,"NameOfView_t",9) ; /* tensor */
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}




/** Explanation of the parameters common to all functions.
 *  el    = pointer to element object
 *  t     = current time
 *  dt    = time increment (i.e. the previous time is t-dt)
 *  u     = pointer to pointer to primary nodal unknowns
 *  p     = p^th interpolation Gauss point of the FE
 *  val   = custom values at the current time
 *  val_n = custom values at the previous time (always an input)
 */
Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
/** On output:
 *  val is initialized with the primary nodal unknowns (strain,pressure,temperature,etc...)
 * 
 *  Return a pointer to val
 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
    
  /* Displacements */
  {
    int dim = Element_GetDimensionOfSpace(el) ;
  
    for(int i = 0 ; i < dim ; i++) {
      val.Displacement[i] = Element_ComputeUnknown(el,u,intfct,p,U_MECH + i) ;
    }
    
    for(int i = dim ; i < 3 ; i++) {
      val.Displacement[i] = 0 ;
    }
  }
    
  /* Strain */
  {
    double* eps =  Element_ComputeLinearStrainTensor(el,u,intfct,p,U_MECH) ;
    
    for(int i = 0 ; i < 9 ; i++) {
      val.Strain[i] = eps[i] ;
    }
      
    Element_FreeBufferFrom(el,eps) ;
  }
    
  /* Pressures */
  {
    val.Pressure1 = Element_ComputeUnknown(el,u,intfct,p,U_P1) ;
    val.Pressure2 = Element_ComputeUnknown(el,u,intfct,p,U_P2) ;
  }
  
  return(&val) ;
}

Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
/** On output:
 *  val = initialize val.
 * 
 *  Return a pointer to val.
 */
{
  return(&val);
}


Values_d* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_t<T>& val)
/** On output:
 *  val = update from the integration of the constitutive law from t-dt to t.
 * 
 *  Return a pointer to val.
 **/
{
  /* Strains */
  double* eps   =  val.Strain ;
  double* eps_n =  val_n.Strain ;
  /* Stresses */
  double* sig_n =  val_n.Stress ;
  /* Pressures */
  double  p1   = val.Pressure1 ;
  double  p2   = val.Pressure2 ;  
  
  /* Compute the secondary variables in terms of the primary ones
   * and backup them in x */
  {
  }
  
  return(&val) ;
}

void MPM_t::SetIndexes(Element_t* el,int* ind)
/** On ouput:
 *  ind = a pointer to an array of ncol integers
 *  ncol = nb of the tangent matrix columns
 *  
 *  ind[k] = index of the k^th unknown in the custom values struct
 */
{
}

void MPM_t::SetIncrements(Element_t* el,double* dui)
/** On ouputs:
 *  dui = a pointer to an array of ncol doubles
 *  ncol = nb of the tangent matrix columns
 *  
 *  dui[k] = arbitrary small increment of the k^th unknown used for 
 *           the numerical derivatives (see operator Differentiate)
 */
{
}

int MPM_t::SetTangentMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
/** On input:
 *  k = the k^th column of the tangent matrix to be filled.
 *  dval = the derivatives of val wrt the k^th primary unknown
 * 
 *  On output:
 *  c = pointer to the matrix to be partially filled (only the column k).
 * 
 *  Return the shift (the size of the matrix: ncols*nrows) if succeeds 
 *  or < 0 if fails.
 *  Exemples:
 *    - for an elastic matrix, shift = 81
 *    - for a poroelastic matrix, shift = 100
 */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  int m = 2 ;
  int dec = (9 + m)*(9 + m) ;
  double* c0 = c + p*dec ;

  /* initialization */
  {      
    for(int i = 0 ; i < dec ; i++) c0[i] = 0. ;
  }


  /* Derivatives with respect to strains */
  if(k == 0) {
    /* Mechanics (tangent stiffness matrix) */
    {
      double* c1 = c0 ;
        
      {
        double young   = 1.e9 ;
        double poisson = 0.25 ;
        double* cel = Elasticity_GetStiffnessTensor(elasty) ;
          
        Elasticity_SetParameters(elasty,young,poisson) ;
        Elasticity_ComputeStiffnessTensor(elasty,cel) ;
      }
      
      {
          CopyElasticTensor(c1) ;
      }
    }
        
    
    /* Coupling matrices */
    {
      double* c00 = c0 + 9*(9 + m) ;
    
      /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
        * and zero for the derivatives w.r.t others */
      /* Coupling matrix */
          
      /* First conservation equation */
      {
        double* c1 = c00 ;

        for(int i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass1 ;
      }
          
      /* Second conservation equation */
      {
        double* c1 = c00 + 9 + m ;

        for(int i = 0 ; i < 3 ; i++) B1(i,i) = dval.Mass2 ;
      }
    }
  }
      
      
  /* Derivatives with respect to unk1 */
  if(k >= 9) {
    /* Mechanics: tangent Biot's type coefficient */
    {
      double* c1 = c0 + 9*k ;

      for(int i = 0 ; i < 9 ; i++) c1[i] = dval.Stress[i] ;
    }
    
    /* General storage matrix */
    {
      double* c00 = c0 + 9*(9 + m) + k ;
          
      /* First conservation equation */
      {
        double* c1 = c00 ;
        
        c1[0] = dval.Mass1 ;
      }
          
      /* Second conservation equation */
      {
        double* c1 = c00 + 9 + m ;
        
        c1[0] = dval.Mass2 ;
      }
    }
  }

  return(dec) ;
}

int MPM_t::SetTransferMatrix(Element_t* el,double const& dt,int const& p,Values_d const& val,double* c)
/** On output:
 *  c = pointer to the transfer matrix to be filled
 * 
 *  Return the shift (the size of the matrix: ncols*nrows)
 *  Example: if there are ndif diffusion process, shift = 9*ndif*ndif
 */
{
  int dec = 6*6 ;
  double* c0 = c + p*dec ;
        
  /* initialisation */
  for(int i = 0 ; i < dec ; i++) c0[i] = 0. ;
  
  /* Permeability tensor 1 */
  {
    double* c1 = c0 ;
      
    c1[0] = val.Permeability1 ;
    c1[4] = c1[0] ;
    c1[8] = c1[0] ;
  }
    
  /* Permeability tensor 2 */
  {
    double* c1 = c0 + 27 ;
      
    c1[0] = val.Permeability2 ;
    c1[4] = c1[0] ;
    c1[8] = c1[0] ;
  }
  
  return(dec) ;
}

