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
#define E_Mech     (0)
#define E_Name1    (dim)
#define E_Name2    (1+dim)

/* Indices of unknowns (generic indices) */
#define U_Mech     E_Mech
#define U_Name1    E_Name1
#define U_Name2    E_Name2


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */

/* Mechanics */
#define U_dis     U_Mech
/* Unknown for equation E_Name1 */
#define U_unk1    U_Name1
/* Unknown for equation E_Name2 */
#define U_unk2    U_Name2



/* Nb of terms per point */
#define NVI   (40)    /*  40 implicit terms per point */
#define NVE   (2)     /*  2  explicit terms per point */
#define NV0   (10)    /*  10 constant terms per point */

/* We define some names for implicit terms
 * (vi and vi_n must be used as pointer in the remaining file) */
#define N_1           (vi      )[0] /* Mass or mole content for E_Name1 */
#define N_1n          (vi_n    )[0] /* Idem but at the previous time step */
#define W_1           (vi   + 1)    /* The fluw (3D vector) */

#define N_2           (vi   + 4)[0] /* Idem for E_Name2 */
#define N_2n          (vi_n + 4)[0]
#define W_2           (vi   + 5)

#define SIG           (vi   + 8)    /* Stress tensor (3D tensor) */
#define SIGn          (vi_n + 8)    /* Idem at the previous time step */

#define PHI           (vi   + 17)[0] /* Porosity */

#define F_MASS        (vi   + 26)

#define SCALAR        (vi   + 27)[0]
#define VECTOR        (vi   + 30)
#define TENSOR        (vi   + 39)


/* We define some names for explicit terms 
 * (ve must be used as pointer in the remaining file) */
#define K_1           (ve[0])  /* Transport coef for E_Name1 */
#define K_2           (ve[1])  /* Idem for E_Name2 */


/* We define some names for constant terms (v0 must be used as pointer below) */
#define SIG0          (v0 + 0)  /* Initial stress tensor (3D tensor) */

#define PR_2          (v0 + 9)[0]



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



/* Intern Functions used below */
static int     ComputeTangentCoefficients(Element_t*,double,double,double*) ;
static int     ComputeTransferCoefficients(Element_t*,double,double*) ;

static double* ComputeVariables(Element_t*,void*,void*,void*,const double,const double,const int) ;
static void    ComputeSecondaryVariables(Element_t*,const double,const double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,const double,const double,double*,const double,const int) ;

#define ComputeTangentStiffnessTensor(...) \
        Plasticity_ComputeTangentStiffnessTensor(plasty,__VA_ARGS__)
#define ReturnMapping(...) \
        Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...) \
        Elasticity_CopyStiffnessTensor(elasty,__VA_ARGS__)
#define CopyTangentStiffnessTensor(...) \
        Plasticity_CopyTangentStiffnessTensor(plasty,__VA_ARGS__)


/* We define some indices for the local variables */
static enum {
I_U     =  0,

I_Unk1,

I_Unk2,

I_EPS,
I_EPS8  = I_EPS + 8,

I_GRD_Unk1,
I_GRD_Unk12 = I_GRD_Unk1 + 2,

I_GRD_Unk2,
I_GRD_Unk22 = I_GRD_Unk2 + 2,

I_SIG,
I_SIG8  = I_SIG + 8,

I_N_1,

I_W_1,
I_W_12 = I_W_1 + 2,

I_N_2,

I_W_2,
I_W_22 = I_W_2 + 2,

I_Last
} ;



/* Locally defined intern variables  */
#define NP               IntFct_MaxNbOfIntPoints
#define NbOfVariables    (I_Last)
static double Variable[NP][NbOfVariables] ;
static double Variable_n[NP][NbOfVariables] ;
static double dVariable[NbOfVariables] ;



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
  if(strcmp(s,"prop1") == 0)        return (0) ;
  else if(strcmp(s,"prop2") == 0)   return (1) ;
  else if(strcmp(s,"prop3") == 0)   return (2) ;
  else if(!strcmp(s,"sig0"))        return (3) ;
  else if(!strncmp(s,"sig0_",5)) {
     int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
     int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;

     return(3 + 3*i + j) ;
  } else return(-1) ;
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
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
  Model_CopyNameOfEquation(model,E_Name1,"name1") ;
  Model_CopyNameOfEquation(model,E_Name2,"name2") ;
  
  /** Names of the main (nodal) unknowns */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_dis + i,name_unk[i]) ;
  }

#if defined (U_unk1)
  Model_CopyNameOfUnknown(model,U_unk1,"unk1") ;
#else
  #error "Ambiguous or undefined unknown"
#endif

#if defined (U_unk2)
  Model_CopyNameOfUnknown(model,U_unk2,"unk2") ;
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

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  /* Skip the rest of code for basic development.
   * For advanced developments find below 
   * some examples of possible operations */
  
  /* Compute some new interpolation functions */
  #if 0
  {
    int dim = Element_GetDimension(el) ;
    int nn  = Element_GetNbOfNodes(el) ;
    IntFct_t* intfct = Element_GetIntFct(el) ;
    char* typ1 = IntFct_GetType(intfct) ;
    char* typ2 = IntFct_GetType(intfct+1) ;
    
    /* Replace "Type1" and "Type2" by existing types */
    if(strcmp(typ1,"Type1") || strcmp(typ2,"Type2")) {
      int i   = IntFcts_AddIntFct(intfcts,nn,dim,"Type1") ;
      int j   = IntFcts_AddIntFct(intfcts,nn,dim,"Type2") ; /* not used */
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
  unsigned short int dim = Element_GetDimensionOfSpace(el) ;
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(load)) == E_Name1) {
      double* r1 = FEM_ComputeSurfaceLoadResidu(fem,intfct,load,t,dt) ;
      int    i ;
      
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
    
    /* Second conservation equation */
    } else if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(load)) == E_Name2) {
      double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;
      Symmetry_t sym = Element_GetSymmetry(el) ;
      Node_t** no = Element_GetPointerToNode(el) ;
      Field_t* field = Load_GetField(load) ;
      char*    load_eqn = Load_GetNameOfEquation(load) ;
      char*    load_type = Load_GetType(load) ;
      Function_t* function = Load_GetFunction(load) ;
      int dim_h  = IntFct_GetDimension(intfct) ;
      int    nf  = IntFct_GetNbOfFunctions(intfct) ;
      int    np  = IntFct_GetNbOfPoints(intfct) ;
      int    neq = NEQ ;
      int    ieq = E_Name2 ;
      double* x[Element_MaxNbOfNodes] ;
      int    i ;
  
      for(i = 0 ; i < nn ; i++) {
        x[i] = Node_GetCoordinate(no[i]) ;
      }

      if(field == NULL) return(0) ;

      /* Position index of the equation */
      if(ieq < 0) {
        arret("FEM_ComputeSurfaceLoadResidu (1): unknown equation") ;
      }

      /* Lin-flux is a user-defined type for a linear dependent flux:
       * Flux = A(x,t) * U */
      if(strncmp(load_type,"Lin-flux",4) == 0) {
        double ft = 1. ;
    
        if(function != NULL) {
          ft = Function_ComputeValue(function,t) ;
        }
    
        if(dim == 1 && nn == 1) {
          double v = u[0][U_Name2] ;
          double radius = x[0][0] ;
      
          r[ieq] = ft*Field_ComputeValueAtPoint(field,x[0],dim)*v ;
      
          if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
          else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
        }
    
        if(dim >= 2) {
          double* rb ;
          double f[3*IntFct_MaxNbOfIntPoints] ;
          int p ;
      
          for(p = 0 ; p < np ; p++) {
            double v = FEM_ComputeUnknown(fem,u,intfct,p,U_Name2) ;
            double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
            double y[3] = {0.,0.,0.,} ;
            for(i = 0 ; i < dim ; i++) {
              int j ;
              for(j = 0 ; j < nf ; j++) y[i] += h[j]*x[j][i] ;
            }
            f[p] = ft*Field_ComputeValueAtPoint(field,y,dim)*v ;
          }
      
          rb = FEM_ComputeBodyForceResidu(fem,intfct,f,1) ;
      
          for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
          FEM_FreeBufferFrom(fem,rb) ;
        }
      }
    
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
  double* ve0 = Element_GetExplicitTerm(el) ;
  double* v0  = Element_GetConstantTerm(el) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  

  /* Compute here vi, ve and v0 for each integration points */

  /* Pre-initialization */
  {
    DataFile_t* datafile = Element_GetDataFile(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
    
      /* storage in vi */
      {
        double* vi  = vi0 + p*NVI ;
        int    i ;
      
        /* How to account for partial initialization? */
        if(DataFile_ContextIsPartialInitialization(datafile)) {
          for(i = 0 ; i < 9 ; i++) SIG0[i] = SIG[i] ;
        } else {
          for(i = 0 ; i < 9 ; i++) SIG0[i] = sig0[i] ;
        }

        /* ... */
      }
    
    
      /* storage in ve */
      {
        double* ve  = ve0 + p*NVE ;

        /* ... */
      }
    }
  }


  
    
  /* Loop on integration points */
  {
    double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;
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
      
        /* ... */
      }
    
    
      /* storage in ve */
      {
        double* ve  = ve0 + p*NVE ;

        /* ... */
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
  double* vi0 = Element_GetPreviousImplicitTerm(el) ;
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
      double* x = ComputeVariables(el,u,u,vi0,t,0,p) ;
    
    
      /* storage in ve */
      {
        double* ve  = ve0 + p*NVE ;

        /* ... */
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
      
        /* ... */
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

  /* Compute here the matrix K(i,j) */
  
  /* Example of a poromechanical problem */
  
  /*
  ** Poromechanics matrix
  */
  {
    int m = 2 ; /* Two coupling terms */
    int n = 81 + m*9 + m*(9 + m) ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTangentCoefficients(el,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,m,E_Mech) ;
    /* The entry c should be stored as 
     * (u for displacement, p for unk1, h for unk2)
     * | Kuu  Kup  Kuh  |
     * | Kpu  Kpp  Kph  |
     * | Khu  Khp  Khh  |
     * i.e. the displacements u are in the positions 0 to dim-1 and
     * p is in the position dim, h in dim+1
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
    int m = 2 ;
    int n = m*m*9 ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTransferCoefficients(el,dt,c) ;
  
    /* Conduction Matrix (dW_1/dUnk1) */
    {
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
      int    i ;
      
      for(i = 0 ; i < nn ; i++) {
        int    j ;
      
        for(j = 0 ; j < nn ; j++) {
          K(E_Name1 + i*NEQ,U_unk1 + j*NEQ) += dt*kc[i*nn + j] ;
        }
      }
    }
  
    /* Conduction Matrix (dW_1/dUnk2) */
    {
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c+9,dec) ;
      int    i ;
      
      for(i = 0 ; i < nn ; i++) {
        int    j ;
      
        for(j = 0 ; j < nn ; j++) {
          K(E_Name1 + i*NEQ,U_unk2 + j*NEQ) += dt*kc[i*nn + j] ;
        }
      }
    }
  
    /* Conduction Matrix (dW_2/dUnk1) */
    {
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c+2*9,dec) ;
      int    i ;
      
      for(i = 0 ; i < nn ; i++) {
        int    j ;
      
        for(j = 0 ; j < nn ; j++) {
          K(E_Name2 + i*NEQ,U_unk1 + j*NEQ) += dt*kc[i*nn + j] ;
        }
      }
    }
  
    /* Conduction Matrix (dW_2/dUnk2) */
    {
      double* kc = FEM_ComputeConductionMatrix(fem,intfct,c+3*9,dec) ;
      int    i ;
      
      for(i = 0 ; i < nn ; i++) {
        int    j ;
      
        for(j = 0 ; j < nn ; j++) {
          K(E_Name2 + i*NEQ,U_unk2 + j*NEQ) += dt*kc[i*nn + j] ;
        }
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
  double* vi0 = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  
  
  /* Initialization */
  {
    double zero = 0. ;
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = zero ;
  }

  if(Element_IsSubmanifold(el)) return(0) ;
  

  /* Compute here the residu R(n,i) */
  


  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  {
    double* vi = vi0 ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  {
    double* vi = vi0 ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Mech + dim - 1) -= -rbf[i] ;
    }
    
  }
  
  
  
  /* 2. First conservation equation */
  
  /* 2.1 Accumulation Terms */
  {
    double* vi = vi0 ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI , vi_n += NVI) g1[i] = N_1 - N_1n ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Name1) -= ra[i] ;
    }
  }
  
  /* 2.2 Transport Terms */
  {
    double* vi = vi0 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_1,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Name1) -= -dt*rf[i] ;
  }
  
  
  
  /* 3. Second conservation equation */
  
  /* 3.1 Accumulation Terms */
  {
    double* vi = vi0 ;
    double* vi_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI , vi_n += NVI) g1[i] = N_2 - N_2n ;
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Name2) -= ra[i] ;
    }
  }
  
  /* 3.2 Transport Terms */
  {
    double* vi = vi0 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_2,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Name2) -= -dt*rf[i] ;
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
  double *vi0 = Element_GetImplicitTerm(el) ;
  double *ve0 = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  /* initialization */
  {
    int    i ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }
  
  {
    int dim = Element_GetDimensionOfSpace(el) ;
    FEM_t* fem = FEM_GetInstance(el) ;
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int m = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    
    /* Nodal unknowns */
    double unk1 = FEM_ComputeUnknown(fem,u,intfct,m,U_Name1) ;
    double unk2 = FEM_ComputeUnknown(fem,u,intfct,m,U_Name2) ;
    
    /* Other quantities */
    double scalar    = 1 ;
    double vector[3] = {1,2,3} ;
    double tensor[9] = {11,12,13,21,22,23,31,32,33} ;
    
    int    i ;
    
    /* We can average some quantities over the element */
    for(i = 0 ; i < np ; i++) {
      double * vi = vi0 + i*NVI ;
      double * ve = ve0 + i*NVE ;
      int j ;
      
      scalar += SCALAR/np ;
      
      for(j = 0 ; j < 3 ; j++) vector[j]  += VECTOR[j]/np ;
      for(j = 0 ; j < 9 ; j++) tensor[j]  += TENSOR[j]/np ;
    }
    
    i = 0  ;
    Result_Store(r + i++,&unk1,"Unknown 1",1) ;
    Result_Store(r + i++,&unk2,"Unknown 2",1) ;
    
    Result_Store(r + i++,&scalar,"NameOfView_x",1) ; /* scalar */
    Result_Store(r + i++,vector ,"NameOfView_v",3) ; /* vector */
    Result_Store(r + i++,tensor ,"NameOfView_t",9) ; /* tensor */
    
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
  double*  vi0   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi0_n = Element_GetPreviousImplicitTerm(el) ;
//  double*  ve0  = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double dui[NEQ] ;
  int m = 2 ;
  int dec = 81 + m*9 + m*(9 + m) ;
  
  
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
      double* vi   = vi0   + p*NVI ;
      double* vi_n = vi0_n + p*NVI ;
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vi0_n,t,dt,p) ;
      double* c0 = c + p*dec ;


      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
      
      
      /* The derivative of equations wrt to unknwons */
      
      /* Loop on the unknowns */
    

      /* Derivatives with respect to strains */
      {
        double deps = 1.e-6 ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,deps,I_EPS) ;
    
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
          double* c00 = c0 + 81 + m*9 ;
    
          /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
           * and zero for the derivatives w.r.t others */
          /* Coupling matrix */
          
          /* First conservation equation */
          {
            double* c1 = c00 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_1] ;
          }
          
          /* Second conservation equation */
          {
            double* c1 = c00 + 9 + m ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_2] ;
          }
        }
      }
      
      
      /* Derivatives with respect to unk1 */
      {
        double  dunk1 = dui[U_Name1] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dunk1,I_Unk1) ;
        
        /* Mechanics: tangent Biot's type coefficient */
        {
          double* dsig = dx + I_SIG ;
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsig[i] ;
        }
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + m*9 + 9 ;
          
          /* First conservation equation */
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_N_1] ;
          }
          
          /* Second conservation equation */
          {
            double* c1 = c00 + 9 + m ;
        
            c1[0] = dx[I_N_2] ;
          }
        }
      }
      
      
      /* Derivatives with respect to unk2 */
      {
        double  dunk2 = dui[U_Name2] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dunk2,I_Unk2) ;
        
        /* Mechanics: tangent Biot's type coefficient */
        {
          double* dsig = dx + I_SIG ;
          double* c1 = c0 + 81 + 9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsig[i] ;
        }
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + m*9 + 9 + 1 ;
          
          /* First conservation equation */
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_N_1] ;
          }
          
          /* Second conservation equation */
          {
            double* c1 = c00 + 9 + m ;
        
            c1[0] = dx[I_N_2] ;
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
  double* ve = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  
  for(p = 0 ; p < np ; p++ , ve += NVE) {
    double* c1 = c + p*dec ;
    int    i ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
    
    /* Permeability tensor */
    c1[0] = K_1 ;
    c1[4] = K_1 ;
    c1[8] = K_1 ;
  }
  
  return(dec) ;
}








double* ComputeVariables(Element_t* el,void* vu,void* vu_n,void* vf_n,const double t,const double dt,const int p)
/** This locally defined function compute the intern variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables. */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ; 
  double** u   = (double**) vu ;
  double** u_n = (double**) vu_n ;
  double*  f_n = (double*)  vf_n ;
  double*  x   = Variable ;
  double*  x_n = Variable_n ;
  /* Variables is a locally defined array of array */
  //Model_t*  model  = Element_GetModel(el) ;
  //double*   x      = Model_GetVariable(model,p) ;

  /*
   * The variables are stored in the array x by using some indexes.
   * Two sets of indexes are used:
   * 1. The first set is used for the primary nodal unknowns. The
   *    notation used for these indexes begins with "U_", e.g. "U_Mech", 
   *    "U_Name1", etc.... E.g if displacements, pressure and temperature 
   *    are the nodal unknowns used in this order at each interpolation
   *    point then "U_Mech = 0, U_Name1 = 3, U_Name2 = 4". 
   *    They are used to extract the value of the primary unknowns 
   *    which will be used to compute the secondary variables.
   * 2. The second set is used for the secondary variables. The 
   *    notation used for these indexes begins with "I_". These
   *    secondary variables are stored in the x at these locations.
   */
  
    
  /* Load the primary variables in x */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Unknown 1 */
    x[I_Unk1] = FEM_ComputeUnknown(fem,u,intfct,p,U_Name1) ;
    
    /* Unknown gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Name1) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_Unk1 + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    
    /* Unknown 2 */
    x[I_Unk2] = FEM_ComputeUnknown(fem,u,intfct,p,U_Name2) ;
    
    /* Unknown gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Name2) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_Unk2 + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,const double t,const double dt,double* x_n,double* x)
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Retrieve the primary variables from x */
  /* Strains */
  double* eps   =  x   + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Stresses */
  double* sig_n =  x_n + I_SIG ;
  /* Pressures */
  double  unk1   = x[I_Unk1] ;
  double  unk1_n = x_n[I_Unk1] ;
  
  
  /* Compute the secondary variables in terms of the primary ones
   * and backup them in x */
  {
  }
}


double* ComputeVariableDerivatives(Element_t* el,const double t,const double dt,double* x,const double dxi,const int i)
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
