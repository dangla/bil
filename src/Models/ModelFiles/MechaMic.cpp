#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM2.h"
#include "FEM.h"
#include "DataSet.h"
#include "Solvers.h"
#include "Modules.h"

#define TITLE "Mechanics from a microstructure (2017)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (dim)

/* Equation index */
#define E_MECH   (0)

/* Indices of unknowns (generic indices) */
#define U_MECH     E_MECH

/* Unknown index */
#define U_DISP     U_MECH



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
  MaterialPointModel_Integrate_t<Values_t> Integrate;
  MaterialPointModel_Initialize_t<Values_t>  Initialize;
  MaterialPointModel_SetTangentMatrix_t<Values_t> SetTangentMatrix;
} ;




/* We define some names for implicit terms */
template<typename T = double>
struct ImplicitValues_t {
  T Displacement[3];
  T Strain[9];
  T Stress[9];
  T BodyForce[3];
  Solutions_t* Solutions;
} ;


template<typename T = double>
struct OtherValues_t {
  int InterpolationPointIndex;
};



/* We define some names for explicit terms */
template<typename T = double>
struct ExplicitValues_t {
} ;



/* We define some names for constant terms (v0 must be used as pointer below) */
template<typename T = double>
struct ConstantValues_t {
  T InitialStress[9];
};



using CI_t = ConstitutiveIntegrator_t<MPM_t>;
//using CI_t = ConstitutiveIntegrator_t<MaterialPointModel_t<Values_t>>;

static MPM_t mpm1;
static MPM_t* mpm = &mpm1;
//static MaterialPointModel_t<Values_t>* mpm = &mpm1;

static CI_t ci(mpm) ;



/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void    GetProperties(Element_t*,double) ;
static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;

/* Parameters */
//static double  gravity ;
//static double  rho_s ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))



int pm(const char* s)
{
         if(!strcmp(s,"gravity")) {
    return(0) ;
  } else if(!strcmp(s,"rho_s")) {
    return(1) ;
  } else if(!strcmp(s,"macro-gradient")) {
    return(2) ;
  } else if(!strncmp(s,"macro-gradient_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(2 + 3*i + j) ;
  } else if(!strcmp(s,"macro-fctindex")) {
    return(11) ;
  } else if(!strncmp(s,"macro-fctindex_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(11 + 3*i + j) ;
  } else if(!strcmp(s,"charlen")) {
    return(20) ;
  } else return(-1) ;
}



double* MacroGradient(Element_t* el,double t)
{
  double* gradient = macrogradient ;
  double  f[9] = {0,0,0,0,0,0,0,0,0} ;
  
  {
    Functions_t* fcts = Material_GetFunctions(Element_GetMaterial(el)) ;
    Function_t*  fct = Functions_GetFunction(fcts) ;
    int nf = Functions_GetNbOfFunctions(fcts) ;
    double* fctindex = &GetProperty("macro-fctindex") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      int idx = floor(fctindex[i] + 0.5) ;
      
      if(0 < idx && idx < nf + 1) {
        Function_t* macrogradfct = fct + idx - 1 ;
      
        f[i] = Function_ComputeValue(macrogradfct,t) ;
      }
    }
  }
  
  {
    double* g = &GetProperty("macro-gradient") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      gradient[i] = g[i] * f[i] ;
    }
  }
  
  return(gradient) ;
}



double* MacroStrain(Element_t* el,double t)
{
  double* strain = macrostrain ;
  double* grd = MacroGradient(el,t) ;
  int i ;
    
  for(i = 0 ; i < 3 ; i++) {
    int j ;
      
    for(j = 0 ; j < 3 ; j++) {
      strain[3*i + j] = 0.5*(grd[3*i + j] + grd[3*j + i]) ;
    }
  }
  
  return(strain) ;
}



void GetProperties(Element_t* el,double t)
{
  //gravity = GetProperty("gravity") ;
  //rho_s   = GetProperty("rho_s") ;
    
  {
    DataSet_t* dataset = Element_GetDataSet(el) ;
    DataFile_t* df = DataSet_GetDataFile(dataset) ;
    
    DataFile_ContextSetToPartialInitialization(df) ;
  }
}



int SetModelProp(Model_t* model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  char i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  for(i = 0 ; i < dim ; i++) {
    char name_eqn[7] ;
    sprintf(name_eqn,"meca_%u",i + 1) ;
    Model_CopyNameOfEquation(model,E_MECH + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%u",i + 1) ;
    Model_CopyNameOfUnknown(model,U_DISP + i,name_unk) ;
  }
  
  Model_GetComputePropertyIndex(model) = &pm ;
  Model_GetComputeMaterialProperties(model) = &GetProperties;
    
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int NbOfProp = 21 ;

  /* Par defaut tout a 0 */
  {
    int i ;
    
    for(i = 0 ; i < NbOfProp ; i++) {
      Material_GetProperty(mat)[i] = 0. ;
    }
  }
  
  /* Pre-initialization */
  {
    //Material_GetProperty(mat)[pm("charlen")] = 1 ;
  }
  

  Material_ScanProperties(mat,datafile,pm) ;


  /* The dataset and solver for the microstructure */
  {
    char* method = Material_GetMethod(mat) ;
    
    if(!strncmp(method,"Microstructure",14)) {
      char* p = strstr(method," ") ;
      char* cellname = p + strspn(p," ") ;
      Options_t* options = Options_Create(NULL) ;
      DataSet_t* dataset = DataSet_Create(cellname,options) ;
      Mesh_t* mesh = DataSet_GetMesh(dataset) ;
      Solvers_t* solvers = Solvers_Create(mesh,options,6) ;
      
      /* Store options in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,options,Options_t,"Options") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }
      
      /* Store dataset in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,dataset,DataSet_t,"DataSet") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }

      
      /* Store solvers in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,solvers,Solvers_t,"Solvers") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }

      /* Initialize the dataset */
      {
        FEM2_t* fem2 = FEM2_GetInstance(dataset,NULL,NULL,NULL) ;
        
        FEM2_InitializeMicrostructureDataSet(fem2) ;
      }
    }
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

  printf("\n") ;
  printf("The system consists in dim equations\n") ;
  printf("\t The equilibrium equations (meca_1,meca_2,meca_3)\n") ;

  printf("\n") ;
  printf("The primary unknowns are\n") ;
  printf("\t The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n") ;
  printf("Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  
  return(0) ;
}



/* The current and the previous solutions */
#define Element_GetCurrentLocalSolutions(EL) \
        ((Solutions_t*) Element_FindCurrentImplicitData(EL,Solutions_t,"Solutions"))

#define Element_GetPreviousLocalSolutions(EL) \
        ((Solutions_t*) Element_FindPreviousImplicitData(EL,Solutions_t,"Solutions"))


/* The dataset */
#define Element_GetDataSet(EL) \
        ((DataSet_t*) Element_FindMaterialData(EL,DataSet_t,"DataSet"))


/* The solvers */
#define Element_GetSolvers(EL) \
        ((Solvers_t*) Element_FindMaterialData(EL,Solvers_t,"Solvers"))




int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
  MaterialPointModel_DefineNbOfInternalValues(MPM_t,el,NbOfIntPoints);


  /** The solutions from the microstructures */
  {
    DataSet_t* dataset = Element_GetDataSet(el) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    const int nsol_micro = 2 ;
    
    /* Store "sols" for the whole history */
    Element_AllocateMicrostructureSolutions(el,mesh,nsol_micro) ;
  }
  
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

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
  
    {
      int i ;
      
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}



int ComputeInitialState(Element_t* el,double t)
{
  #if 1
  double* vim0  = Element_GetImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el,t) ;
      
  ci.Set(el,t,0,u,vim0,u,vim0) ;
  
  {
    int i = ci.ComputeInitialStateByFEM();
    
    return(i);
  }
  #else
  int i = MaterialPointModel_ComputeInitialStateByFEM(MPM_t,el,t);
  
  return(i);
  #endif
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  #if 1
  double* vim0   = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el,t) ;
    
  ci.Set(el,t,dt,u_n,vim_n,u,vim0) ;
  
  {
    int i = ci.ComputeImplicitTermsByFEM();
    
    return(i);
  }
  #else
  int i = MaterialPointModel_ComputeImplicitTermsByFEM(MPM_t,el,t,dt);
  
  return(i);
  #endif
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
  #if 1
  #define K(i,j)    (k[(i)*ndof + (j)])
  double*  vi   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;


  /* Initialization */
  {
    for(int i = 0 ; i < ndof*ndof ; i++) k[i] = 0 ;
  }


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el,t) ;
      
  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;

  {
    double* kp = ci.ComputeTangentStiffnessMatrixByFEM() ;

    for(int i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  
  return(0) ;
  #undef K
  #else
  int i = MaterialPointModel_ComputeTangentStifnessMatrixByFEM(MPM_t,el,t,dt,k);
  
  return(i);
  #endif
}




int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
  #if 1
  #define R(n,i)    (r[(n)*NEQ+(i)])
  double*  vi   = Element_GetCurrentImplicitTerm(el) ;
  double*  vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;

  /* Initialization */
  for(int i = 0 ; i < ndof ; i++) r[i] = 0 ;

  if(Element_IsSubmanifold(el)) return(0) ;

  GetProperties(el,t) ;
      
  ci.Set(el,t,dt,u_n,vi_n,u,vi) ;

  {
    int istress = Values_Index(Stress[0]);
    int ibforce = Values_Index(BodyForce[0]);
    double* rw = ci.ComputeMechanicalEquilibiumResiduByFEM(istress,ibforce);
    
    for(int i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_MECH + j) -= rw[i*dim + j] ;
    }
  }
  #undef R
  
  return(0) ;
  #else
  /* Initialization */
  {
    int ndof = Element_GetNbOfDOF(el) ;
    
    for(int i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }
  
  MaterialPointModel_ComputeMechanicalEquilibriumResiduByFEM(MPM_t,el,t,dt,r,E_MECH,Stress,BodyForce);

  return(0);
  #endif
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 4 ;
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;

  //if(Element_IsSubmanifold(el)) return(0) ;

  /* Initialization */
  {
    for(int i = 0 ; i < NbOfOutputs ; i++) {
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
    /* Displacement */
    double* pdis = Element_ComputeDisplacementVector(el,u,intfct,p,U_DISP) ;
    double dis[3] = {0,0,0} ;
    /* strains */
    double *eps =  Element_ComputeLinearStrainTensor(el,u,intfct,p,U_DISP) ;
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    
    /* Total displacements and strains */
    {
      for(int i = 0 ; i < dim ; i++) {
        dis[i] = pdis[i] ;
      }

      if(ItIsPeriodic) {
        for(int i = 0 ; i < 3 ; i++) {
          int j ;
        
          for(j = 0 ; j < 3 ; j++) {
            dis[i] += MacroGradient(el,t)[3*i + j] * s[j] ;
          }
        }
        
        for(int i = 0 ; i < 9 ; i++) {
          eps[i] += MacroStrain(el,t)[i] ;
        }
      }
    }
    
    /* Averaging */
    for(p = 0 ; p < np ; p++) {
      CustomValues_t<double,ImplicitValues_t>* val1 = (CustomValues_t<double,ImplicitValues_t>*) vim0 ;

      for(int j = 0 ; j < 9 ; j++) sig[j] += val1[p].Stress[j]/np ;
    }
      
    {
      int i = 0 ;
      
      Result_Store(r + i++,dis   ,"Displacements",3) ;
      Result_Store(r + i++,sig   ,"Stresses",9) ;
      Result_Store(r + i++,pdis  ,"Perturbated-displacements",3) ;
      Result_Store(r + i++,eps   ,"Strains",9) ;
      
      if(i != NbOfOutputs) arret("ComputeOutputs") ;
    }
  }
  
  return(NbOfOutputs) ;
}



int MPM_t::SetTangentMatrix(Element_t* el,double const& t,double const& dt,int const& p,Values_d const& val,Values_d const& dval,int const& k,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
  DataSet_t* dataset = Element_GetDataSet(el) ;
  DataFile_t* df = DataSet_GetDataFile(dataset) ;
  Solvers_t* solvers = Element_GetSolvers(el) ;
  Solutions_t* sols   = Element_GetCurrentLocalSolutions(el) ;
  int    dec = 81 ;

  if(k == 0) {    
    {
      double* c0 = c + p*dec ;
      
      {
        Session_Open() ;
        Message_SetVerbosity(0) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure matrix (tangent matrix): %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        {
          FEM2_t* fem2 = FEM2_GetInstance(dataset,solvers,NULL,sols+p) ;
          int i = FEM2_HomogenizeTangentStiffnessTensor(fem2,t,dt,c0) ;
        
          if(i < 0) {
            Message_FatalError("ComputeTangentCoefficients: something went wrong") ;
          }
        }
        
        Session_Close() ;
      }
      
#if 0
      {
        int i ;
      
        printf("\n") ;
        printf("4th rank stiffness tensor:\n") ;
      
        for(i = 0 ; i < 9 ; i++) {
          int j = i - (i/3)*3 ;
        
          printf("C%d%d--:",i/3 + 1,j + 1) ;
        
          for (j = 0 ; j < 9 ; j++) {
            printf(" % e",c0[i*9 + j]) ;
          }
        
          printf("\n") ;
        }
      }
#endif
    }
  }
  
  return(dec) ;
}



Values_d* MPM_t::SetInputs(Element_t* el,const double& t,const int& p,double const* const* u,Values_d& val)
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
  
  if(ItIsPeriodic) {
    for(int i = 0 ; i < 9 ; i++) {
      val.Strain[i]   += MacroStrain(el,t)[i] ;
    }
  }
  
  val.InterpolationPointIndex = p ;
  
  return(&val) ;
}



Values_d* MPM_t::Integrate(Element_t* el,const double& t,const double& dt,Values_d const& val_n,Values_d& val)
{
  /* Strains */
  double* eps   =  val.Strain ;
  double* eps_n =  val_n.Strain ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = val.Stress ;
    double* sig_n = val_n.Stress ;
    
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
      
      /* Compute here the microstructure */
      {
        DataSet_t* dataset = Element_GetDataSet(el) ;
        Solvers_t* solvers = Element_GetSolvers(el) ;
        DataFile_t* df = DataSet_GetDataFile(dataset) ;

        
        Session_Open() ;
        Message_SetVerbosity(99) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure: %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        {
          int p = val.InterpolationPointIndex;
          Solutions_t* sols_n = val_n.Solutions;
          Solutions_t* sols = Element_GetCurrentLocalSolutions(el) ;
          FEM2_t* fem2 = FEM2_GetInstance(dataset,solvers,sols_n+p,sols+p) ;
          int j = FEM2_ComputeHomogenizedStressTensor(fem2,t,dt,deps,sig) ;
          
          val.Solutions = sols;
        
          if(j < 0) {
            Message_FatalError("ComputeSecondaryVariables: something went wrong") ;
          }
        }
        
        Session_Close() ;
      }

    }
  }
  
  
  /* Backup body force */
  {
    {
      double gravity = GetProperty("gravity") ;
      double rho_s   = GetProperty("rho_s") ;
      int dim = Element_GetDimensionOfSpace(el) ;
      double* f_mass = val.BodyForce ;
      
      for(int i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s)*gravity ;
    }
  }
  
  return(&val);
}




Values_d* MPM_t::Initialize(Element_t* el,double const& t,Values_d& val)
{
  DataSet_t* dataset = Element_GetDataSet(el) ;
  DataFile_t* df = DataSet_GetDataFile(dataset) ;
  DataFile_t* datafile = Element_GetDataFile(el) ;
    
  if(DataFile_ContextIsPartialInitialization(datafile)) {
    DataFile_ContextSetToPartialInitialization(df) ;
  } else {
    DataFile_ContextSetToFullInitialization(df) ;
  }
    
  {
    Solutions_t* sols = Element_GetCurrentLocalSolutions(el) ;
    
    val.Solutions = sols;
  }
  
  return(&val) ;
}
