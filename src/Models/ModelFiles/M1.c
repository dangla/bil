#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"
#include "FVM.h"

#define TITLE "Richards Equation (1D)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ (1)
#define NVE (1)
#define NVI (3)

#define E_liq (0)
#define I_p_l (0)


/* Value of the nodal unknown (u and el must be used as pointers below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)


/* We define some names for nodal unknowns */
#define P_l(n)          (UNKNOWN(n,I_p_l))


/* We define some names for implicit terms (f must be used as pointer below) */
#define M_l(n)   (f[(n)])
#define W_l      (f[(2)])
#define M_ln(n)  (f_n[(n)])


/* We define some names for explicit terms (va must be used as pointer below) */
#define K_l      (va[(0)])


/* Curves */
#define SATURATION(x)    (Curve_ComputeValue(Element_GetCurve(el),x))
#define DSATURATION(x)   (Curve_ComputeDerivative(Element_GetCurve(el),x))
#define RELATIVEPERM(x)  (Curve_ComputeValue(Element_GetCurve(el) + 1,x))


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])



/* Fonctions */
static int    pm(const char *s) ;
//static int    c1(Element_t*,double*) ;
//static int    k1(Element_t*,double*) ;
static int    TangentCoefficients(Element_t*,double,double*) ;


/* Parametres */
static double gravite,phi,rho_l,k_int,mu_l,p_g,schema ;
#if SharedMS_APIis(OpenMP)
  #pragma omp threadprivate(gravite,phi,rho_l,k_int,mu_l,p_g,schema)
#endif


int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"p_g") == 0) return (5) ;
  else if(strcmp(s,"schema") == 0) return (6) ;
  else return(-1) ;
}


int SetModelProp(Model_t *model)
{
  Model_CopyShortTitle(model,TITLE) ;
  Model_GetNbOfEquations(model) = NEQ ;
  Model_CopyNameOfEquation(model,E_liq,"liq") ;
  Model_CopyNameOfUnknown(model,I_p_l,"p_l") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int dim = Material_GetDimension(mat) ;
  int  n_donnees = 7 ;

  if(dim > 1) arret("ReadMatProp : dimension > 1 non prevue") ;
  
  Material_GetProperty(mat)[pm("schema")] = -1 ;
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(n_donnees) ;
}


int PrintModelProp(Model_t *model,FILE *ficd)
/* Print Model Properties */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  fprintf(stdout,"\n\n\
L\'inconnue est la pression de liquide : p_l.\n\
Exemple de donnees :\n\n") ;

  fprintf(ficd,"gravite = -9.81  # La gravite\n") ;
  fprintf(ficd,"phi = 0.38       # La porosite\n") ;
  fprintf(ficd,"rho_l = 1000     # La masse volumique du fluide\n") ;
  fprintf(ficd,"k_int = 8.9e-12  # La permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001     # La viscosite du fluide\n") ;
  fprintf(ficd,"p_g = 1.e5       # La pression du gaz\n") ;
  fprintf(ficd,"courbes = billes # Le nom du fichier p_c S_l k_rl\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  /*
    Donnees
  */
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  p_g     = GetProperty("p_g") ;

  /* MASSE LIQUIDE */
  for(i=0;i<2;i++) {
    double pc  = p_g - P_l(i) ;
    double sl  = SATURATION(pc) ;
    
    M_l(i) = rho_l*phi*sl ;
  }
  
  /* COEFFICIENTS DE TRANSFERT */
  {
    double pl   = (P_l(0) + P_l(1))*0.5 ;
    double pc   = p_g - pl ;
    double k_rl = RELATIVEPERM(pc) ;
    
    K_l = rho_l*k_int/mu_l*k_rl ;
  }
  
  /* FLUX */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    
    W_l = - K_l*grd_p_l + K_l*rho_l*gravite ;
  }
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  /*
    Donnees
  */
  gravite = GetProperty("gravite") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  p_g     = GetProperty("p_g") ;
  schema  = GetProperty("schema") ;
  
  /*
    COEFFICIENTS DE TRANSFERT
  */
  {
    double pl = (P_l(0) + P_l(1))*0.5 ;
    double pc ;
    double k_rl ;
    
    if(schema == 1) {
      double x1 = Element_GetNodeCoordinate(el,1)[0] ;
      double x0 = Element_GetNodeCoordinate(el,0)[0] ;
      double dx = x1 - x0 ;
      double grdh = (P_l(1) - P_l(0))/dx - rho_l*gravite ;
      
      pl = (grdh > 0) ? P_l(1) : P_l(0) ;
    }
    
    pc   = p_g - pl ;
    k_rl = RELATIVEPERM(pc) ;
    K_l  = rho_l*k_int/mu_l*k_rl ;
  }
  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  /*
    Donnees
  */
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  p_g     = GetProperty("p_g") ;

  /* masse de fluide */
  for(i=0;i<2;i++) {
    double pc = p_g - P_l(i) ;
    double sl  = SATURATION(pc) ;
    
    M_l(i) = rho_l*phi*sl ;
  }

  /* flux */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    
    W_l = - K_l*grd_p_l + K_l*rho_l*gravite ;
  }
      
  #if 0
  {
    int id = SharedMS_CurrentThreadId ;
        
    if(id == 0) {
      int nthreads = omp_get_num_threads() ;
          
      printf("Actual number of threads: %d\n",nthreads) ;
    }
        
    printf("Current thread id: %d\n",id) ;
  }
  #endif
  
  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  int    nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  double c[Element_MaxNbOfNodes*Element_MaxNbOfNodes*NEQ*NEQ] ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Data
  */
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  p_g     = GetProperty("p_g") ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,1) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  
  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume ;
  double surf ;
  int    i ;
  double zero = 0. ;

  /* initialisation */
  for(i = 0 ; i < nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  volume = FVM_ComputeCellVolumes(fvm) ;
  {
    double *surface = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = surface[1] ;
  }

  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  r[0] -= volume[0]*(M_l(0) - M_ln(0)) + dt*surf*W_l ;
  r[1] -= volume[1]*(M_l(1) - M_ln(1)) - dt*surf*W_l ;
  return(0) ;
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *va = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;
  int    nso = 3 ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Donnees
  */
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  p_g     = GetProperty("p_g") ;

  /* initialisation */
  for(i = 0 ; i < nso ; i++) Result_SetValuesToZero(r+i) ;

  /* quantites exploitees */
  {
    int    j = FVM_FindLocalCellIndex(fvm,s) ;
    /* pression */
    double pl =  P_l(j) ;
    /* saturation */
    double pc = p_g - pl ;
    double sl  = SATURATION(pc) ;
    /* flux */
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    double wlx = - K_l*grd_p_l + K_l*rho_l*gravite ;
    double wl[3]  = {0,0,0} ;
    
    Result_Store(r + 0,&pl,"pression-liquide",1) ;
    wl[0] = wlx ;
    Result_Store(r + 1,wl,"flux-liquide",3) ; 
    Result_Store(r + 2,&sl,"saturation",1) ;
  }

  return(nso) ;
}



int c1(Element_t *el,double *c)
/*
**  Mass matrix (c) and shift (dec)
*/
{
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nn = Element_GetNbOfNodes(el) ;
  int    dec = 1 ;
  int    i ;
  /*
    Data
  */
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  p_g     = GetProperty("p_g") ;
  
  for(i = 0 ; i < nn ; i++) {
    double pc      = p_g - P_l(i) ;
    double dslsdpc = DSATURATION(pc) ;
    
    c[i] = rho_l*phi*(-dslsdpc) ;
  }
  
  return(dec) ;
}



int k1(Element_t *el,double *c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  int    nn = Element_GetNbOfNodes(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int    dec = 1 ;
  int    i ;

  /* Initialization */
  for(i = 0 ; i < nn*nn ; i++) c[i] = 0. ;
    
  /* Permeability coefficient */
  for(i = 0 ; i < nn ; i++) {
    int j ;
    for(j = i + 1 ; j < nn ; j++) {
      c[i*nn + j] = K_l ;
      c[j*nn + i] = K_l ;
    }
  }

  return(dec) ;
}



int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  int    nn = Element_GetNbOfNodes(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double x1    = Element_GetNodeCoordinate(el,1)[0] ;
  double x0    = Element_GetNodeCoordinate(el,0)[0] ;
  double dij   = x1 - x0 ;
  double dtdij = dt/dij ;
  int    dec = 1 ;
  int    i ;

  /* Initialization */
  for(i = 0 ; i < nn*nn ; i++) c[i] = 0. ;
    
  /* Permeability coefficient */
  for(i = 0 ; i < 2 ; i++) {
    double pc      = p_g - P_l(i) ;
    double dslsdpc = DSATURATION(pc) ;
    int j = 1 - i;
    
    /* Content term at node i */
    c[i*nn + i] = rho_l*phi*(-dslsdpc) ;
    
    /* Transfer term from node i to node j */
    c[i*nn + j] = dtdij*K_l ;
  }

  return(dec) ;
}
