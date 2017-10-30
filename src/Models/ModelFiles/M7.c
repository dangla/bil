#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Common.h"

/* The Finite Element Method */
#include "FEM.h"

#define TITLE   "Unsaturated Poroelasticity"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ   (1+dim)

#define NVI   (17)
#define NVE   (1)
#define NV0   (0)

#define E_liq (dim)
#define E_mec (0)

#define I_p_l (dim)
#define I_u   (0)

#define M_l     (vim[0])
#define W_l     (vim + 1)
#define SIG     (vim + 4)
#define F_MASS  (vim + 13)
#define PHI     (vim[16])

#define M_ln    (vim_n[0])

#define K_l     (vex[0])


/* Functions */
static int    pm(const char*) ;
static double EquivalentPressure(double,double,Curve_t*) ;
static int    c7(FEM_t*,double*) ;
static int    k7(FEM_t*,double*) ;
/* Curves */
#define EQPRESSURE(x,y)  (EquivalentPressure(x,y,Element_GetCurve(el)))
#define SATURATION(x)    (Curve_ComputeValue(Element_GetCurve(el),x))
#define DSATURATION(x)   (Curve_ComputeDerivative(Element_GetCurve(el),x))
#define RELATIVEPERM(x)  (Curve_ComputeValue(Element_GetCurve(el) + 1,x))

/* Parameters */
static double gravite ;
static double young,poisson,b,N ;
static double k_int,mu_l ;
static double phi0,rho_l,rho_s ;
static double p_g,p_l0,sig0_11,sig0_22,sig0_33 ;

#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"young") == 0) return (1) ;
  else if(strcmp(s,"poisson") == 0) return (2) ;
  else if(strcmp(s,"phi") == 0) return (3) ;
  else if(strcmp(s,"rho_l") == 0) return (4) ;
  else if(strcmp(s,"p_l0") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_l") == 0) return (7) ;
  else if(strcmp(s,"b") == 0) return (8) ;
  else if(strcmp(s,"N") == 0) return (9) ;
  else if(strcmp(s,"p_g") == 0) return (10) ;
  else if(strcmp(s,"rho_s") == 0) return (11) ;
  else if(strcmp(s,"sig0_11") == 0) return (12) ;
  else if(strcmp(s,"sig0_22") == 0) return (13) ;
  else if(strcmp(s,"sig0_33") == 0) return (14) ;
  else return(-1) ;
}

int SetModelProp(Model_t *model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_liq,"liq") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn[i]) ;
  }
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,I_p_l,"p_l") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,I_u + i,name_unk[i]) ;
  }
  
  return(0) ;
}

int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/** Read the material properties in the stream file ficd */
{
  int  NbOfProp = 15 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/** Print the model characteristics */
{
  printf(TITLE) ;
  
  if(!ficd) return(0) ;

  printf("\n\n\
The system of equations consists in 1+dim equations:\n\
\t 1. the water mass balance equation (p_l)\n\
\t 2. the mechanical equilibrium (u_1,u_2,u_3)\n") ;

  printf("\n\
Input data examples\n\n") ;

  fprintf(ficd,"gravite = 0       # gravite\n") ;
  fprintf(ficd,"rho_s = 0         # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 0.833333  # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.25    # coefficient de Poisson\n") ;
  fprintf(ficd,"phi = 0.3         # porosite\n") ;
  fprintf(ficd,"rho_l = 1         # masse volumique du fluide\n") ;
  fprintf(ficd,"p_l0 = 0          # pression initiale du fluid\n") ;
  fprintf(ficd,"k_int = 1         # permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 1          # viscosite du liquide\n") ;
  fprintf(ficd,"p_g = 0           # pression de gaz\n") ;
  fprintf(ficd,"b = 1             # coefficient de Biot\n") ;
  fprintf(ficd,"N = 0             # compressibilite des pores\n") ;
  fprintf(ficd,"sig0_11 = 0       # contrainte initiale 11\n") ;
  fprintf(ficd,"sig0_22 = 0       # contrainte initiale 22\n") ;
  fprintf(ficd,"sig0_33 = 0       # contrainte initiale 33\n") ;
  fprintf(ficd,"courbes = My_file # Nom du fichier : p_c S_l k_rl\n") ;

  return(0) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
/** Define some properties attached to each element */
{
  IntFct_t *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t *fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;

  {
    double *r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
  
    /* hydraulic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_liq) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
      
    /* other */
    } else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 */ 
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double dmu,lame,pp0 ;
  int    p ;
  double un = 1.,deux = 2. ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  gravite = GetProperty("gravite") ;
  young   = GetProperty("young") ;
  poisson = GetProperty("poisson") ;
  phi0    = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;
  b       = GetProperty("b") ;
  N       = GetProperty("N") ;
  sig0_11 = GetProperty("sig0_11") ;
  sig0_22 = GetProperty("sig0_22") ;
  sig0_33 = GetProperty("sig0_33") ;

  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  pp0     = EQPRESSURE(p_l0,p_g) ;

  /* Loop on the integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    int i ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    
    /* strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
      
    /* stresses */
    double sig[9] ;
    
    for(i = 0 ; i < 9 ; i++) sig[i] = dmu*eps[i] ;
    
    sig[0] += sig0_11 + lame*tre - b*(pp - pp0) ;
    sig[4] += sig0_22 + lame*tre - b*(pp - pp0) ;
    sig[8] += sig0_33 + lame*tre - b*(pp - pp0) ;
    
    {
      /* porosite */
      double dphi = b*tre + N*(pp - pp0) ;
      double phi  = phi0 + dphi ;
      
      /* saturation */
      double sl  = SATURATION(pc) ;
      
      /* liquid mass */
      double m_l = rho_l*phi*sl ;
      
      /* transfert coefficient */
      double k_l = rho_l*k_int/mu_l*RELATIVEPERM(pc) ;
      
      /* flux */
      double *gpl = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_l) ;
      double w_l[3] ;
      
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i] ;
      w_l[dim-1] += k_l*rho_l*gravite ;
      
      /* storage in vim */
      M_l = m_l ;
      for(i = 0 ; i < 3 ; i++) W_l[i]    = w_l[i] ;
      for(i = 0 ; i < 9 ; i++) SIG[i]    = sig[i] ;
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
      F_MASS[dim-1] = (rho_s+m_l)*gravite ;
      PHI = phi ;
    }
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/** Compute the explicit terms */
{
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vex += NVE) {
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    
    /* pressures */
    double pl  = FEM_ComputePreviousUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    
    /* transfert coefficient */
    double k_l = rho_l*k_int/mu_l*RELATIVEPERM(pc) ;
    
    /* storage in vex */
    K_l = k_l ;
  }
  
  return(0) ;
}




int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/** Compute the implicit terms */
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double dmu,lame,pp0 ;
  int    p ;
  double un = 1.,deux = 2. ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  gravite = GetProperty("gravite") ;
  young   = GetProperty("young") ;
  poisson = GetProperty("poisson") ;
  phi0    = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;
  b       = GetProperty("b") ;
  N       = GetProperty("N") ;
  sig0_11 = GetProperty("sig0_11") ;
  sig0_22 = GetProperty("sig0_22") ;
  sig0_33 = GetProperty("sig0_33") ;

  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  pp0     = EQPRESSURE(p_l0,p_g) ;


  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    int i ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    
    /* strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
      
    /* stresses */
    double sig[9] ;
      
    for(i = 0 ; i < 9 ; i++) sig[i] = dmu*eps[i] ;
    
    sig[0] += sig0_11 + lame*tre - b*(pp - pp0) ;
    sig[4] += sig0_22 + lame*tre - b*(pp - pp0) ;
    sig[8] += sig0_33 + lame*tre - b*(pp - pp0) ;
      
    {
      /* porosity */
      double dphi = b*tre + N*(pp - pp0) ;
      double phi  = phi0 + dphi ;
        
      /* saturation */
      double sl  = SATURATION(pc) ;
        
      /* liquid mass */
      double m_l = rho_l*phi*sl ;
        
      /* transfert coefficient */
      double k_l = K_l ;
        
      /* flux */
      double *gpl = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_l) ;
      double w_l[3] ;
        
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i] ;
      w_l[dim-1] += k_l*rho_l*gravite ;
        
      /* storage in vim */
      M_l = m_l ;
      for(i = 0 ; i < 3 ; i++) W_l[i] = w_l[i] ;
      for(i = 0 ; i < 9 ; i++) SIG[i] = sig[i] ;
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
      F_MASS[dim-1] = (rho_s+m_l)*gravite ;
      PHI = phi ;
    }
  }

  return(0) ;
}




int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  char *method = Material_GetMethod(Element_GetMaterial(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i,j,dec ;
  double c[MAX_PGAUSS*100] ;
  double zero = 0. ;

  
  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
  ** Poromechanic Matrix
  */
  dec = c7(fem,c) ;
  {
    double *kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  
  /*
  ** Conduction Matrix
  */
  dec = k7(fem,c) ;
  {
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_liq + i*NEQ,I_p_l + j*NEQ) += dt*kc[i*nn + j] ;
    }
  }
  
  /* elements P2P1 */
  if(strstr(method,"P2P1")) {
    FEM_TransformMatrixFromDegree2IntoDegree1(fem,I_p_l,E_liq,k) ;
  }
  return(0) ;
#undef K
}




int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  char *method = Material_GetMethod(Element_GetMaterial(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < nn*NEQ ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /* 1. Mechanics */
  /* 1.1 Stresses */
  {
    double *vim = vim_1 ;
    double *rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }
  }
  
  /* 1.2 Body forces */
  {
    double *vim = vim_1 ;
    double *rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_mec + dim - 1) -= -rbf[i] ;
  }
  
  /* 2. Hydraulic */
  /* 2.1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    double *ra ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_l - M_ln ;
    
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= ra[i] ;
  }
  
  /* 2.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_l,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= -dt*rf[i] ;
  }

  /* elements P2P1 */
  if(strstr(method,"P2P1")) {
    FEM_TransformResiduFromDegree2IntoDegree1(fem,E_liq,r) ;
  }
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 7 ;
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double zero = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /* initialization */
  {
    int    i,j ;
    for(i = 0 ; i < NbOfOutputs ; i++) for(j = 0 ; j < 9 ; j++) {
      Result_GetValue(r + i)[j] = zero ;
    }
  }

  {
    /* Interpolation functions at s */
    double *h_s = FEM_ComputeIsoShapeFctInActualSpace(fem,s) ;
    
    /* pressures */
    double pl  =  FEM_ComputeCurrentUnknown(fem,h_s,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    double u[3] = {0,0,0} ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double phi = 0 ;
    int    i ;
      
    for(i = 0 ; i < dim ; i++) {
      u[i] = FEM_ComputeCurrentUnknown(fem,h_s,nn,I_u + i) ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      int j ;
      for(j = 0 ; j < dim ; j++) w_l[j] += W_l[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
	
      phi += PHI/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&pl,"pression-liquide",1) ;
    Result_Store(r + i++,u,"deplacements",3) ;
    Result_Store(r + i++,w_l,"flux-liquide",3) ;
    Result_Store(r + i++,sig,"contraintes",9) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&sl,"saturation",1) ;
    Result_Store(r + i++,&pp,"pression-pi",1) ;
      
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }

  return (NbOfOutputs) ;
}

double EquivalentPressure(double pl,double pg,Curve_t *cb)
{
  int    j ;
  double pc,sl,sg,u ;
  int    np = 3 ;
  double a[3] = {0.93246951420,0.66120938646,0.23861918608} ;
  double w[3] = {0.17132449237,0.36076157304,0.46791393457} ;

  pc = pg - pl ;
  sl = Curve_ComputeValue(cb,pc) ;
  sg = 1. - sl ;
  u  = 0. ;
  for(j = 0 ; j < np ; j++) {
    u += w[j]*(Curve_ComputeValue(cb,pc*0.5*(1 + a[j]))
       +       Curve_ComputeValue(cb,pc*0.5*(1 - a[j]))) ;
  }
  u  = u*pc*0.5 - sl*pc ;
  return(sl*pl + sg*pg - 2./3.*u) ;
}

int c7(FEM_t *fem,double *c)
/*
**  Poroelastic matrix (c) and shift (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  Element_t *el = FEM_GetElement(fem) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  double dmu,lame,mu ;
  int    dec = 100 ;
  int    p ;
  double zero = 0.,un = 1.,deux = 2. ;

  /*
    Input data
  */
  young   = GetProperty("young") ;
  poisson = GetProperty("poisson") ;
  phi0    = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  p_g     = GetProperty("p_g") ;
  b       = GetProperty("b") ;
  N       = GetProperty("N") ;

  dmu     = young/(un+poisson) ;
  mu      = dmu/deux ;
  lame    = dmu*poisson/(un-deux*poisson) ;

  for(p = 0 ; p < np ; p++) {
    int i,j ;
    double *c1 = c + p*dec ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    
    /* pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    double dslsdpc = DSATURATION(pc) ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;

    /* Mechanics */
    for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
      C1(i,i,j,j) += lame ;
      C1(i,j,i,j) += mu ;
      C1(i,j,j,i) += mu ;
    }
    c1 += 81 ;
    for(i = 0 ; i < 3 ; i++) B1(i,i) = -b*(sl + pc*dslsdpc/3) ;

    c1 += 9 ;
    /* Hydraulic */
    for(i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*sl*b ;
    c1 += 9 ;
    c1[0] = -rho_l*phi0*dslsdpc + rho_l*sl*N*(sl + pc*dslsdpc/3) ;
  }
  
  return(dec) ;

#undef C1
#undef B1
}

int k7(FEM_t *fem,double *c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  double zero = 0. ;

  for(p = 0 ; p < np ; p++ , vex += NVE) {
    int i ;
    double *c1 = c + p*dec ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */
    c1[0] = K_l ;
    c1[4] = K_l ;
    c1[8] = K_l ;
  }
  
  return(dec) ;
}
