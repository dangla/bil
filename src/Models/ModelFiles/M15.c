#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"
#include "FEM.h"

#define TITLE "Unsaturated soils"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ   (1+dim)

/* Nb of terms per point */
#define NVI   (27)
#define NVE   (1)

/* Indices of equations */
#define E_liq (dim)
#define E_mec (0)
/* Indices of unknowns */
#define I_p_l (dim)
#define I_u   (0)

/* We define some names for implicit terms */
#define M_L           (vim[0])
#define W_L           (vim + 1)
#define SIG           (vim + 4)
#define F_MASS        (vim + 13)
#define DEF_P         (vim + 16)
#define P_CO          (vim[25])
#define CRIT          (vim[26])

#define M_L_n         (vim_n[0])
#define SIG_n         (vim_n + 4)
#define DEF_P_n       (vim_n + 16)
#define P_CO_n        (vim_n[25])

/* We define some names for explicit terms */
#define K_L           (vex[0])

/* Fonctions */
static int    pm(const char *) ;
static double pie(double,double,crbe_t*) ;
static double dpiesdpl(double,double,crbe_t*) ;
static int    c15(FEM_t*,double *) ;
static int    k15(FEM_t*,double*) ;
static double rn15(double *,double *,double *,double *,Element_t*) ;
static double cr15(double *,double,double *,double *,double *,Element_t*) ;

/* Curves */
#define SATURATION_CURVE       (Element_GetCurve(el))
#define RELATIVEPERM_CURVE     (Element_GetCurve(el) + 1)
#define CAPIHARDENING_CURVE    (Element_GetCurve(el) + 2)

#define EQPRESSURE(x,y)  (pie(x,y,SATURATION_CURVE))
#define DEQPRESSURE(x,y) (dpiesdpl(x,y,SATURATION_CURVE))
#define SATURATION(x)    (Curve_ComputeValue(SATURATION_CURVE,x))
#define DSATURATION(x)   (Curve_ComputeDerivative(SATURATION_CURVE,x))
#define RELATIVEPERM(x)  (Curve_ComputeValue(RELATIVEPERM_CURVE,x))
#define CAPIHARDENING(x) (Curve_ComputeValue(CAPIHARDENING_CURVE,x))

/* Parametres */
static double gravite,kappa,mu,phi0,k_int,mu_l,rho_l,p_l0,p_g,rho_s,lambda,m,sig0,p_co0 ;

#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"kappa") == 0) return (1) ;
  else if(strcmp(s,"mu") == 0) return (2) ;
  else if(strcmp(s,"phi") == 0) return (3) ;
  else if(strcmp(s,"rho_l") == 0) return (4) ;
  else if(strcmp(s,"p_l0") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_l") == 0) return (7) ;
  else if(strcmp(s,"lambda") == 0) return (8) ;
  else if(strcmp(s,"M") == 0) return (9) ;
  else if(strcmp(s,"p_g") == 0) return (10) ;
  else if(strcmp(s,"rho_s") == 0) return (11) ;
  else if(strcmp(s,"p_co") == 0) return (12) ;
  else if(strcmp(s,"sig0") == 0) return (13) ;
  else return(-1) ;
}

int SetModelProp(Model_t *model)
/** Set the model properties 
 *  Return 0 */
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
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 14 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}


int PrintModelProp(Model_t *model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;

  printf("\n\
The system consists in (1 + dim) equations\n\
\t 1. The mass balance equation for water (liq)\n\
\t 2. The equilibrium equation (meca_1,meca_2,meca_3)\n") ;

  printf("\n\
The primary unknowns are\n\
\t 1. The liquid pressure (p_l)\n\
\t 2. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;

  fprintf(ficd,"gravite = 0    # gravity\n") ;
  fprintf(ficd,"rho_s = 2000   # mass density of the dry soil\n") ;
  fprintf(ficd,"kappa = 0.004  # slope of the elastic unloading line\n") ;
  fprintf(ficd,"lambda = 0.037 # slope of the virgin consolidation line\n") ;
  fprintf(ficd,"mu = 1.e8      # shear modulus\n") ;
  fprintf(ficd,"M = 1.2        # slope of the critical state line\n") ;
  fprintf(ficd,"p_co = 18000   # initial preconsolidation pressure\n") ;
  fprintf(ficd,"sig0 = -1000   # initial mean stress\n") ;
  fprintf(ficd,"phi = 0.25     # porosity\n") ;
  fprintf(ficd,"rho_l = 1000   # mass density of the liquid\n") ;
  fprintf(ficd,"p_l0 = 0       # initial liquid pressure\n") ;
  fprintf(ficd,"k_int = 1e-20  # intrinsic permeability\n") ;
  fprintf(ficd,"mu_l = 0.001   # liquid viscosity\n") ;
  fprintf(ficd,"p_g = 0        # gas pressure\n") ;
  fprintf(ficd,"Curves = my_file # file name p_c S_l k_rl l\n") ;
  
  return(0) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  int dim = Element_GetDimensionOfSpace(el) ;
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
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    p ;
  
  double pp_0 ;
  
  /* We can skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Donnees
  */
  gravite = GetProperty("gravite") ;
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  p_co0   = GetProperty("p_co") ;
  phi0    = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;
  sig0    = GetProperty("sig0") ;
  pp_0    = EQPRESSURE(p_l0,p_g) ;
  
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    int i ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* Pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    
    /* strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* porosity */
    double phi  = phi0 + tre ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    
    /* stresses */
    double sig[9] ;
    
    /* Compute effective stress */
    for(i = 0 ; i < 9 ; i++) sig[i] = 2*mu*eps[i] ;
    
    sig[0] += sig0 - 2*mu*tre/3. + pp_0 ;
    sig[4] += sig0 - 2*mu*tre/3. + pp_0 ;
    sig[8] += sig0 - 2*mu*tre/3. + pp_0 ;
    
    {
      /* Preconsolidation pressure  */
      double ppi = p_co0*CAPIHARDENING(pc) ;
      
      /* criterion */
      double dfsds[3][3],dgsds[3][3],hm ;
      double crit = cr15(sig,ppi,dfsds[0],dgsds[0],&hm,el) ;
      
      /* Total stresses */
      sig[0] -= pp ;
      sig[4] -= pp ;
      sig[8] -= pp ;
      
      {
        /* liquid mass */
        double m_l = rho_l*phi*sl ;
        
        /* coefficient de transfert */
        double k_l = rho_l*k_int/mu_l*RELATIVEPERM(pc) ;
        
        /* flux */
        double *gpl = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_l) ;
        double w_l[3] ;
        
        /* Compute liquid flow */
        for(i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i] ;
        w_l[dim-1] += k_l*rho_l*gravite ;
        
        /* storage in vim */
        M_L = m_l ;
        for(i = 0 ; i < 3 ; i++) W_L[i]    = w_l[i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]    = sig[i] ;
        for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
        F_MASS[dim - 1] = (rho_s + m_l)*gravite ;
        for(i = 0 ; i < 9 ; i++) DEF_P[i]  = 0. ;
        P_CO = p_co0 ;
        CRIT = crit ;
      }
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
    K_L = k_l ;
  }
  return(0) ;
}




int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/** Compute the implicit terms */
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
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
  gravite = GetProperty("gravite") ;
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  p_co0   = GetProperty("p_co") ;
  phi0    = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;
  
  
  
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vim_n += NVI , vex += NVE) {
    int i ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* Pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    double pl_n  = FEM_ComputePreviousUnknown(fem,h,nn,I_p_l) ;
    double pp_n  = EQPRESSURE(pl_n,p_g) ;
    
    /* strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* porosity */
    double phi  = phi0 + tre ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    
    /* incremental strains */
    double *deps =  FEM_ComputeIncrementalLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double trde = deps[0] + deps[4] + deps[8] ;
    
    /* stresses, plastic strains and preconsolidation pressure */
    double sig[9] ;
    {
      /* stresses */
      double sig_n[9] ;
      double sigm_n,sigm ;
      
      /* Compute effective stress */
      for(i = 0 ; i < 9 ; i++) sig_n[i] = SIG_n[i] ;
      
      sig_n[0] += pp_n ;
      sig_n[4] += pp_n ;
      sig_n[8] += pp_n ;
      
      sigm_n  = (sig_n[0] + sig_n[4] + sig_n[8])/3. ;
      sigm    = sigm_n*exp(-trde/((1. - phi0)*kappa)) ;
      
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] + 2*mu*deps[i] ;
      
      sig[0] += sigm - sigm_n - 2*mu*trde/3. ;
      sig[4] += sigm - sigm_n - 2*mu*trde/3. ;
      sig[8] += sigm - sigm_n - 2*mu*trde/3. ;
      
      /* critere */
      {
        double deps_p[9],eps_p[9] ;
        double h_c  = CAPIHARDENING(pc) ;
        double p_co = P_CO_n ;
        double ppi  = p_co*h_c ;
        double crit = rn15(sig,sig_n,&ppi,deps_p,el) ;
        
        /* plastic strains */
        for(i = 0 ; i < 9 ; i++) eps_p[i] = DEF_P_n[i] ;
        
        if(crit > 0.) {
          for(i = 0 ; i < 9 ; i++) eps_p[i] += deps_p[i] ;
          p_co = ppi/h_c ;
        }
        
        /* storage in vim */
        CRIT = crit ;
        for(i = 0 ; i < 9 ; i++) DEF_P[i]  = eps_p[i] ;
        P_CO = p_co ;
      }
      
      /* Total stresses */
      sig[0] -= pp ;
      sig[4] -= pp ;
      sig[8] -= pp ;
      
      /* storage in vim */
      for(i = 0 ; i < 9 ; i++) SIG[i]    = sig[i] ;
    }
    {
      /* liquid mass */
      double m_l = rho_l*phi*sl ;
      /* transfert coefficient  */
      double k_l = K_L ;
      /* Compute liquid flow */
      double *gpl = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_l) ;
      double w_l[3] ;
      
      for(i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i] ;
      
      w_l[dim - 1] += k_l*rho_l*gravite ;
      
      /* storage in vim */
      M_L = m_l ;
      for(i = 0 ; i < 3 ; i++) W_L[i]    = w_l[i] ;
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
      F_MASS[dim - 1] = (rho_s + m_l)*gravite ;
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
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i,j ;
  double c[IntFct_MaxNbOfIntPoints*100] ;
  double zero = 0. ;

  
  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;


  /*
  ** Poromechanics Matrix
  */
  {
    int dec = c15(fem,c) ;
    double *kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  /*
  ** Conduction Matrix
  */
  {
    int dec = k15(fem,c) ;
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_liq + i*NEQ,I_p_l + j*NEQ) += dt*kc[i*nn + j] ;
    }
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
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

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
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_L - M_L_n ;
    
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= ra[i] ;
  }
  
  /* 2.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_L,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= -dt*rf[i] ;
  }
  
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 8 ;
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
  
  /*
    input data
  */
  phi0    = GetProperty("phi") ;

  {
    /* Interpolation functions at s */
    double *h_s = Element_ComputeIsoShapeFctInActualSpace(el,s) ;
    /* pressures */
    double pl  =  FEM_ComputeCurrentUnknown(fem,h_s,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    /* saturation */
    double sl  = SATURATION(pc) ;
    /* strains */
    double eps[9] = {0,0,0,0,0,0,0,0,0} ;
    double tre = 0 ;
    double u[3] = {0,0,0} ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0},sigm ;
    double e ;
    int    i ;
      
    for(i = 0 ; i < dim ; i++) {
      u[i] = FEM_ComputeCurrentUnknown(fem,h_s,nn,I_u + i) ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      /* interpolation functions */
      double *h  = IntFct_GetFunctionAtPoint(intfct,i) ;
      double *dh = IntFct_GetFunctionGradientAtPoint(intfct,i) ;
      double *def = FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
      int j ;
      for(j = 0 ; j < dim ; j++) w_l[j] += W_L[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps[j] += def[j]/np ;
    }
    
    sigm  = (sig[0] + sig[4] + sig[8])/3 ;
    tre  = eps[0] + eps[4] + eps[8] ;
    e = tre/(1. - phi0) ;
      
    i = 0 ;
    Result_Store(r + i++,&pl,"pression-liquide",1) ;
    Result_Store(r + i++,u,"deplacements",3) ;
    Result_Store(r + i++,w_l,"flux-liquide",3) ;
    Result_Store(r + i++,sig,"contraintes",9) ;
    Result_Store(r + i++,&sl,"saturation",1) ;
    Result_Store(r + i++,&pp,"pression-pi",1) ;
    Result_Store(r + i++,&e,"indice-des-vides",1) ;
    Result_Store(r + i++,&sigm,"sig_m",1) ;
  }
  
  return(NbOfOutputs) ;
}

double pie(double pl,double pg,crbe_t *cb)
{
  int    i ;
  double pc,sl,sg,u ;
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double pc1 = Curve_GetXRange(cb)[0] ;
  double pc2 = Curve_GetXRange(cb)[1] ;
  double *sat = Curve_GetYValue(cb) ;
  double dpc = (pc2 - pc1)/n_i ;
  double zero = 0.,un = 1. ;

  pc  = pg - pl ;
  
  /* U */
  if(pc < pc1) {
    u = zero ;
    sl = sat[0] ;
    
  } else {
    if(pc > pc2) pc = pc2 ;
    u  = zero ;
    for(i = 0 ; pc1 + (i + 1)*dpc < pc ; i++) u += sat[i] + sat[i+1] ;
    u *= dpc*0.5 ;
    u += sat[0]*pc1 ;
    sl  = sat[i] + (sat[i+1] - sat[i])/dpc*(pc - pc1 - i*dpc) ;
    u += (sat[i] + sl)*0.5*(pc - pc1 - i*dpc) - sl*pc ;
  }
  
  /* pi */
  sg  = un - sl ;
  return(sl*pl + sg*pg - 2./3.*u) ;
}


double dpiesdpl(double pl,double pg,crbe_t *cb)
{
  double pc  = pg - pl ;
  double sl  = Curve_ComputeValue(cb,pc) ;
  double dslsdpc = Curve_ComputeDerivative(cb,pc) ;
  double dusdpc  = - pc*dslsdpc ;
  
  return(sl - dusdpc/3) ;
}


int c15(FEM_t *fem,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  Element_t *el = FEM_GetElement(fem) ;
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  
  int    dec = 100 ;
  int    p ;
  double zero = 0. ;
  
  /*
    Input data
  */
  gravite = GetProperty("gravite") ;
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  p_co0   = GetProperty("p_co") ;
  phi0    = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l   = GetProperty("rho_l") ;
  rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  p_g     = GetProperty("p_g") ;
  
  for(p = 0 ; p < np ; p++ , vim += NVI , vex += NVE) {
    int i,j ;
    double *c1 = c + p*dec ;
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    
    /* pressures */
    double pl  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = EQPRESSURE(pl,p_g) ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    double dslsdpc = DSATURATION(pc) ;
    
    /* preconsolidation pressure */
    double h_c  = CAPIHARDENING(pc) ;
    double ppi  = P_CO*h_c ;
    
    /* effective stresses */
    double sig[9] ;
    
    for(i = 0 ; i < 9 ; i++) sig[i] = SIG[i] ;
    
    sig[0] += pp ;
    sig[4] += pp ;
    sig[8] += pp ;

    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    { /* Mechanics */
      double sigm  = (sig[0] + sig[4] + sig[8])/3. ;
      double lame = -sigm/((1. - phi0)*kappa) - 2*mu/3. ;
      double crit = CRIT ;
      
      for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
        C1(i,i,j,j) += lame ;
        C1(i,j,i,j) += mu ;
        C1(i,j,j,i) += mu ;
      }
      
      /* criterion */
      if(crit >= 0.) {
        double dfsds[3][3],dgsds[3][3],fc[3][3],cg[3][3],fcg ;
        double hm ;
        
        cr15(sig,ppi,dfsds[0],dgsds[0],&hm,el) ;
        
        fcg = hm ;
        for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
          int k,l ;
          
          fc[i][j] = 0. ;
          cg[i][j] = 0. ;
          
          for(k = 0 ; k < 3 ; k++) for(l = 0 ; l < 3 ; l++) {
            fc[i][j] += dfsds[k][l]*C1(l,k,i,j) ;
            cg[i][j] += C1(i,j,k,l)*dgsds[l][k] ;
            fcg += dfsds[j][i]*C1(i,j,k,l)*dgsds[l][k] ;
          }
        }
        
        fcg = 1./fcg ;
        
        for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
          int k,l ;
          
          for(k = 0 ; k < 3 ; k++) for(l = 0 ; l < 3 ; l++) {
            C1(i,j,k,l) -= cg[i][j]*fc[k][l]*fcg ;
          }
        }
      }
      
      c1 += 81 ;
      for(i = 0 ; i < 3 ; i++) B1(i,i) = - DEQPRESSURE(pl,p_g) ;
    }
    
    c1 += 9 ;
    /* Hydraulic */
    for(i = 0 ; i < 3 ; i++) B1(i,i) = rho_l*sl ;
    c1 += 9 ;
    c1[0] = -rho_l*phi0*dslsdpc ;
  }
  return(dec) ;
  
#undef C1
#undef B1
}



int k15(FEM_t *fem,double *c)
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
    c1[0] = K_L ;
    c1[4] = K_L ;
    c1[8] = K_L ;
  }
  return(dec) ;
}


double cr15(double *sig,double pc,double *dfsds,double *dgsds,double *hm,Element_t *el)
/* Critere de Cam-Clay */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit,m2,v ;
  int    i ;
  /*
    Donnees
  */
  kappa   = GetProperty("kappa") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  phi0    = GetProperty("phi") ;
  m2      = m*m ;
  v       = 1./(lambda - kappa) ;
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Les gradients
  */
  for(i = 0 ; i < 9 ; i++) {
    double dev = sig[i] - p*id[i] ;
    
    dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    dgsds[i] = dfsds[i] ;
  }
  
  /* Le module d'ecrouissage */
  *hm = v/(1 - phi0)*p*(2*p + pc)*pc ;
  return(crit) ;
}


double rn15(double *sig,double *sig_n,double *p_co,double *eps_p,Element_t *el)
/* Critere de Cam-Clay : return mapping algorithm (sig,p_co,eps_p)
   (d apres Borja & Lee 1990, modifie par Dangla)
*/
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double tol = 1.e-8 ;
  double p_t,q_t,p,q,pc,crit,m2,v,a ;
  double dl ;
  int    i ;
  /*
    Donnees
  */
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  phi0    = GetProperty("phi") ;
  m2      = m*m ;
  v       = 1./(lambda - kappa) ;
  
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  pc   = *p_co ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Algorithme de projection (closest point projection)
    Une seule boucle iterative pour le calcul de p, racine de
    q*q/m2 + p*(p + pc) = 0
    Les autres variables (pc,q,dl) sont donnees explicitement par p.
  */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double dlsdp  = ((1 - phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*dlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Les variables (pc,dl,q) sont explicites en p */
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
        printf("pas de convergence (rn15)") ;
        exit(0) ;
      }
    }
  }
  
  /*
    Les contraintes et deformations plastiques
  */
  a = 1./(1 + 6*mu/m2*dl) ;
  
  for(i = 0 ; i < 9 ; i++) {
    double dev      = a*(sig[i] - p_t*id[i]) ;
    double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
    sig[i]   = p*id[i] + dev ;
    eps_p[i] = dl*dfsds ;
  }
  
  /* La pression de consolidation */
  *p_co = pc ;
  return(crit) ;
}

