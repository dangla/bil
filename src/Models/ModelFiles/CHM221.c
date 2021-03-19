#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#include "FEM.h"

#define TITLE   "Chemo-hydro-mechanical model"
#define AUTHORS "G.MELOT"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ     (2 + dim)
/* Nb of implicit terms */
#define NVI     (52)
/* Nb of explicit terms */
#define NVE     (5)
/* Nb of constant terms */
#define NV0     (0)

/* Equation index */
#define E_mec     (0)
#define E_mass    (dim)
#define E_salt    (1 + dim)


/* Unknown index */
#define I_u     (0)
#define I_p     (dim)
#define I_w_s   (1 + dim)


/* Implicit terms index */
#define W_l		(vim)			// Liquid flux (water + dissolved salt)
#define W_salt	(vim + 3)		// Salt flux (dissolved salt + crystal salt)

#define M_total	(vim[6])		// Mass of (water + dissolved salt + crystal salt)
#define M_salt	(vim[7])		// Mass of dissolved salt + crystal salt

#define Phi_c   (vim[8])		// Volumetric fraction of crystal salt
#define Phi_l 	(vim[9])		// Porosity
#define Rho_l 	(vim[10])		// Density of liquide phases ( not equal to Rho_w)
#define Rho_w	(vim[11])		// Density of water
#define Dm_s	(vim[12])		// Quantity of crystal salt dissolved

#define SIG	    (vim + 13)
#define EPSI_e	(vim + 22)
#define EPSI_vp	(vim + 31)

#define M_W_water (vim[40])		// Mass of water transported through each element in each time step
#define M_W_salt  (vim[41])		// Mass of salt transported through each element in each time step

#define SIG_M	(vim[42])
#define SIG_D   (vim + 43)


/* Implicite variables at the time n */
#define M_totaln  	(vim_n[6])
#define M_saltn  	(vim_n[7])
#define Phi_c_n     (vim_n[8])
#define Phi_ln 	    (vim_n[9])
#define Rho_ln 	    (vim_n[10])

#define SIG_n 		(vim_n + 13)
#define EPSI_vpn	(vim_n + 31)

#define M_W_watern 	(vim_n[40])
#define M_W_saltn 	(vim_n[41])

#define SIG_M_n		(vim_n[42])
#define SIG_D_n   	(vim_n + 43)

/* Explicit terms index */
#define K_l	(vex[0])
#define D	(vex[1])
#define Tau	(vex[2])

#define M_Wc_water (vex[3])		// Mass of water transported through each element from the begining
#define M_Wc_salt  (vex[4])		// Mass of salt transported through each element from the begining
////////////////////////
/* Physical constant */
#define RT        (2436.)   /* product of R=8.3143 and T=293 (Pa.m3/mol) */ 

/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

/* Functions */
static int     pm(const char*) ;
static void    GetProperties(Element_t*) ;

static int    kp2(FEM_t*,double*,double) ;
static int    klp2(FEM_t*,double*) ;
static int    klws2(FEM_t*,double*) ;
static int    ksp2(FEM_t*,double*) ;
static int    ksws2(FEM_t*,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,int) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static void  ComputeSecondaryVariables(Element_t*,double,double*) ;

/************************/
/** Material parameters */
/************************/
/**********************************************************************/
//Transport parameters
static double  k_i;				//Initial intrinsic permeability
static double  d_i;				//Initial diffusion coefficient
static double  tau_i;			//Initial efficiency coefficient
static double  phi_l_i;			//Initial porosity
//Water property
static double  rho_w_i;			//Liquid density
static double  viscosity_w;		//Dynamic viscosity of water		
static double  compressibility_w;	//Compressibility of water
//Salt property
static double  rho_c;			//Crystal salt density
static double  beta;			//Dissolution rate
static double  sigma_c;			//Crystal specific surface	
static double  Molarm_s;		//Solute molar mass
static double  w_s_sat;			//Solute mass fraction at saturation
static double  v_m_s;			//Specific volume of disolve salt
//Elasticity law parameters
static double  young;			//Young modulus
static double  poisson;			//Poisson coefficient
//Creep law parameters 
static double  eta_v;			//Bitumen volumetric viscosity (Maxwell law)
static double  eta_d;			//Bitumen deviatoric viscosity (Maxwell law)
/**********************************************************************/
/** Initial condition */
/**********************************************************************/
static double sig0_11;			//Initial longitudinal stress
static double sig0_22;			//Initial stress
static double sig0_33;			//Initial stress 
static double w_s_i;			//Initial solute mass fraction
static double phi_c_i;			//Initial crystal volumetric fraction
static double p_i;				//Initial pore water pressure
/**********************************************************************/
static double b;				//Biot coefficient
/**********************************************************************/
/** Coefficients dependency  */
/**********************************************************************/
static double D_dep;			// 0 -> D = cte / 1 -> D = f(phi_l) 
static double K_dep;			// 0 -> K = cte / 1 -> K = f(phi_l) 
static double Tau_dep;			// 0 -> Tau = cte / 1 -> Tau = f(phi_l)
static double rho_l_dep;			// 0 -> rho_l = cte = p / 1 -> rho_l = f(w_s)
/**********************************************************************/
static double SampleSurface;	//Sample exchange surface [m²]
static double ElementSize; 		//Element size [m]  

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) // creat the table of all variables


#define NbOfVariables     (95)
static double Variable[NbOfVariables] ;

#define I_U            	(0)
#define I_P        	   	(3)
#define I_P_n			(4)		
#define I_W_S          	(5)
#define I_W_S_n			(6)

#define I_EPS          	(7)
#define I_EPS_n			(16)

#define I_SIG          	(25)
#define I_SIG_n			(34)
#define I_SIG_M			(43)
#define I_SIG_M_n		(44)
#define I_SIG_D			(45)
#define I_SIG_D_n		(54)

#define I_PHI_C			(63)
#define I_PHI_C_n		(64)
#define I_PHI_L			(65)

#define I_RHO_W			(66)
#define I_RHO_L			(67)

#define I_M_TOT			(68)
#define I_M_SALT		(69)

#define I_K_L			(70)
#define I_D				(71)
#define I_TAU			(72)

#define I_GRD_P      	(73)
#define I_GRD_W_S      	(76)

#define I_W_F          	(79)
#define I_W_D          	(82)
#define I_W_O          	(85)
#define I_W_L          	(88)
#define I_W_SALT       	(91)
#define I_RHO_L_n       (94)








////////////////////////////////////////////////////////////////////////

int pm(const char* s)	// this function is used to read model parameters 
{
       if(!strcmp(s,"k_i")) 		return(0) ;//if s egal "gravity" then return 0
  else if(!strcmp(s,"d_i"))   		return(1) ;//...	
  else if(!strcmp(s,"tau_i")) 		return(2) ;//...
  else if(!strcmp(s,"phi_l_i"))   	return(3) ;//...
  else if(!strcmp(s,"rho_w_i"))  	return(4) ;//...	
  else if(!strcmp(s,"viscosity_w"))	return(5) ;//...
  else if(!strcmp(s,"compressibility_w"))return(6) ;//...
  else if(!strcmp(s,"rho_c"))  		return(7) ;//...	
  else if(!strcmp(s,"beta")) 		return(8) ;//...
  else if(!strcmp(s,"sigma_c"))   	return(9) ;//...
  else if(!strcmp(s,"Molarm_s"))   	return(10);//...	
  else if(!strcmp(s,"w_s_sat"))		return(11);//...
  else if(!strcmp(s,"v_m_s"))		return(12);//...
  else if(!strcmp(s,"young"))	  	return(13);//...
  else if(!strcmp(s,"poisson"))		return(14);//...	
  else if(!strcmp(s,"eta_v"))		return(15);//...
  else if(!strcmp(s,"eta_d"))		return(16);//...
  else if(!strcmp(s,"sig0_11"))	   	return(17);//...
  else if(!strcmp(s,"sig0_22"))	   	return(18);//...
  else if(!strcmp(s,"sig0_33"))	   	return(19);//...
  else if(!strcmp(s,"w_s_i"))	   	return(20);//...
  else if(!strcmp(s,"phi_c_i"))	   	return(21);//...
  else if(!strcmp(s,"p_i"))	   		return(22);//...
  else if(!strcmp(s,"b"))	   		return(23);//...
  else if(!strcmp(s,"D_dep"))	   	return(24);//...
  else if(!strcmp(s,"K_dep"))	   	return(25);//...
  else if(!strcmp(s,"Tau_dep"))	   	return(26);//...
  else if(!strcmp(s,"rho_l_dep"))	return(27);//... 
  else if(!strcmp(s,"SampleSurface"))return(28);//...
  else if(!strcmp(s,"ElementSize")) return(29);//... 
  else return(-1);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void GetProperties(Element_t* el)// copy the value of material properties
{
  k_i   		= GetProperty("k_i");
  d_i  			= GetProperty("d_i");
  tau_i			= GetProperty("tau_i");
  phi_l_i		= GetProperty("phi_l_i");
  rho_w_i		= GetProperty("rho_w_i");
  viscosity_w	= GetProperty("viscosity_w");
  compressibility_w	= GetProperty("compressibility_w");
  rho_c			= GetProperty("rho_c");
  beta			= GetProperty("beta");
  sigma_c		= GetProperty("sigma_c");
  Molarm_s		= GetProperty("Molarm_s");
  w_s_sat		= GetProperty("w_s_sat");
  v_m_s			= GetProperty("v_m_s");
  young			= GetProperty("young");
  poisson		= GetProperty("poisson");
  eta_v			= GetProperty("eta_v");
  eta_d			= GetProperty("eta_d");
  sig0_11		= GetProperty("sig0_11");
  sig0_22		= GetProperty("sig0_22");
  sig0_33		= GetProperty("sig0_33");
  w_s_i			= GetProperty("w_s_i");
  phi_c_i		= GetProperty("phi_c_i");
  p_i			= GetProperty("p_i");
  b             = GetProperty("b");
  D_dep         = GetProperty("D_dep");
  K_dep         = GetProperty("K_dep");
  Tau_dep       = GetProperty("Tau_dep");
  rho_l_dep     = GetProperty("rho_l_dep");
  SampleSurface	= GetProperty("SampleSurface");
  ElementSize	= GetProperty("ElementSize");
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int SetModelProp(Model_t *model)
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */  
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn[i]) ;
  }
  Model_CopyNameOfEquation(model,E_mass,"Mass") ;
  Model_CopyNameOfEquation(model,E_salt,"Salt") ;

  /** Names of the main (nodal) unknowns */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,I_u + i,name_unk[i]) ;
  }
  Model_CopyNameOfUnknown(model,I_p,"p") ;
  Model_CopyNameOfUnknown(model,I_w_s,"w_s") ;

  
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Reading of material properties in file ficd */
{
  int  NbOfProp = 30;
  Material_ScanProperties(mat,datafile,pm) ;

  return(NbOfProp) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Following it is the informations given when you ask about this model with the commande "bil -model CHM2" 
int PrintModelChar(Model_t* model,FILE* ficd)
/** Print the model characteristics */
{
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mechanical Equilibrium		(mec)\n") ;
  printf("\t- Total mass conservation		(mass)\n") ;
  printf("\t- Salt mass conservation		(salt)\n") ;


  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Displacements					(u_1,u_2,u_3) \n") ;
  printf("\t- Pore water pressure			(p) \n") ;
  printf("\t- Solute mass fraction			(w_s) \n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"k_i = 4e-25 			# (m2)    Initial intrinsic permeability\n") ;
  fprintf(ficd,"d_i = 3e-14 			# (m2/s)  Initial diffusion coefficient\n") ;
  fprintf(ficd,"tau_i = 0.9996  			# (-)     Initial efficiency coefficient\n") ;
  fprintf(ficd,"phi_l_i = 0.01   		# (-)     Initial porosity\n") ;
  fprintf(ficd,"rho_w_i = 1000  		# (kg/m3) Water density \n") ;
  fprintf(ficd,"viscosity_w = 1e-3 		# (Pa.s)  Dynamic viscosity of water\n") ;
  fprintf(ficd,"compressibility_w = 5e-10 	# (1/Pa)  Compressibility of water\n") ;
  fprintf(ficd,"rho_c = 2260  			# (kg/m3) Crystal salt density\n") ;
  fprintf(ficd,"beta = 1e-5  			# (kg/(s.m3)) Dissolution rate\n") ;
  fprintf(ficd,"sigma_c = 148   		# (m2/m3) Crystal specific surface\n") ;
  fprintf(ficd,"Molarm_s = 85e-3   		# (kg/mol) Solute molar mass\n") ;
  fprintf(ficd,"w_s_sat = 0.47  		# (-)     Solute mass fraction at saturation\n") ;
  fprintf(ficd,"v_m_s = 1.689e-3		# (m3/kgl) Specific volume of dissolved salt\n") ;
  fprintf(ficd,"young = 1e8   			# (Pa)    Young modulus\n") ;
  fprintf(ficd,"poisson = 0.33  		# (-)     Poisson's coefficient\n") ;
  fprintf(ficd,"eta_v = 3.33e14			# (Pa.s)  Bitumen volumetric viscosity (Maxwell law)\n") ;
  fprintf(ficd,"eta_d = 1e15			# (Pa.s)  Bitumen deviatoric viscosity (Maxwell law)\n") ;
  fprintf(ficd,"sig0_11 = 0     		# (Pa)    Initial longitudinal stress\n") ;
  fprintf(ficd,"sig0_22 = 0    			# (Pa)    Initial stress\n") ;
  fprintf(ficd,"sig0_33 = 0     		# (Pa)    Initial stress\n") ;
  fprintf(ficd,"w_s_i = 0.47  			# (-)     Initial solute mass fraction\n") ;
  fprintf(ficd,"phi_c_i = 0.16  		# (-)     Initial crystal volumetric fraction\n") ;
  fprintf(ficd,"p_i = 1e5   			# (Pa)    Initial pore water pressure\n") ;
  fprintf(ficd,"b = 1     			# (-)     Coefficient de Biot\n") ;
  fprintf(ficd,"D_dep = 2     			# (-)     Diffusion coef. dependency with phi_lt\n") ;
  fprintf(ficd,"K_dep = 2     			# (-)     Permeability dependency with phi_l\n") ;
  fprintf(ficd,"Tau_dep = 0     		# (-)     Osmosis coef. dependency with phi_l\n") ;
  fprintf(ficd,"rho_l_dep = 1     		# (-)     Liquid density dependency with w_s\n") ;
  fprintf(ficd,"SampleSurface = 19.63e-4 # (m2)    Sample exchange surface\n") ;
  fprintf(ficd,"ElementSize = 2e-5 		# (m)    Element size\n") ;
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 
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
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_mass) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }   
    /* Salt */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_salt) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
    // other //  
     else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		INITIAL STATE		*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int ComputeInitialState(Element_t* el)
{
/** Compute the initial state i.e. the constant terms,the explicit terms,the implicit terms. */
  double* vim0 = Element_GetImplicitTerm(el) ;
  double *vex0 = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    m ;
  // à garder ???////////////////////
  //int nn = IntFct_GetNbOfNodes(intfct) ;
  //int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  //FEM_t *fem = FEM_GetInstance(el) ;
 /////////////////////////////
 
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;

  /* Pre-initialization */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
	  
	  /* Storage in vim */
	  {
		  double* vim = vim0 + m*NVI ;
		  int i ;
		  
		  for(i = 0 ; i < 9 ; i++) SIG[i] = 0 ;
		  SIG[0] = sig0_11 ;
		  SIG[4] = sig0_22 ;
		  SIG[8] = sig0_33 ;
		  
		  Phi_c = phi_c_i ;
		  
		  SIG_M = (1./3.)* (SIG[0]+SIG[4]+SIG[8]) ;
		  for(i = 0 ; i < 9 ; i++) SIG_D[i] = SIG[i] ;
		  SIG_D[0] += - SIG_M ;
		  SIG_D[4] += - SIG_M ;
		  SIG_D[8] += - SIG_M ; 
	  }
  }
  
  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,0,m) ;

    /* Storage in vim */
    {
		double* vim = vim0 + m*NVI ;
		int i ;
		
		M_total	= x[I_M_TOT] ;
		M_salt	= x[I_M_SALT] ;
		
		for(i = 0 ; i < 3 ; i++) W_l[i] = x[I_W_L + i] ;
		for(i = 0 ; i < 3 ; i++) W_salt[i] = x[I_W_SALT + i] ;
		
		Rho_l = x[I_RHO_L] ;		
		Rho_w = x[I_RHO_W] ;		
		// manque Dm_s ... utile ??
		
		for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
		for(i = 0 ; i < 9 ; i++) SIG_D[i] = x[I_SIG_D + i] ;
		SIG_M = x[I_SIG_M] ;
		
		Phi_l = x[I_PHI_L] ;
		Phi_c = x[I_PHI_C] ;
		
		/* For post-treatment of M_salt_leached and m_water_absorbed */
		M_W_water = 0 ;
		M_W_salt = 0 ;	
	}
	
    /* storage in vex */
    {
		double * vex = vex0 + m*NVE ;
		/** transfert coefficient */
		double k_l		= k_i ;
		if (K_dep == 2) k_l = k_i*(1+3*phi_l_i/(1-phi_l_i));
		double d		= d_i ;
		if (D_dep == 2) d = d_i*(1+3*phi_l_i/(1-phi_l_i));
		double tau		= tau_i ;
		if (Tau_dep == 2) tau = tau_i;
		else if (Tau_dep == 3) tau = tau_i*(1-3*0.1*phi_l_i/(2+phi_l_i)) ;
		else if (Tau_dep == 4) tau = tau_i*(1-3*0.5*phi_l_i/(2+phi_l_i)) ;
		else if (Tau_dep == 5) tau = 0.964 ;
		else if (Tau_dep == 6) tau = tau_i*(1-3*0.8*phi_l_i/(2+phi_l_i)) ;
		
		K_l = k_l ;
		D = d ;
		Tau = tau ;
	}
  }
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		EXPLICITES TERMS		*///////////////////////////
////////////////////////////////////////////////////////////////////////
int  ComputeExplicitTerms(Element_t *el,double t)
/* Explicites terms  */
{
  double *vex = Element_GetExplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;   
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;   
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    m ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;
  
  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++ , vex += NVE, vim_n  += NVI) {
	  /* Variables */
    double* x = ComputeVariables(el,u,u,vim_n,0,m) ;
    
    /* porosite */
    double phi_l  = x[I_PHI_L];
    
    /* Volumetric deformation */
    double* eps =  x + I_EPS ;
    double tre   = eps[0] + eps[4] + eps[8] ;

    /* transfert coefficients */
    double k_l, d, tau;

    //		Permeabilite 		// 
    ////////////////////////////////////////////////////////////////////
     
    if (K_dep == 0) k_l = k_i;																						// 0_ Constant
    else if (K_dep == 1) k_l = k_i*(pow(((phi_l/(1+tre))/phi_l_i),3)*(pow(((1-phi_l_i)/(1-(phi_l/(1+tre)))),2))) ;  // 1_ UPC
    else if (K_dep == 2) k_l = k_i *(1+3*phi_l/(1+tre-phi_l)) ;														// 2_ Maxwell
   
    //		Coefficient d efficacite osmotique 		// 
    ////////////////////////////////////////////////////////////////////   

    double n_tau = 1 ;
    double A_tau = 30 ;
    double B_tau = 0.3;
   
    if (Tau_dep == 0) tau = tau_i;										// 0_ Constant
    else if (Tau_dep == 2) tau = tau_i*(phi_l_i*(1+tre)/phi_l) ;		// 2_ UPC
    else if (Tau_dep == 3) tau = tau_i*(1-3*0.1*phi_l/(2+phi_l)) ;		// 3_ Maxwell modifie : a = 0.1
    else if (Tau_dep == 4) tau = tau_i*(1-3*0.5*phi_l/(2+phi_l)) ;		// 4_ Maxwell modifie : a = 0.5 
    else if (Tau_dep == 5) {											// 5_ Constant jusqu a un seuil puis decroissance rapide 
		if ((phi_l/(1+tre)) < 0.55) tau = tau_i ;
		else tau = tau_i*(phi_l_i/(phi_l/(1+tre)-0.55)) ;
		}
	else if (Tau_dep == 6) tau = tau_i*(1-3*0.8*phi_l/(2+phi_l)) ;		// 6_ Maxwell modifie : a = 0.8
    
    else if (Tau_dep == 7) tau = tau_i*(1-1/(1+exp(-A_tau*(phi_l/(1+tre)-B_tau)))) ; 	// 7_ Signoide parametre A_tau (vitesse de decroissance 30 ) B_tau (point de decroissance 0.3)  
    
    //		Coefficient de diffusion		// 
    ////////////////////////////////////////////////////////////////////
    
    if (D_dep == 0) d = d_i;											// 0_ Constant.
    else if (D_dep == 1) d = d_i*(1+3*phi_l/(1+tre-phi_l))*(1-tau) ;	// 1_ Actuel : 		Tau_dep = 0   -    	D_i = 3e-14  m2/s
    else if (D_dep == 2) d = d_i*(1+3*phi_l/(1+tre-phi_l)) ;			// 2_ Novateur : 	Tau_dep = ..  - 	D_i = ....
    else if (D_dep == 3) d = d_i*(1-tau) ; 								// 3_ UPC : 		Tau_dep = 2   - 	D_i = 1.6e-16 m2/s
    else d		= d_i*phi_l/phi_l_i ;
    
    ////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //printf("\n valeur tre %.2f | valeur tau  %.2f | valeur phi_l %.2f ",tre,tau,phi_l);
    //printf("\n valeur w_l  %.2f | valeur w_salt %.2f ",w_l[0]*1e30,w_salt[0]*1e30);
    /////////////////////////////////////////////////////////////////////////////////////
  
    /* storage in vex */
    {
    K_l = k_l ;
    D = d ;
    Tau = tau ;
	}
    /* For post-treatment of M_salt_leached and m_water_absorbed */
	M_Wc_salt += M_W_saltn ;
	M_Wc_water += M_W_watern;
  }

  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		IMPLICITE TERMS	 variables internes*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the implicit terms */
{
  double* vim 	= Element_GetImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    m ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;

  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++ , vim  += NVI) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,dt,m) ;
    
    /* Storage in vim */
    {
		int i ;
		
		M_total	= x[I_M_TOT] ;
		M_salt	= x[I_M_SALT] ;
		
		for(i = 0 ; i < 3 ; i++) W_l[i] = x[I_W_L + i] ;
		for(i = 0 ; i < 3 ; i++) W_salt[i] = x[I_W_SALT + i] ;
		
		Rho_l = x[I_RHO_L] ;		
		Rho_w = x[I_RHO_W] ;		
		// manque Dm_s par rapport à CHM21 ... utile ??
		
		for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
		for(i = 0 ; i < 9 ; i++) SIG_D[i] = x[I_SIG_D + i] ;
		SIG_M = x[I_SIG_M] ;
		
		Phi_l = x[I_PHI_L] ;
		Phi_c = x[I_PHI_C] ;
		
		/* For post-treatment of M_salt_leached and m_water_absorbed */
        M_W_water = (W_l[0] - W_salt[0])*dt ;
        M_W_salt = W_salt[0]*dt ;
	}
  }
  return(0) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		COMPUTE MATRIX		*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
// //  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;													// number of degree of freedom = nn (number of nodes) * NEQ (number of equation)
// //  char *method = Material_GetMethod(Element_GetMaterial(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i,j,dec ;														
  double c[MAX_PGAUSS*121] ;											// c a table containing MAX_PGAUSS*121 double
  double zero = 0. ;

  
  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /***** Poromechanic Matrix *****/
  dec = kp2(fem,c,dt) ;			// Calculation of Cijkl,  121 terms store in a vector c
  {
    double *kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,2) ;
    for(i = 0 ; i < ndof*ndof ; i++) {				
		k[i] = kp[i] ;		// 36 terms
    }
  }
  /***** Conduction Matrix liquid (dW_l/dp) *****/
  dec = klp2(fem,c) ;
  {
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_mass + i*NEQ,I_p + j*NEQ) -= dt*kc[i*nn + j] ;
    }
  }
  /***** Conduction Matrix liquid (dW_l/dw_s) *****/
  dec = klws2(fem,c) ;
  {
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_mass + i*NEQ,I_w_s + j*NEQ) -= dt*kc[i*nn + j] ;
    }
  }
  /***** Conduction Matrix salt (dW_salt/dp) *****/
  dec = ksp2(fem,c) ;
  {
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_salt + i*NEQ,I_p + j*NEQ) -=  dt*kc[i*nn + j] ;
    }
  }
  /***** Conduction Matrix salt (dW_salt/dw_s) *****/
  dec = ksws2(fem,c) ;
  {
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_salt + i*NEQ,I_w_s + j*NEQ) -= dt*kc[i*nn + j] ;
    }
  } 
  return(0) ;
#undef K
}
////////////////////////////////////////////////////////////////////////
////////////*		COMPUTE RESIDU		*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double *vim_n1 = Element_GetPreviousImplicitTerm(el) ;
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


  /***** 1. Mechanics *****/
  /** 1.1 Stresses */
  {
    double *vim = vim_1 ;
    double *rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_mec) -= rw[i] ;		// = residu = -Eq(u)
  }

  /***** 2. Total mass conservation *****/
  /** 2.1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double *vim_n = vim_n1 ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    double *ra ;
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_total - M_totaln ;
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    for(i = 0 ; i < nn ; i++) R(i,E_mass) -= ra[i] ;
  }
  /** 2.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_l,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_mass) -= -dt*rf[i] ; 
  }
  /** 2.3 Elements P2P1 */
  if(strstr(method,"P2P1")) {
    FEM_TransformResiduFromDegree2IntoDegree1(fem,E_mass,r) ;
  }
  
  /***** 3. Salt mass conservation *****/
  /** 3.1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double *vim_n = vim_n1;
    double g1[IntFct_MaxNbOfIntPoints] ;
    double *ra ;
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_salt - M_saltn ;
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    for(i = 0 ; i < nn ; i++) R(i,E_salt) -= ra[i] ;
  }
  /** 3.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_salt,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_salt) -= -dt*rf[i] ; 
  }
  return(0) ;
#undef R	
}
////////////////////////////////////////////////////////////////////////
////////////*		COMPUTE OUTPUTS		*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 17 ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
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
    
    /* unknown */
    double p  =  FEM_ComputeCurrentUnknown(fem,h_s,nn,I_p) ;
    double w_s =  FEM_ComputeCurrentUnknown(fem,h_s,nn,I_w_s) ;
    
    double dis[3] = {0,0,0} ;
    int i,j;
    for(i = 0 ; i < dim ; i++) {										// Displacement calculation
      dis[i] = FEM_ComputeCurrentUnknown(fem,h_s,nn,I_u + i) ;
    }
    
    /* VI */
    double w_l[3] = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double phil = 0 ;
    double phil_Eulerian = 0 ;
    double phic = 0;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double tre = 0 ; 
    
    /* For post-treatment of M_salt_leached and m_water_absorbed */   
    double mwcsalt = 0 ;
    double mwcwater = 0 ;
    double msaltleached = 0 ;
    double mwaterabsorbed = 0 ;
    ////////////////////////////////////////////////////////////////////
    double msalt_element = 0 ;
    double mwater_element = 0 ;
    double indicateur =0 ;
    
    /* VE = k,d and tau */
    double k = 0 ;
    double d = 0 ;
    double tau = 0 ;
    
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {		// np = 2; permet de moyenner sur chaque élément 

      for(j = 0 ; j < dim ; j++) w_l[j] += W_l[j]/np ;
      for(j = 0 ; j < dim ; j++) w_salt[j] += W_salt[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      phil += Phi_l/np ;
      phic += Phi_c/np ;
      mwcsalt += M_Wc_salt/np ;
      mwcwater += M_Wc_water/np ;
          ///////////////////////////////////////////////////////////////////////////////////////////////////
      //~ printf("\n valeur 1 %.2f | valeur deux  %.2f | valeur 3 %i ",Phi_l,phil,np);
     ///////////////////////////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////////////////////////
      //~ printf("\n valeur i %i | valeur np  %i | valeur 3 %i ",i,np,2);
     ///////////////////////////////////////////////////////////////////////////////////////////////////


      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,i,I_u) ;
      tre   += (eps[0] + eps[4] + eps[8])/np ;
      
      //~ phil_Eulerian += (phil / (1+tre))/np ;
      
      k += K_l/np ;
      d += D/np ;
      tau += Tau/np ;
      ////////////////////////////////////////////////
      msalt_element += M_salt* SampleSurface * ElementSize /np ;	// msalt_element [kg] = M_Salt [kg/m3] x SampleSurface [m2] x ElementSize [m] / np [-]
      mwater_element += (M_total - M_salt)* SampleSurface * ElementSize /np ;
    }
        for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {		// np = 2; permet de moyenner sur chaque élément 
      phil_Eulerian += (phil / (1+tre))/np ;
          ///////////////////////////////////////////////////////////////////////////////////////////////////
      //~ printf("\n valeur 1 %.2f | valeur deux  %.2f | valeur 3 %.2f | boucle spéciale pr phil_Eulerian ",Phi_l,phil,phil_Eulerian);
     ///////////////////////////////////////////////////////////////////////////////////////////////////

    }
    msaltleached = -mwcsalt * SampleSurface *1e6 ;
    mwaterabsorbed = mwcwater * SampleSurface *1e6;
      
    i = 0 ;
    Result_Store(r + i++,dis,"					Displacements	u_1			u_2			u_3",3) ;
    Result_Store(r + i++,&p,"Liquide pressure",1) ;
    Result_Store(r + i++,&w_s,"Solute mass fraction",1) ;
    
    Result_Store(r + i++,w_l,"Liquide fluxe",3) ;
    Result_Store(r + i++,w_salt,"Solute fluxe",3) ;
    Result_Store(r + i++,sig,"Stress",9) ;
    Result_Store(r + i++,&phil,"Lagrangian porosity",1) ;
    Result_Store(r + i++,&phic,"Crystal volumetric fraction",1) ;
    Result_Store(r + i++,&msaltleached,"msaltleached",1) ;
    Result_Store(r + i++,&mwaterabsorbed,"mwaterabsorbed",1) ;
    
    Result_Store(r + i++,&k,"k",1) ;
    Result_Store(r + i++,&d,"d",1) ;
    Result_Store(r + i++,&tau,"tau",1) ;
    
    Result_Store(r + i++,&tre,"tre",1) ;
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     //~ printf("\n juste avant la sauvegarde de la valeur  %.2f ",phil_Eulerian);
     ///////////////////////////////////////////////////////////////////////////////////////////////////

    Result_Store(r + i++,&phil_Eulerian,"Eulerian porosity",1) ;   
    /////////////////////////
    Result_Store(r + i++,&msalt_element,"msalt_element",1) ;
    Result_Store(r + i++,&mwater_element,"mwater_element",1) ;
           
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }

  return (NbOfOutputs) ; 
}
////////////////////////////////////////////////////////////////////////
////////////*		MATRIX TERMS		*///////////////////////////////
////////////////////////////////////////////////////////////////////////
int kp2(FEM_t *fem,double *c,double dt)
/* Poro-viscoelastic matrix (c) and shift (dec)*/
// C1 correspond au tenseur meca alors que B1 est la matrice globale 
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  Element_t *el = FEM_GetElement(fem) ;
  double *vim = Element_GetCurrentImplicitTerm(el);
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int    dec = 121 ;
  int    m ;
  double zero = 0. ;

  /** Input data **/
  GetProperties(el) ;
    

  for(m = 0 ; m < np ; m++ , vim += NVI) {
    int i,j ;
    double *c1 = c + m*dec ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,m) ;
    /* unknown */
    double w_s  =  FEM_ComputeCurrentUnknown(fem,h,nn,I_w_s) ;
    
    
    double** u   = Element_ComputePointerToNodalUnknowns(el) ;
    double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,m,I_u) ;
    double tre   = eps[0] + eps[4] + eps[8] ;	// Matrix volumetric strain
    
    double K_s = young/(3.*(1. - 2.*poisson)) ;
    double G_s = young/(2.*(1. + poisson)) ;
    
    
    double  K   = K_s * 4 * (1-(Phi_l/(1+tre))) * G_s / (3*(Phi_l/(1+tre))*K_s + 4*G_s) ;
    double  G   = G_s * (1-(Phi_l/(1+tre))) * (9*K_s + 8*G_s) / (9 * K_s * (1+2*(Phi_l/(1+tre))/3) + 8 * G_s * (1+3*(Phi_l/(1+tre))/2) );
    
    double eta_v_mt = eta_v * 4 * (1-(Phi_l/(1+tre))) * eta_d / (3*(Phi_l/(1+tre))*eta_v + 4*eta_d) ;
    double eta_d_mt = eta_d * (1-(Phi_l/(1+tre))) * (9*eta_v + 8*eta_d) / (9 * eta_v * (1+2*(Phi_l/(1+tre))/3) + 8 * eta_d * (1+3*(Phi_l/(1+tre))/2) );
    
    /** initialisation **/
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;

    /** Mechanics **/		// Maxwell isotropic visco-elastic tensor
    for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
		C1(i,i,j,j) += K / (1. + (K*dt)/eta_v_mt) ;
		C1(i,i,j,j) += - (1./3.)* 2.*G / (1. + (2.*G*dt)/eta_d_mt) ;
		
		C1(i,j,i,j) += 2.*G / (1. + (2.*G*dt)/eta_d_mt) ;

    }
//C1=+81 signifie qu'on part de 81, puis de 81+9, etc ...
// ordres ddl : def gradp gradoms et en vertical : contraintes masses totale et masse de sel  
    c1 += 81 ;
    // Derivative of sigma_ii by p (index 81-89)
    for(i = 0 ; i < 3 ; i++) B1(i,i) = -b ;
    c1 += 9 ;
    // Derivative of sigma_ii by w_s (index 90-98)		(?)
    for(i = 0 ; i < 3 ; i++) B1(i,i) = 0 ;
    c1 += 9 ;
    
    /** Total mass **/	// Derivative of M_total by epsi_ii (index 99-107)
    for(i = 0 ; i < 3 ; i++) B1(i,i) = Rho_l * b ;  
    c1 += 9 ;

    /** Derivative of M_total by p*/ //	(with incompressibility = 0) (index 108)
    if (rho_l_dep == 0) c1[0] = Phi_l * compressibility_w * Rho_w ;
    else 				c1[0] = Phi_l * compressibility_w * Rho_l;  //Phi_l *compressibility_w * Rho_w/(1-w_s);// ((1-w_s)*Rho_l*Rho_l*compressibility_w/Rho_w) ;
    c1 += 1;

    /** Derivative of M_total by w_s*/ //		(= drho_l/dw_s) (index 109)								(= dphi_l/dw_s)												(= dphi_c/dw_s)
    if (rho_l_dep == 0 ) c1[0] = 											   		- (Rho_l * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat)) 	+ rho_c * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat) ;
    else        					c1[0] = - Phi_l  *rho_w_i*(1e-3-v_m_s)*Rho_l  	- (Rho_l * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat)) 	+ rho_c * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat) ; 
    c1 += 1;

    /** Mass of salt **/	// Derivative of M_salt by epsi_ii (index 110-118)
    for(i = 0 ; i < 3 ; i++) B1(i,i) = Rho_l*w_s * b ; 
    c1 += 9 ;
    
    /** Derivative of M_salt by p */ // 	(with incompressibility = 0) (index 119)
    if (rho_l_dep == 0 ) c1[0] = Phi_l * w_s * compressibility_w * Rho_w;
    else c1[0] = 	             Phi_l * w_s * compressibility_w * Rho_l;
    c1 += 1;
    
    /** Derivative of M_salt by w_s	*/ //   (= dphi_c/dw_s) (index 120)										(= d phi_l / dw_s)											(= d rho_l / dw_s)
    if (rho_l_dep == 0 ) c1[0] = + rho_c * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat) + Phi_l*Rho_l - (Rho_l*w_s* ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat));
    else                 c1[0] = + rho_c * ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat) + Phi_l*Rho_l - (Rho_l*w_s* ((1/rho_c)*dt*Phi_c*sigma_c*beta/w_s_sat)) - w_s * Phi_l *rho_w_i*(1e-3-v_m_s)*Rho_l ; 
    ///////////////////////////////////////////////////////////////////////////////////////////////////
     //printf("\n valeur 1 %.2f | valeur deux  %.2f | valeur 3 %.2f ",c1[0],Phi_l,Rho_l);
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     
  }
  
  return(dec) ;

#undef C1
#undef B1
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int klp2(FEM_t *fem,double *c)
/* Conduction matrix (c) and shift (dec) */
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el);
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    m ;
  double zero = 0. ;

  for(m = 0 ; m < np ; m++ , vim_n += NVI , vex += NVE) {
    int i ;
    double *c1 = c + m*dec ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */ 	// (= dW_l/dp)
    c1[0] = - Rho_ln*K_l/viscosity_w ;
    c1[4] = - Rho_ln*K_l/viscosity_w ;
    c1[8] = - Rho_ln*K_l/viscosity_w ;
  }
  
  return(dec) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int klws2(FEM_t *fem,double *c)
/* Conduction matrix (c) and shift (dec) */
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el);
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    m ;
  double zero = 0. ;

  for(m = 0 ; m < np ; m++ , vim_n += NVI , vex += NVE) {
    int i ;
    double *c1 = c + m*dec ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */ 	// (= dW_l/dws)
    c1[0] = Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ;
    c1[4] = Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ;
    c1[8] = Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ; 
  }
  
  return(dec) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int ksp2(FEM_t *fem,double *c)
/* Conduction matrix (c) and shift (dec) */
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el);
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;    
  int    dec = 9 ;
  int    m ;
  double zero = 0. ;

  for(m = 0 ; m < np ; m++ , vim_n += NVI , vex += NVE) {
    int i ;
    double *c1 = c + m*dec ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,m) ;
    /* unknown */
    double w_sn  =  FEM_ComputePreviousUnknown(fem,h,nn,I_w_s) ;    
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */ 	// (= dW_salt/dp)
    c1[0] = - (1-Tau)*Rho_ln*w_sn*K_l/viscosity_w ;
    c1[4] = - (1-Tau)*Rho_ln*w_sn*K_l/viscosity_w ;
    c1[8] = - (1-Tau)*Rho_ln*w_sn*K_l/viscosity_w ;
  }
  
  return(dec) ;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int ksws2(FEM_t *fem,double *c)
/* Conduction matrix (c) and shift (dec) */
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el);
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;      
  int    dec = 9 ;
  int    m ;
  double zero = 0. ;

  for(m = 0 ; m < np ; m++ , vim_n += NVI , vex += NVE) {
    int i ;
    double *c1 = c + m*dec ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,m) ;
    /* unknown */
    double w_sn  =  FEM_ComputePreviousUnknown(fem,h,nn,I_w_s) ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */ 	// (= dW_salt/dws)
    c1[0] = - D*Rho_ln   +   (1-Tau)*w_sn*Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ;
    c1[4] = - D*Rho_ln   +   (1-Tau)*w_sn*Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ;
    c1[8] = - D*Rho_ln   +   (1-Tau)*w_sn*Rho_ln*K_l*RT*Rho_ln*Tau/(viscosity_w*Molarm_s) ;
  }
  
  return(dec) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		COMPUTE VARIABLES		*///////////////////////////
////////////////////////////////////////////////////////////////////////
double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double dt,int m)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
//  double*   x      = Model_GetVariable(model,p) ;
  double*   x      = Variable ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,m,I_u + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,m,I_u) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Pressure */
    x[I_P] = FEM_ComputeUnknown(fem,u,intfct,m,I_p) ;
    
    /* Pressure gradient */
    {
      double* grd_p = FEM_ComputeUnknownGradient(fem,u,intfct,m,I_p) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P + i] = grd_p[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd_p) ;
    }
    
    /* Salt mass fraction */
    x[I_W_S] = FEM_ComputeUnknown(fem,u,intfct,m,I_w_s) ;
    
    /* Salt mass fraction gradient */
    {
      double* grd_w_s = FEM_ComputeUnknownGradient(fem,u,intfct,m,I_w_s) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_W_S + i] = grd_w_s[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd_w_s) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,m,I_u) ;
      double* vim_n = f_n + m*NVI ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS_n   + i] = eps_n[i] ;
        x[I_SIG_n   + i] = SIG_n[i] ;
        x[I_SIG_D_n	+ i] = SIG_D_n[i] ;
      }
      
      x[I_SIG_M_n] = SIG_M_n ;
      /* Crystal salt volumique fraction at previous time step */
      x[I_PHI_C_n] = Phi_c_n ;
      /* liquide density at previous time step */
      x[I_RHO_L_n] = Rho_ln ;
      
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
    
    /* Pressure at previous time step */
    x[I_P_n] = FEM_ComputeUnknown(fem,u_n,intfct,m,I_p) ;
    
    /* Salt mass fraction at previous time step */
    x[I_W_S_n] = FEM_ComputeUnknown(fem,u_n,intfct,m,I_w_s) ;
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + m*NVE ;
      
      x[I_K_L]  = K_l ;
      x[I_D]    = D ;
      x[I_TAU]  = Tau ;
    }
  }
    
  ComputeSecondaryVariables(el,dt,x) ;
  
  return(x) ;
}
////////////////////////////////////////////////////////////////////////
////////////*		COMPUTE SECONDARY VARIABLES	(Flux, eq complémentaires etc.)	*///////////////////
////////////////////////////////////////////////////////////////////////
void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  double *vim = Element_GetCurrentImplicitTerm(el);
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps =  x + I_EPS ;
  double* eps_n =  x + I_EPS_n ;
  /* Pressure */
  double  p   = x[I_P] ;
  //double  p_n = x[I_P_n] ;
  /* Solute mass fraction */
  double w_s	= x[I_W_S] ;
  double w_s_n	= x[I_W_S_n] ;
    
  double  deps[9],trde,deps_d[9] ;
  int    i ;
      
  /* Incremental deformations */
  for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
  trde = deps[0] + deps[4] + deps[8] ;
  
  /* Backup mass flow */
  {
    /* Porosity */
    /** Strain */
    double tre   = eps[0] + eps[4] + eps[8] ;	// Matrix volumetric strain
    double dphi_l_epsi = b*tre ;				// Porosity volumetric strain
    /** crystal salt dissolution */
    double phi_c_n = x[I_PHI_C_n] ;
    double phi_c = phi_c_n   +   (dt * phi_c_n * sigma_c * beta * ((w_s/w_s_sat)-1)) / rho_c ;
    double dphi_l_salt = phi_c_i - phi_c ;
    
    x[I_PHI_C] = phi_c ;
    
    /** Porosity evolution */
    double phi_l = (phi_l_i + dphi_l_epsi + dphi_l_salt);
    //~ if (phi_l < 1e-6) phi_l = 1e-6;			// A priori pas utile
    x[I_PHI_L] = phi_l ;

    /* Fluide density */
    double rho_w,rho_l ;
    rho_w = rho_w_i*(exp(compressibility_w*(p-p_i))) ;
    if (rho_l_dep == 0 ) rho_l = rho_w ;
    else rho_l = rho_w*(exp(-rho_w_i*(1e-3-v_m_s)*w_s));
    
    x[I_RHO_W] = rho_w ;
    x[I_RHO_L] = rho_l ;
    
    /* mass contents */
    double m_total = rho_l*phi_l     + rho_c*phi_c ;
    double m_salt  = rho_l*phi_l*w_s + rho_c*phi_c ;
    x[I_M_TOT] = m_total ;
    x[I_M_SALT] = m_salt ;
    
    /* Transfer coefficient */
    double k_l	= x[I_K_L] ;
    double d	= x[I_D] ;
    double tau	= x[I_TAU] ;
    
    /* Gradients */
    double* grd_p 	= x + I_GRD_P ;
    double* grd_w_s = x + I_GRD_W_S ;
    
    /* Fluxes */
    double* w_f 	= x + I_W_F ;
    double* w_d 	= x + I_W_D ;
    double* w_o 	= x + I_W_O ;
    double* w_l 	= x + I_W_L ;
    double* w_salt 	= x + I_W_SALT ;
    
    double rho_l_n = x[I_RHO_L_n] ;
    
    for(i = 0 ; i < 3 ; i++){
		  w_f[i] = - d*rho_l_n*grd_w_s[i] ;  // modif Rho_ln en rho_l   !!!!!!
		  w_d[i] = - (k_l/viscosity_w)*grd_p[i] ;
		  w_o[i] =   (k_l/viscosity_w)*(RT*rho_l/Molarm_s)*tau*grd_w_s[i] ;
		  
		  w_l[i] = rho_l*(w_d[i] + w_o[i]) ;
		  w_salt[i] = w_f[i] + (1-tau)*w_s_n*  rho_l_n*(w_d[i]+w_o[i]) ;    		// modif Rho_ln en rho_l       !!!!!!!!!                  printf("\n tau=%.2f et 1-tau=%.2f ",tau,1-tau);
		      ///////////////////////////////////////////////////////////////////////////////////////////////////
     //printf("\n valeur w_f %.2f | valeur w_d  %.2f | valeur w_o %.2f ",w_f[0]*1e30,w_d[0]*1e30,w_o[0]*1e30);
     //printf("\n valeur w_l  %.2f | valeur w_salt %.2f ",w_l[0]*1e30,w_salt[0]*1e30);
     ///////////////////////////////////////////////////////////////////////////////////////////////////
	}
	/* Backup stresses, spheric stress and deviatoric stresses */// Mori Tanaka  
  {
	double K_s = young/(3.*(1. - 2.*poisson)) ;
    double G_s = young/(2.*(1. + poisson)) ;
    
    double  K   = K_s * 4 * (1-(phi_l/(1+tre))) * G_s / (3*(phi_l/(1+tre))*K_s + 4*G_s) ;
    double  G   = G_s * (1-(phi_l/(1+tre))) * (9*K_s + 8*G_s) / (9 * K_s * (1+2*(phi_l/(1+tre))/3) + 8 * G_s * (1+3*(phi_l/(1+tre))/2) );
  
    double eta_v_mt = eta_v * 4 * (1-(phi_l/(1+tre))) * eta_d / (3*(phi_l/(1+tre))*eta_v + 4*eta_d) ;
    double eta_d_mt = eta_d * (1-(phi_l/(1+tre))) * (9*eta_v + 8*eta_d) / (9 * eta_v * (1+2*(phi_l/(1+tre))/3) + 8 * eta_d * (1+3*(phi_l/(1+tre))/2) );
    
    //~ if (phi_l < 1e-6) {
		//~ poisson = 0.5;
		//~ eta_v_mt =1e30;
		    //~ double K_s = young/(3.*(1. - 2.*poisson)) ;
		    //~ double G_s = young/(2.*(1. + poisson)) ;
		    //~ double  K   = K_s * 4 * (1-(phi_l/(1+tre))) * G_s / (3*(phi_l/(1+tre))*K_s + 4*G_s) ;
		    //~ double  G   = G_s * (1-(phi_l/(1+tre))) * (9*K_s + 8*G_s) / (9 * K_s * (1+2*(phi_l/(1+tre))/3) + 8 * G_s * (1+3*(phi_l/(1+tre))/2) );
		    //~ double eta_v_mt = eta_v * 4 * (1-(phi_l/(1+tre))) * eta_d / (3*(phi_l/(1+tre))*eta_v + 4*eta_d) ;
		    //~ double eta_d_mt = eta_d * (1-(phi_l/(1+tre))) * (9*eta_v + 8*eta_d) / (9 * eta_v * (1+2*(phi_l/(1+tre))/3) + 8 * eta_d * (1+3*(phi_l/(1+tre))/2) );
		//~ }
    
    double* sig   = x + I_SIG ;
    //double* sig_n = x + I_SIG_n ;
    
    double sig_m = x[I_SIG_M] ;
    double sig_m_n = x[I_SIG_M_n] ;
    
    double* sig_d   = x + I_SIG_D ;
    double* sig_d_n = x + I_SIG_D_n ;
    
    {     
      /* Volumetric stress */
      sig_m = (trde + sig_m_n/K)*K/(1 + (K*dt)/eta_v_mt) ;
      x[I_SIG_M] = sig_m ;
      /* Deviatoric stress */
      /** Deviatoric strains */
	  for(i = 0 ; i < 9 ; i++) deps_d[i] = deps[i] ;
	  deps_d[0] += -1./3. * trde ;
	  deps_d[4] += -1./3. * trde ;
	  deps_d[8] += -1./3. * trde ;
	  /** Deviatoric stresses */
	  for(i = 0 ; i < 9 ; i++) sig_d[i] = (2*G/(1 + (2*G*dt)/eta_d_mt))*(deps_d[i] + sig_d_n[i]/(2*G)) ;
	  
	  /* Total stresses */
	  for(i = 0 ; i < 9 ; i++) sig[i] = sig_d[i] ;
      sig[0] += sig_m + sig0_11 - b*(p - p_i) ;
      sig[4] += sig_m + sig0_22 - b*(p - p_i) ;
      sig[8] += sig_m + sig0_33 - b*(p - p_i) ;
      
    }
  }
  }
}
