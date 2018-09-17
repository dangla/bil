#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Message.h"
#include "DataFile.h"
#include "Tools/Math.h"
#include "Plasticity.h"

typedef double    (Plasticity_Criterion_t)(Plasticity_t*,const double*,const double) ;


/* Drucker-Prager */
static int    pmDruckerPrager(const char* s) ;
static double Plasticity_CriterionDruckerPrager(Plasticity_t*,const double*,const double,const double,const double,const double,const double,const double) ;
static double Plasticity_ReturnMappingDruckerPrager(Plasticity_t*,double*,double*,double*,const double,const double,const double,const double,const double) ;


/* Cam-Clay */
static int    pmCamClay(const char* s) ;
static double Plasticity_CriterionCamClay(Plasticity_t*,const double*,const double,const double,const double,const double,const double) ;
static double Plasticity_ReturnMappingCamClay(Plasticity_t*,double*,double*,double*,const double,const double,const double,const double,const double) ;




Plasticity_t*  (Plasticity_Create)(void)
{
  Plasticity_t* plasty = (Plasticity_t*) malloc(sizeof(Plasticity_t)) ;
  
  assert(plasty) ;
  
  /* Allocation of space for the code name of the model */
  {
    size_t sz = Plasticity_MaxLengthOfKeyWord*sizeof(char) ;
    char* name = (char*) malloc(sz) ;
    
    assert(name) ;
    
    Plasticity_GetCodeNameOfModel(plasty) = name ;
  }
  
  /* Allocation of space for the yield function gradient */
  {
    size_t sz = 9*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetYieldFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the potential function gradient */
  {
    size_t sz = 9*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetPotentialFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the hardening modulus */
  {
    size_t sz = sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetHardeningModulus(plasty) = c ;
  }
  
  /* Allocation of space for Fji*Cijkl */
  {
    size_t sz = 9*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetFjiCijkl(plasty) = c ;
  }
  
  /* Allocation of space for Cijkl*Glk */
  {
    size_t sz = 9*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetCijklGlk(plasty) = c ;
  }
  
  /* Allocation of space for the tangent stiffness tensor */
  {
    size_t sz = 81*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Plasticity_GetTangentStiffnessTensor(plasty) = c ;
  }
  
  /* Elasticity */
  Plasticity_GetElasticity(plasty) = Elasticity_Create() ;
  
  return(plasty) ;
}



void  (Plasticity_Delete)(void* self)
{
  Plasticity_t** pplasty = (Plasticity_t**) self ;
  Plasticity_t*  plasty  = *pplasty ;
  
  {
    char* name = Plasticity_GetCodeNameOfModel(plasty) ;
    free(name) ;
  }
  
  {
    double* c = Plasticity_GetYieldFunctionGradient(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetPotentialFunctionGradient(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetHardeningModulus(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetFjiCijkl(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetCijklGlk(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    free(c) ;
  }
  
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    Elasticity_Delete(&elasty) ;
  }
  
  free(*pplasty) ;
  *pplasty = NULL ;
}







int pmDruckerPrager(const char* s)
{
       if(!strcmp(s,"cohesion"))     return (0) ;
  else if(!strcmp(s,"friction"))     return (1) ;
  else if(!strcmp(s,"dilatancy"))    return (2) ;
  else if(!strcmp(s,"alpha"))        return (4) ;
  else if(!strcmp(s,"gamma_R"))      return (5) ;
  else return(-1) ;
}







int pmCamClay(const char* s)
{
       if(!strcmp(s,"kappa"))        return (0) ;
  else if(!strcmp(s,"lambda"))       return (1) ;
  else if(!strcmp(s,"dilatancy"))    return (2) ;
  else if(!strcmp(s,"M"))            return (4) ;
  else if(!strcmp(s,"phi"))          return (5) ;
  else return(-1) ;
}






#if 0
void (Plasticity_ScanProperties)(Plasticity_t* plasty,DataFile_t* datafile,int (*pm)(const char*))
/** Read the plastic properties in the stream file ficd */
{
  FILE* ficd = DataFile_GetFileStream(datafile) ;
  int    nd = Plasticity_GetNbOfProperties(plasty) ;
  short int    cont = 1 ;
  
  if(!ficd) return ;

  while(cont) {
    char   mot[Plasticity_MaxLengthOfKeyWord] ;
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    
    if(!line) break ;
    
    sscanf(line," %[^= ]",mot) ;

    /* Reading some curves */
    if(!strncmp(mot,"Courbes",6) || !strncmp(mot,"Curves",5)) {
      Curves_t* curves = Plasticity_GetCurves(plasty) ;
      Curves_FreeBuffer(curves) ;
      Curves_ReadCurves(curves,line) ;
      
      if(Curves_GetNbOfCurves(curves) > Plasticity_MaxNbOfCurves) {
        arret("Plasticity_ScanProperties (2) : trop de courbes") ;
      }
      
    /* Reading the material properties and storing through pm */
    } else if(pm) {
      char   *p = strchr(line,'=') ;

      /* We assume that this is a property as long as "=" is found */
      if(p) {
        int i = (*pm)(mot) ;
        
        if(i >= 0) {
        
          sscanf(p+1,"%lf",Plasticity_GetProperty(plasty) + i) ;
          nd = (nd > i + 1) ? nd : i + 1 ;
        
        } else {
        
          Message_RuntimeError("%s is not known",mot) ;
          
        }
      
      /* ... otherwise we stop reading */
      } else {
        
        /* go out */
        cont = 0 ;
      }
      
    } else {
      break ;
    }
    
  }

  Plasticity_GetNbOfProperties(plasty) = nd ;
  return ;
}
#endif







double Plasticity_UpdateElastoplasticTensor(Plasticity_t* plasty,double* c)
/** Update the 4th rank elastoplastic tensor in c.
 *  Inputs are: 
 *  dfsds is the yield function gradient
 *  dgsds is the potential function gradient
 *  hm    is the hardening modulus
 *  Outputs are:
 *  fc(k,l) = dfsds(j,i) * C^el(i,j,k,l)
 *  cg(i,j) = C^el(i,j,k,l) * dgsds(l,k)
 *  fcg = hm + dfsds(j,i) * C^el(i,j,k,l) * dgsds(l,k))
 *  Tensor c is then updated as 
 *  C(i,j,k,l) = C^el(i,j,k,l) + cg(i,j) * fc(k,l) / fcg
 *  Return the inverse of fcg: 1/fcg */
{
#define CEP(i,j)  ((c)[(i)*9+(j)])
#define CEL(i,j)  ((cel)[(i)*9+(j)])
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double* cel    = Elasticity_GetStiffnessTensor(elasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double  fcg    = hm[0] ;
  //double* fc     = Plasticity_GetFjiCijkl(plasty) ;
  //double* cg     = Plasticity_GetCijklGlk(plasty) ;
  
  /* Tangent elastoplastic tensor */
  {
      double  fc[9] ;
      double  cg[9] ;
      int i ;
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        fc[i] = 0. ;
        cg[i] = 0. ;
          
        for(j = 0 ; j < 9 ; j++) {
              
          fc[i] += dfsds[j]*CEL(j,i) ;
          cg[i] += CEL(i,j)*dgsds[j] ;
          fcg   += dfsds[i]*CEL(i,j)*dgsds[j] ;
        }
      }
        
      if(fcg > 0.) {
        fcg = 1./fcg ;
      } else {
            
        printf("\n") ;
        printf("dfsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dfsds[i]) ;
        }
        printf("\n") ;
        printf("dgsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dgsds[i]) ;
        }
        printf("\n") ;
        printf("fcg = %e\n",fcg) ;
        printf("\n") ;
        
        arret("Plasticity_UpdateElastoplasticTensor") ;
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          CEP(i,j) = CEL(i,j) - cg[i]*fc[j]*fcg ;
        }
      }
  }
  
  return(fcg) ;
#undef CEP
}



double Plasticity_CriterionDruckerPrager(Plasticity_t* plasty,const double* sig,const double gam_p,const double young,const double poisson,const double af,const double ad,const double cohesion)
/** Drucker-Prager criterion. 
 *  Inputs are: 
 *  the stresses (sig), the cumulative plastic shear strain (gam_p). 
 *  Parameters are:
 *  the friction angle (af), the dilatancy angle (ad) and the cohesion.
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 *  Return the value of the yield function. */
{
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9] ;
  double crit ;
  double cc ;
  
  
  /*
    Cohesion
  */
  {
    //double c1 = (gam_p < gam_R) ? 1 - (1 - alpha)*gam_p/gam_R : alpha ;
    double c1 = 1 ;
    
    cc = cc0*c1*c1 ;
  }
  
  /*
    Yield function
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q + ff*p - cc ;
  
  /*
    Deviatoric stresses
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dev[i] = sig[i] - p*id[i] ;
    }
  }
  
  /*
    Yield function gradient
  */
  if(q > 0.) {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = 1.5*dev[i]/q + id[i]*ff/3. ;
    }
    
  } else {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = id[i]*ff/3. ;
    }
  }
  
  /*
    Potential function gradient
  */
  
  /* Elastic case */
  if(crit <= 0.) {
    if(q > 0.) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    } else {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = id[i]*dd/3. ;
      }
    }
  }
  
  /* Plastic case */
  if(crit > 0.) {
    double k   = young/(1. - 2.*poisson)/3. ;
    double dmu = young/(1.+poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    /* Flow regime at the notch apex */
    } else {
      double dl = (ff*p - cc)/(k*ff*dd) ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = dev[i]/(dmu*dl) + id[i]*dd/3. ;
      }
    }
  }
  
  /* Hardening modulus */
  *hm = 0. ;
  /*
  if(gam_p < gam_R) {
    *hm = -2.*(1.-alpha)/gam_R*(1.-(1.-alpha)*gam_p/gam_R)*cc0 ;
    *hm *= sqrt(2*j2(dgsds)) ;
  }
  */
  
  Plasticity_GetCriterionValue(plasty) = crit ;
  return(crit) ;
}



double Plasticity_ReturnMappingDruckerPrager(Plasticity_t* plasty,double* sig,double* eps_p,double* gam_p,const double young,const double poisson,const double af,const double ad,const double cohesion)
/** Drucker-Prager return mapping.
 *  Parameters are:
 *  the Young modulus (young),
 *  the Poisson's ratio (poisson),
 *  the friction angle (af), 
 *  the dilatancy angle (ad),
 *  the cohesion (cohesion).
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the cumulative plastic shear strain (gam_p).
 *  Return the value of the yield function. */
{
  double dmu     = young/(1.+poisson) ;
  double mu      = dmu/2. ;
  double k       = young/(1. - 2.*poisson)/3. ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,sdev_t[9] ;
  double crit ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*j2(sig)) ;
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      sdev_t[i] = sig[i] - p_t*id[i] ;
    }
  }
  
  /*
    Criterion
  */ 
  {
    //double c1 = ((*gam_p) < gam_R) ? 1 - (1 - alpha)*(*gam_p)/gam_R : alpha ;
    double c1 = 1 ;
    double cc = cc0*c1*c1 ;
    
    crit = q_t + ff*p_t - cc ;
  }
  
  /*
    Return mapping: update plastic strains and stresses
  */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double p = p_t ;
    double q = q_t ;
    double gam_pn = *gam_p ;
    double gam_p1 = *gam_p ;
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      double fcrit = crit ;
      int    nf    = 0 ;
      double dl    = 0 ;
      double tol   = 1.e-10 ;
      
      //while(fabs(fcrit) > tol*cc0) {
      while(fabs(fcrit/q) > tol) {
        double dqsdl ;
        double dpsdl ;
        double dccsdl ;
        double c1 ;
        double cc ;
        int    i ;
        
	      /* Plastic strain increments */
	      for(i = 0 ; i < 9 ; i++) {
	        deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
	      }
        
	      /* Cumulative plastic shear strain */
	      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
        
	      /* p, q, cc */
	      //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
	      c1 = 1 ;
	      cc = cc0*c1*c1 ;
	      q  = q_t - dl*3*mu ;
	      p  = p_t - dl*k*dd ;
        
	      /* dqsdl, dpsdl, dccsdl */
	      dqsdl = -3*mu ;
	      dpsdl = -k*dd ;
	      dccsdl = 0. ;
        /*
	      if(gam_p1 < gam_R) {
	        dccsdl = -2*(1 - alpha)/gam_R*c1*cc0 ;
	        dccsdl *= 1.5*sqrt(2*j2(sdev_t))/q_t ;
	      }
        */
        
        /* Criterion */
        fcrit = q + ff*p - cc ;
        
	      /* dl */
        {
          double df = dqsdl + ff*dpsdl - dccsdl ;
          
	        dl   -= fcrit/df ;
        }
        
	      if(nf++ > 20) {
	        printf("No convergence (Plasticity_ReturnMappingDruckerPrager)") ;
	        exit(0) ;
	      }
      }
      
      /* Stresses */
      {
        int    i ;
        
        for(i = 0 ; i < 9 ; i++) {
	        sig[i] = sdev_t[i]*q/q_t + p*id[i] ;
        }
      }
    }
    
    /* Flow regime at the notch apex */
    if(q <= 0.) {
      double c1 ;
      double cc ;
      double dl ;
      int    i ;
      
      /* Deviatoric plastic strain increments */
      for(i = 0 ; i < 9 ; i++) deps_p[i]  = sdev_t[i]/dmu ;
        
      /* Cumulative plastic shear strain */
      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
      
      /* p, q, cc */
      //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
      c1 = 1 ;
      cc = cc0*c1*c1 ;
      p  = cc/ff ;
      q  = 0. ;
      
      /* dl */
      dl   = (ff*p_t - cc)/(k*ff*dd) ;
      
      /* Plastic strain increments and stresses */
      for(i = 0 ; i < 9 ; i++) {
	      deps_p[i] = sdev_t[i]/dmu + dl*id[i]*dd/3. ;
	      sig[i]    = p*id[i] ;
      }
    }
    
      
    /* Total plastic strains */
    {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        eps_p[i] += deps_p[i] ;
      }
    }
    (*gam_p) = gam_p1 ;
  }
  
  return(crit) ;
}





double Plasticity_CriterionCamClay(Plasticity_t* plasty,const double* sig,const double pc,const double m,const double kappa,const double lambda,const double phi0)
/** Cam-Clay criterion */
{
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double m2      = m*m ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Gradients
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
    
      dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
      dgsds[i] = dfsds[i] ;
    }
  }
  
  /* The hardening modulus */
  {
    double v = 1./(lambda - kappa) ;
    
    *hm = v/(1 - phi0)*p*(2*p + pc)*pc ;
  }
  return(crit) ;
}



double Plasticity_ReturnMappingCamClay(Plasticity_t* plasty,double* sig,double* eps_p,double* p_co,const double m,const double kappa,const double lambda,const double phi0,const double mu)
/** Cam-Clay return mapping. Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (p_co),
 *  the porosity or the void ratio (phi0,e0).
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the hardening variable (p_co).
 *  Return the value of the yield function. 
 *  Algorithm from Borja & Lee 1990 modified by Dangla. */
{
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,p,q,pc,crit ;
  double dl ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  pc   = *p_co ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + p*(p + pc) = 0
     for p. The other variables (pc,q,dl) are expressed with p.
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.e-8 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double dlsdp  = ((1 - phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*dlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
	      printf("no convergence (ReturnMapping_CamClay)") ;
	      exit(0) ;
      }
    }
  }
  
  /*
    Plastic stresses and strains
  */
  
  {
    double a = 1./(1 + 6*mu/m2*dl) ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev      = a*(sig[i] - p_t*id[i]) ;
      double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
      sig[i]   = p*id[i] + dev ;
      eps_p[i] = dl*dfsds ;
    }
  }
  
  /* Consolidation pressure */
  *p_co = pc ;
  return(crit) ;
}



