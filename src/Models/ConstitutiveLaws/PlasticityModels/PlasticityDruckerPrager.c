#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include "Message.h"
#include "Math_.h"
#include "Plasticity.h"
#include "autodiff.h"

static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityDruckerPrager_CT ;
static Plasticity_ReturnMapping_t                    PlasticityDruckerPrager_RM ;
static Plasticity_YieldFunction_t                    PlasticityDruckerPrager_YF ;
static Plasticity_FlowRules_t                        PlasticityDruckerPrager_FR ;
static Plasticity_SetParameters_t                    PlasticityDruckerPrager_SP ;
extern Plasticity_SetModelProp_t                     PlasticityDruckerPrager_SetModelProp ;


        
#define Plasticity_GetFrictionAngle(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetDilatancyAngle(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetCohesion(PL) \
        Plasticity_GetParameter(PL)[2]
        
#define Plasticity_GetCohesionFactorCurve(PL) \
        (Curves_GetNbOfCurves(Plasticity_GetCurves(PL)) ? \
        Curves_GetCurve(Plasticity_GetCurves(PL)) : NULL)
        
        



void PlasticityDruckerPrager_SetModelProp(Plasticity_t* plasty)
{

  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityDruckerPrager_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityDruckerPrager_RM ;
    Plasticity_GetYieldFunction(plasty)                 = PlasticityDruckerPrager_YF ;
    Plasticity_GetFlowRules(plasty)                     = PlasticityDruckerPrager_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityDruckerPrager_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 0 ;
  }
  
}



void PlasticityDruckerPrager_SP(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  {
    Plasticity_GetFrictionAngle(plasty)            = va_arg(args,double) ;
    Plasticity_GetDilatancyAngle(plasty)           = va_arg(args,double) ;
    Plasticity_GetCohesion(plasty)                 = va_arg(args,double) ;
    Curve_t* cofac                                 = va_arg(args,Curve_t*) ;
    int i = Curves_Append(Plasticity_GetCurves(plasty),cofac) ;
    
    if(cofac && i != 0) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    {
      Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
      double k = Elasticity_GetBulkModulus(elasty) ;
      
      //Plasticity_GetHardeningVariable(plasty)[0] = 0 ;
      
      //Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 0 ;
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-8*k ;
    }
    
  }

  va_end(args) ;
}



double* (PlasticityDruckerPrager_CT)(Plasticity_t* plasty,const double* sig,const double* hardv,const double* plambda)
/** Drucker-Prager criterion. 
 * 
 *  Inputs are: 
 *  the stresses (sig),
 *  the cumulative plastic shear strain (gam_p = hardv[0], not used here). 
 * 
 *  Parameters are:
 *  the friction angle (af),
 *  the dilatancy angle (ad) and the cohesion.
 * 
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  double* yield  = Plasticity_GetCriterionValue(plasty) ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  Curve_t* cofac = Plasticity_GetCohesionFactorCurve(plasty) ;
  double gam_p   = hardv[0] ;
  double fac     = (cofac) ? Curve_ComputeValue(cofac,gam_p) : 1. ;
  double cc      = cc0*fac ;
  
  double dlambda = (plambda) ? plambda[0] : 0 ;
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9],devn[9] ;
  double crit ;

  
  /*
    Yield function
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
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
  
  /* The derivative of q wrt sig_ij: devn_ij */
  {
    int    i ;
    
    if(q > 0.) {
      for(i = 0 ; i < 9 ; i++) {
        devn[i] = 1.5 * dev[i] / q ;
      }
    } else {
      for(i = 0 ; i < 9 ; i++) {
        devn[i] = 0 ;
      }
    }
  }
  
  /*
    Yield function gradient
  */
  {
    int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        dfsds[i] = devn[i] + id[i]*ff/3. ;
      }
  }
  
  /*
    Potential function gradient
  */
  {
    int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = devn[i] + id[i]*dd/3. ;
      }
  }
  
  /* Plastic case */
  #if 0
  if(crit > 0.) {
    double k   = young/(3 - 6*poisson) ;
    double dmu = young/(1 + poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = devn[i] + id[i]*dd/3. ;
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
  #endif
  
  /* Hardening modulus */
  {
    if(q > 0) {
      double dfac = (cofac) ? Curve_ComputeDerivative(cofac,gam_p) : 0 ;
      double dcc  = cc0*dfac ;
      
      hm[0] = dcc * sqrt(1.5) ; /* choice 1. */
    } else {
      hm[0] = 0 ;
    }
  }
  
  /*
   * Consistent tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
     
    if(q > 0 && dlambda > 0) {
      double bulk    = Elasticity_GetBulkModulus(elasty) ;
      double shear   = Elasticity_GetShearModulus(elasty) ;
      double g0 = shear ;
      double k0 = bulk  ;
      double g1 = g0 * q / (3*g0*dlambda + q) ;
      double k1 = k0 ;
   
      {
        double mu      = g1 ;
        double lame    = k1 - 2./3.*mu ;
        int    i ;

        for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
        for(i = 0 ; i < 3 ; i++) {
          int j ;
      
          #define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
          for(j = 0 ; j < 3 ; j++) {
            C(i,i,j,j) += lame ;
            C(i,j,i,j) += mu ;
            C(i,j,j,i) += mu ;
          }
          #undef C
        }
      }
       
      {
        //double D = q - 3*g1*dlambda ;
        double D = q*q / (3*g0*dlambda + q) ;
        //double bnn = 4*g1*g1*dlambda/D ;
        double bnn = 4*g0*g0*dlambda / (3*g0*dlambda + q) ;
        int    i ;
        
        for(i = 0 ; i < 9 ; i++) {
          int j ;
         
          #define C1(i,j)    (c[(i)*9+(j)])
          for(j = 0 ; j < 9 ; j++) {
            C1(i,j) += bnn*devn[i]*devn[j] ;
          }
          #undef C1
        }
      }
       
      Plasticity_UpdateElastoplasticTensor(plasty,c) ;
    } else {
      double bulk    = Elasticity_GetBulkModulus(elasty) ;
      double h       = hm[0] ;
      /* 1/K1 = 1/K + ff*dd/H -> K1 = K/(1+ff*dd*K/H) */
      double k1      = bulk ;
   
      {
        double mu      = 0 ;
        double lame    = bulk - 2./3.*mu ;
        int    i ;

        for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
        for(i = 0 ; i < 3 ; i++) {
          int j ;
      
          #define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
          for(j = 0 ; j < 3 ; j++) {
            C(i,i,j,j) += lame ;
            C(i,j,i,j) += mu ;
            C(i,j,j,i) += mu ;
          }
          #undef C
        }
      }
    }
  }
  
  yield[0] = crit ;
  
  return(yield) ;
}


#if 0
double* (PlasticityDruckerPrager_RM)(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Drucker-Prager return mapping.
 * 
 *  Parameters are:
 *  the Young modulus (young),
 *  the Poisson's ratio (poisson),
 *  the friction angle (af), 
 *  the dilatancy angle (ad),
 *  the cohesion (cohesion).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the cumulative plastic shear strain (gam_p = hardv[0]).
 *  the plastic multiplier in plasty.
 * 
 *  Return the value of the yield function. */
{
  double* yield  = Plasticity_GetCriterionValue(plasty) ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = dmu/2. ;
  double k       = young/(3 - 6*poisson) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  double gam_p   = hardv[0] ;

  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t ;
  double q_t ;
  double crit ;
  double dl = 0 ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  
  /*
    Criterion
  */ 
  {
    //double c1 = ((gam_p) < gam_R) ? 1 - (1 - alpha)*(gam_p)/gam_R : alpha ;
    double c1 = 1 ;
    double cc = cc0*c1*c1 ;
    
    crit = q_t + ff*p_t - cc ;
  }
  
  /*
   * Directions of the plastic strain rates
   * --------------------------------------
   *
   * They are given by the potential function gradient of g(sig) = q + dd*p
   * The potential function gradient dg_ij is the subdifferential:
   * 
   * dg_ij = {y_ij:  y_kl (sig'_kl - sig_kl) <= G(sig') - G(sig) for any sig'_kl}
   * 
   * i.e.
   * 
   * dg_ij = 1/3 dd Id_ij + n_ij with n_ij a deviator defined as:
   * 
   * if q > 0: n_ij = 3/2 s_ij/q 
   * if q = 0: n_ij = {y_ij:   y_kl s'_lk <= q' for any deviator s'_ij}
   */

  /*
   * Return mapping: update plastic strains and stresses
   */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double p = p_t ;
    double q = q_t ;
    double gam_pn = gam_p ;
    double gam_p1 = gam_p ;
    double sdev_t[9] ;
    
    /* Deviatoric trial stresses */
    {
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        sdev_t[i] = sig[i] - p_t*id[i] ;
      }
    }
    
    /* General laws:
     * ------------
     * 
     * deps_v^p = dl * dd
     * de_ij^p  = dl * n_ij
     * 
     * p        = p_t    - dl * k * dd
     * s_ij     = s_ij_t - dl * 2mu * n_ij
     * 
     * gam_p is the cumulative plastic shear strain rate:
     * 
     * gam_p    = gam_pn + dl * sqrt(n_kl n_lk)
     */
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      double fcrit = crit ;
      int    nf    = 0 ;
      double tol   = 1.e-10 ;
      
      //while(fabs(fcrit) > tol*cc0) {
      while(fabs(fcrit/q) > tol) {
        double dqsdl ;
        double dpsdl ;
        double dccsdl ;
        double c1 ;
        double cc ;
        int    i ;
        
        /*
         * If q > 0 then
         * -------------
         *     n_ij  = 3/2 s_ij / q
         * 
         *     p     = p_t - k * dd * dl
         *     s_ij (1 + dl * 3mu / q) = s_ij_t
         *     q     = q_t - dl * 3mu
         * 
         *     We note that
         *     n_ij  = 3/2 s_ij_t / q_t
         *     sqrt(n_kl n_lk) = sqrt(3/2)
         * 
         *     and so
         *     gam_p = gam_pn + dl * sqrt(3/2)
         * 
         *     cc    = cc(gam_p)
         * 
         *     This is valid provided that
         *     q > 0 i.e. q_t > 3mu * dl
         * 
         *     The consistency equation, F = q + ff*p - cc = 0, is solved
         *     for the plastic multiplier, dl:
         * 
         *     dl += - F/dF
         * 
         *     with dF = dq/ddl + ff*dp/ddl - cc'*dgam_p/ddl i.e.
         * 
         *     dF = - 3*mu - ff*k*dd - cc'*sqrt(3/2)
         * 
         */
        
        /* Plastic strain increments */
        for(i = 0 ; i < 9 ; i++) {
          deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
        }
        
        /* Cumulative plastic shear strain:
         * 
         * gamma_p = gamma_pn + sqrt(3/2) * dl
         * cc(gamma_p) = cc(gamma_pn) + cc'(gamma_pn) * sqrt(3/2) * dl
         * h = cc'(gamma_pn) * sqrt(3/2)
         */
        gam_p1 = gam_pn + sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(deps_p)) ;
        
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
        * cc'    = - cc0*2*c1*(1 - alpha)/gam_R
          h      = cc' * sqrt(3/2) ;
          dccsdl = h ;
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
          Message_FatalError("PlasticityDruckerPrager_RM: no convergence") ;
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
      int    i ;
      
      /* Deviatoric plastic strain increments */
      for(i = 0 ; i < 9 ; i++) deps_p[i]  = sdev_t[i]/dmu ;
        
      /* Cumulative plastic shear strain */
      gam_p1 = gam_pn + sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(deps_p)) ;
      
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
    hardv[0] = gam_p1 ;
  }
  
  /* Plastic multiplier */
  Plasticity_GetPlasticMultiplier(plasty)[0] = dl ;
  
  yield[0] = crit ;
  
  return(yield) ;
}
#endif



#if 1
double* (PlasticityDruckerPrager_RM)(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Drucker-Prager return mapping.
 * 
 *  Parameters are:
 *  the Young modulus (young),
 *  the Poisson's ratio (poisson),
 *  the friction angle (af), 
 *  the dilatancy angle (ad),
 *  the cohesion (cohesion).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the cumulative plastic shear strain (gam_p = hardv[0]).
 *  the plastic multiplier in plasty.
 * 
 *  Return the value of the yield function. */
{
  double* yield  = Plasticity_GetCriterionValue(plasty) ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = dmu/2. ;
  double k       = young/(3 - 6*poisson) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  Curve_t* cofac = Plasticity_GetCohesionFactorCurve(plasty) ;
  double gam_p   = hardv[0] ;
  double fac     = (cofac) ? Curve_ComputeValue(cofac,gam_p) : 1. ;
  double cc      = cc0*fac ;

  double p       = (sig[0] + sig[4] + sig[8])/3. ;
  double q       = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  double crit    = q + ff*p - cc ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double dl = 0 ;
  
  /*
   * Directions of the plastic strain rates
   * --------------------------------------
   *
   * They are given by the potential function gradient of g(sig) = q + dd*p
   * The potential function gradient dg_ij is the subdifferential:
   * 
   * dg_ij = {y_ij:  y_kl (sig'_kl - sig_kl) <= G(sig') - G(sig) for any sig'_kl}
   * 
   * i.e.
   * 
   * dg_ij = 1/3 dd Id_ij + n_ij with n_ij a deviator defined as:
   * 
   * if q > 0: n_ij = 3/2 s_ij/q 
   * if q = 0: n_ij = {y_ij:   y_kl s'_lk <= q' for any deviator s'_ij}
   */
  
  /*
   * Return mapping: update plastic strains and stresses
   */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    /* Trial stresses */
    double sdev_t[9] ;
    double p_t = p ;
    double q_t = q ;
    double gam_pn = gam_p ;
  
    /* Deviatoric trial stresses */
    {
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        sdev_t[i] = sig[i] - p*id[i] ;
      }
    }
    
    /* General laws:
     * ------------
     * 
     * deps_v^p = dl * dd
     * de_ij^p  = dl * n_ij
     * 
     * p        = p_t    - dl * k * dd
     * s_ij     = s_ij_t - dl * 2mu * n_ij
     * 
     * The hardening variable is:
     *   1. either the cumulative plastic shear strain rate:
     *      gam_p = gam_pn + dl * sqrt(n_kl n_lk)
     * 
     *   2. or the volumetric plastic strain:
     *      gam_p = gam_pn + dl * dd
     * 
     */
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      int iter = 0 ;
      int niter = 20 ;
      int convergencenotattained = 1 ;
      
      dl = 0 ;
      
      while(convergencenotattained) {
        /*
         * If q > 0 then
         * -------------
         *     n_ij  = 3/2 s_ij / q
         * 
         *     p     = p_t - k * dd * dl
         *     s_ij (1 + dl * 3mu / q) = s_ij_t
         *     q     = q_t - dl * 3mu
         * 
         *     We note that
         *     n_ij  = 3/2 s_ij_t / q_t
         *     n_kl n_lk = 3/2
         * 
         *     and so
         *     1. gam_p = gam_pn + dl * sqrt(3/2)
         *     2. gam_p = gam_pn + dl * dd
         * 
         *     cc    = cc(gam_p)
         * 
         *     This is valid provided that
         *     q > 0 i.e. q_t > 3mu * dl
         * 
         *     The consistency equation, F = q + ff*p - cc = 0, is solved
         *     for the plastic multiplier, dl:
         * 
         *     dl += - F/dF
         * 
         *     with dF = dq/ddl + ff*dp/ddl - cc'*dgam_p/ddl i.e.
         * 
         *     dF = - 3*mu - ff*k*dd - cc'*sqrt(3/2)
         * 
         */
        
        /* Plastic strain increments */
        {
          int i ;
          
          for(i = 0 ; i < 9 ; i++) {
            deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
          }
        }
        
        /* Hardening variable: */
        gam_p = gam_pn + dl * sqrt(1.5) ; // choice 1.
        //gam_p = gam_pn + dl * dd ; // choice 2.
        
        /* Cohesion */
        fac   = (cofac) ? Curve_ComputeValue(cofac,gam_p) : 1. ;
        cc    = cc0*fac ;
        
        /* Stresses */
        q     = q_t - dl*3*mu ;
        p     = p_t - dl*k*dd ;
        
        /* Update dl */
        {
          double fcrit = q + ff*p - cc ;
          double dfac = (cofac) ? Curve_ComputeDerivative(cofac,gam_p) : 0 ;
          double dcc  = cc0*dfac ;
          double dqsdl = - 3*mu;
          double dpsdl = - k*dd;
          double dccsdl = dcc * sqrt(1.5) ; // choice 1.
          //double dccsdl = dcc * dd ; // choice 2.
          double df = dqsdl + ff*dpsdl - dccsdl ;
          double ddl = (df != 0) ? -fcrit/df : 0 ;
          
          dl += ddl ;
      
          /* Convergence check */
          {
            double tol = 1.e-6 ;
        
            convergencenotattained = 0 ;
          
            if(fabs(ddl) > tol*fabs(dl)) {
              convergencenotattained = 1 ;
            }
          }
        }
        
        if(iter++ > niter) {
          Message_FatalError("PlasticityDruckerPrager_RM: no convergence") ;
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
      
    /* Flow regime at the notch apex 
     * (at the output of the previous calculation q can be <= 0)  */
    if(q <= 0) {
      int iter = 0 ;
      int niter = 20 ;
      int convergencenotattained = 1 ;
      
      dl = 0 ;
      
      while(convergencenotattained) {
        /*
         * If q = 0 (s_ij = 0) then
         * ------------------------
         *     n_ij s'_ji <= q' for any deviator s'_ij
         * 
         *     p    = p_t - dl * k * dd
         *     s_ij = 0
         *     q    = 0
         * 
         *     s_ij_t = dl * 2*mu * n_ij
         * 
         *     The product n_ij s'_ji is maximum for s'_ji = s_ij_t
         * 
         *     so that
         *     n_ij s_ij_t <= q_t
         *     dl          >= q_t/3G
         * 
         *     de_ij^p = s_ij_t / (2*mu)
         *     q_t     = dl * 3mu * sqrt(2/3) * sqrt(n_kl n_lk)
         * 
         *     1. gam_p = gam_pn + sqrt(3/2) * q_t / (3*mu)
         *     2. gam_p = gam_pn + dl * dd
         * 
         *     cc    = cc(gam_p)
         * 
         *     The consistency equation, F = ff*p - cc = 0, is solved
         *     for the plastic multiplier, dl:
         * 
         *     dl += - F/dF
         * 
         *     with dF = ff*dp/ddl - cc'*dgam_p/ddl i.e.
         * 
         *     dF = - ff*k*dd
         */
        
        /* Hardening variable: */
        gam_p = gam_pn + sqrt(1.5) * q_t / (3*mu) ; // choice 1.
        //gam_p = gam_pn + dl * dd ; // choice 2.
        
        /* Cohesion */
        fac   = (cofac) ? Curve_ComputeValue(cofac,gam_p) : 1. ;
        cc    = cc0*fac ;
        
        /* Plastic strain increments */
        {
          int i ;
          
          for(i = 0 ; i < 9 ; i++) {
            deps_p[i] = sdev_t[i]/(2*mu) + dl * id[i]*dd/3. ;
          }
        }
      
        /* Stresses */
        q  = 0 ;
        p  = p_t - dl*k*dd ;
        
        /* Update dl */
        {
          double fcrit = ff*p - cc ;
          double dfac = (cofac) ? Curve_ComputeDerivative(cofac,gam_p) : 0 ;
          double dcc  = cc0*dfac ;
          double dpsdl = - k*dd;
          double dccsdl = 0 ; /* choice 1. */
          //double dccsdl = dcc * dd ; // choice 2.
          double df = ff*dpsdl - dccsdl ;
          double ddl = (df != 0) ? -fcrit/df : 0 ;
          
          dl += ddl ;
      
          /* Convergence check */
          {
            double tol = 1.e-6 ;
        
            convergencenotattained = 0 ;
          
            if(fabs(ddl) > tol*fabs(dl)) {
              convergencenotattained = 1 ;
            }
          }
        }
        
        if(iter++ > niter) {
          Message_FatalError("PlasticityDruckerPrager_RM: no convergence") ;
        }
      }
      
      /* Stresses */
      {
        int    i ;
        
        for(i = 0 ; i < 9 ; i++) {
          sig[i] = p*id[i] ;
        }
      }
    }
  
    if(dl < 0) {
      Message_FatalError("PlasticityDruckerPrager_RM: negative plastic multiplier") ;
    }
      
    /* Total plastic strains */
    {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        eps_p[i] += deps_p[i] ;
      }
    }
    
    /* Hardening variable */
    hardv[0] = gam_p ;
  }
  
  /* Plastic multiplier */
  Plasticity_GetPlasticMultiplier(plasty)[0] = dl ;
  
  yield[0] = crit ;
  
  return(yield) ;
}
#endif




double* (PlasticityDruckerPrager_YF)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the value of the yield function. */
{
  size_t SizeNeeded = sizeof(double) ;
  double* yield   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  //double gam_p   = hardv[0] ;
  double p       = (stress[0] + stress[4] + stress[8])/3. ;
  double q       = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
  
  {
    //double c1 = (gam_p < gam_R) ? 1 - (1 - alpha)*gam_p/gam_R : alpha ;
    double c1 = 1 ;
    double cc = cc0*c1*c1 ;
  
    yield[0] = q + ff*p - cc ;
  }

  return(yield) ;
}




double* (PlasticityDruckerPrager_FR)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Drucker-Prager criterion. 
 * 
 *  Inputs are: 
 *  the stresses (stress),
 *  the cumulative plastic shear strain (gam_p = hardv[0], not used here). 
 * 
 *  Parameters are:
 *  the friction angle (af),
 *  the dilatancy angle (ad) and the cohesion.
 * 
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  size_t SizeNeeded = (9+1)*(sizeof(double)) ;
  double* flow   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  //double gam_p   = hardv[0] ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double dev[9] ;
  double p    = (stress[0] + stress[4] + stress[8])/3. ;
  double q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
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
  crit = q + ff*p - cc ;
  
  /*
    Deviatoric stresses
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dev[i] = stress[i] - p*id[i] ;
    }
  }

  
  /*
   * Potential function: G(sig) = q + dd*p
   * Potential function gradient dg_ij is the subdifferential:
   * 
   * dg_ij = {y_ij:  y_kl (sig'_kl - sig_kl) <= G(sig') - G(sig) for any sig'_kl}
   * 
   * i.e.
   * 
   * dg_ij = 1/3 dd Id_ij + n_ij with n_ij a deviator defined as:
   * 
   * n_ij = 3/2 s_ij/q   if q > 0
   * n_ij = {y_ij:   y_kl s'_lk <= q' for any deviator s'_ij}  if q = 0
   */
  
  /* Elastic case */
  if(crit <= 0.) {
    if(q > 0.) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        flow[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    } else {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        flow[i] = id[i]*dd/3. ;
      }
    }
  }
  
  /* Plastic */
  if(crit > 0.) {
    double k   = young/(3 - 6*poisson) ;
    double dmu = young/(1 + poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        flow[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    /* Flow regime at the notch apex */
    } else {
      double dl = (ff*p - cc)/(k*ff*dd) ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        flow[i] = dev[i]/(dmu*dl) + id[i]*dd/3. ;
      }
    }
  }
  
   
  return(flow) ;
}


        
#undef Plasticity_GetFrictionAngle
#undef Plasticity_GetDilatancyAngle
#undef Plasticity_GetCohesion
