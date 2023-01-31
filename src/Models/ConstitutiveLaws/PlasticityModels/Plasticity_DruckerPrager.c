static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_CTDruckerPrager ;

static Plasticity_ReturnMapping_t                    Plasticity_RMDruckerPrager ;


double Plasticity_CTDruckerPrager(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
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
  //double gam_p   = hardv[0] ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9],devn[9] ;
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
    double k   = young/(3 - 6*poisson) ;
    double dmu = young/(1 + poisson) ;
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
  hm[0] = 0. ;
  /*
  if(gam_p < gam_R) {
    hm[0] = -2.*(1.-alpha)/gam_R*(1.-(1.-alpha)*gam_p/gam_R)*cc0 ;
    hm[0] *= sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(dgsds)) ;
  }
  */
  
  Plasticity_GetCriterionValue(plasty) = crit ;
  
  /*
   * Consistent tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
     
    if(dlambda > 0 && q > 0.) {
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
        double D = q - 3*g1*dlambda ;
        double bnn = 4*g1*g1*dlambda/D ;
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
    }
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
   
  return(crit) ;
}



double Plasticity_RMDruckerPrager(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
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
  double p_t,q_t,sdev_t[9] ;
  double crit ;
  double dl = 0 ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
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
    //double c1 = ((gam_p) < gam_R) ? 1 - (1 - alpha)*(gam_p)/gam_R : alpha ;
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
    double gam_pn = gam_p ;
    double gam_p1 = gam_p ;
    
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
        
        /* Plastic strain increments */
        for(i = 0 ; i < 9 ; i++) {
          deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
        }
        
        /* Cumulative plastic shear strain */
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
          dccsdl = -2*(1 - alpha)/gam_R*c1*cc0 ;
          dccsdl *= 1.5*sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(sdev_t))/q_t ;
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
          Message_FatalError("Plasticity_RMDruckerPrager: no convergence") ;
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
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}
