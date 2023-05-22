static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityCamClayOffset_CT ;
static Plasticity_ReturnMapping_t                    PlasticityCamClayOffset_RM ;
static Plasticity_YieldFunction_t                    PlasticityCamClayOffset_YF ;
static Plasticity_FlowRules_t                        PlasticityCamClayOffset_FR ;
static Plasticity_SetParameters_t                    PlasticityCamClayOffset_SP ;
static Plasticity_SetModelProp_t                     PlasticityCamClayOffset_SetModelProp ;


#define Plasticity_GetSlopeSwellingLine(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetSlopeVirginConsolidationLine(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetSlopeCriticalStateLine(PL) \
        Plasticity_GetParameter(PL)[2]
        
#define Plasticity_GetInitialPreconsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[3]
        
#define Plasticity_GetInitialVoidRatio(PL) \
        Plasticity_GetParameter(PL)[4]
        
        



void PlasticityCamClayOffset_SetModelProp(Plasticity_t* plasty)
{
  
  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityCamClayOffset_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityCamClayOffset_RM ;
    Plasticity_GetYieldFunction(plasty)                 = PlasticityCamClayOffset_YF ;
    Plasticity_GetFlowRules(plasty)                     = PlasticityCamClayOffset_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityCamClayOffset_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 2 ;
  }
  
}


void PlasticityCamClayOffset_SP(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialVoidRatio(plasty)                = va_arg(args,double) ;
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      
      //Plasticity_GetHardeningVariable(plasty)[0] = pc ;
      Plasticity_GetHardeningVariable(plasty)[0] = log(pc) ;
      
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6 ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[1] = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
    }
    
  }

  va_end(args) ;
}



double (PlasticityCamClayOffset_CT)(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
/** Modified Cam-Clay criterion with offset 
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 *  the isotropic tensile strength (ps=hardv[1]),
 * 
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 *  c     = tangent stiffness tensor
 * 
 *  Return the value of the yield function. 
 **/
{
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double pc      = exp(hardv[0]) ;
  double ps      = hardv[1] ;
  double m2      = m*m ;
  double beta    = 1 ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + (p - ps)*(p + pc) ;
  
  /*
    Gradients
    ---------
    dp/dsig_ij = 1/3 delta_ij
    dq/dsig_ij = 3/2 dev_ij/q 
    df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q 
    df/dp      = 2*p + pc - ps
    df/dq      = 2*q/m2
    
    df/dsig_ij = 1/3 (2*p + pc - ps) delta_ij + (3/m2) dev_ij 
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
    
      dfsds[i] = (2*p + pc - ps)*id[i]/3 + 3/m2*dev ;
      dgsds[i] = (2*p + pc - ps)*id[i]/3 + 3*beta/m2*dev ;
    }
  }
  
  
  /* The hardening modulus */
  /* df = (df/dsig_ij)*dsig_ij + (df/da)*da
   * given da = dl*h then df = (df/dsig_ij)*dsig_ij - dl*H
   * with H the hardening modulus: H = - (df/da) * h
   * On the other hand with a = ln(pc)
   * d(a) = - (1 + e0)*v*deps_p = - (1 + e0)*v*dl*(dg/dp)
   * i.e. h = - (1 + e0)*v*(dg/dp)
   */
  {
    double v = 1./(lambda - kappa) ;
    double v1 = (1 + e0)*v ;
    double h = - v1*(2*p + pc - ps) ;
    double dpcda    = pc ;
    double dfda     = (p - ps)*dpcda ;
    
    hm[0] = - dfda * h ;
  }
  
  /*
   * Tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;

    /* 
     * Using a = ln(pc) as hardening variable instead of pc.
     * So h(p,a) = - v1*(2*p + pc(a) - ps)
     */
    if(dlambda > 0) {
      double v        = 1./(lambda - kappa) ;
      double v1       = (1 + e0)*v ;
      double h        = - v1*(2*p + pc - ps) ;
      double dpcda    = pc ;
      double dhda     = - v1*dpcda ;
      double dhdp     = - v1*2 ;
      double dfda     = (p - ps)*dpcda ;
      double ddgdpda  = dpcda ;
      double dlambda1 = dlambda / (1 - dlambda*dhda) ;
      Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
      double bulk    = Elasticity_GetBulkModulus(elasty) ;
      double shear   = Elasticity_GetShearModulus(elasty) ;
      double g0 = shear ;
      double k0 = bulk  ;
      double g1 = g0 / (6*g0*dlambda*beta/m2 + 1) ;
      double k1 = k0 / (2*k0*dlambda1 + 1) ;
   
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
        int    i ;
    
        for(i = 0 ; i < 9 ; i++) {
          dfsds[i] += dlambda1*dfda*dhdp*id[i]/3 ;
          dgsds[i] += dlambda1*h*ddgdpda*id[i]/3 ;
        }
      }
      
      hm[0] /= (1 - dlambda*dhda) ;
    }
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  
  return(crit) ;
}



double (PlasticityCamClayOffset_RM)(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Modified Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Dangla.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 *  The tensile strength offsetting (ps=hardv[1]).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the pre-consolidation pressure (pc=hardv[0]).
 * 
 *  Return the value of the yield function. 
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double bulk    = Elasticity_GetBulkModulus(elasty) ;
  double mu      = Elasticity_GetShearModulus(elasty) ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double pc      = exp(hardv[0]) ;
  double ps      = hardv[1] ;
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  double beta    = 1 ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  double p_t,q_t ;
  double dl ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + (p - ps)*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + (p - ps)*(p + pc) = 0
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
    double v1     = (1+e0)*v ;
    double klub   = 1/bulk ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      /*
       * Flow rule
       * ---------
       * g(p,q,a)   = q*q*beta/m2 + (p - ps)*(p + pc(a)) = 0
       * dg/dp      = 2*p + pc - ps
       * dg/dq      = (2*beta/m2) q
       * dp/dsig_ij = 1/3 delta_ij
       * dq/dsig_ij = 3/2 dev_ij/q
       * dg/dsig_ij = 1/3 (2*p + pc - ps) delta_ij + (3*beta/m2) dev_ij 
       * deps_p     = dl (2*p + pc - ps)
       * deij_p     = dl (3*beta/m2) dev_ij
       */
      double dfsdp  = 2*p + pc - ps ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p - ps ;
      double dpcsdp = v1*klub*pc ;
      double d2fsdp2  = 2 + dpcsdp ;
      double ddlsdp = (-klub - dl*d2fsdp2)/dfsdp ;
      double dqsdp  = -6*mu*beta*ddlsdp*q/(m2 + 6*mu*beta*dl) ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      /* Plastic multiplier (dl):
       * ------------------------
       * deps_e = (p - p_n) / bulk ;
       * deps   = (p_t - p_n) / bulk ;
       * deps_p = (p_t - p) / bulk = dl (2*p + pc - ps)
       * Hence 
       * dl = (p_t - p) / (bulk*(2*p + pc - ps))
       * Pre-consolidation pressure (pc):
       * --------------------------------
       * deps_p = - (lambda - kappa)/(1+e0) * ln(pc/pc_n)
       * Hence using the above relation of deps_p
       * ln(pc/pc_n) = (p - p_t) * (1+e0) / (bulk*(lambda-kappa))
       * Deviatoric behavior (q):
       * ------------------------
       * dev_ij = dev_ij_t - 2 mu deij_p 
       *        = dev_ij_t - 6 mu beta / m2 dl dev_ij
       * Hence 
       * dev_ij = dev_ij_t / (1 + 6 mu beta / m2 dl)
       * q      = q_t / (1 + 6 mu beta / m2 dl)
       */

      pc     = pc_n*exp((p-p_t)*v1*klub) ;
      dl     = (p_t - p)*klub/(2*p + pc - ps) ;
      q      = q_t*m2/(m2 + 6*mu*beta*dl) ;
      fcrit  = q*q/m2 + (p - ps)*(p + pc) ;
      
      if(nf++ > 20) {
        Message_FatalError("PlasticityCamClayOffset_RM: no convergence") ;
      }
    }
  }
  
  /*
    Stresses and plastic strains
  */
  
  {
    double a = m2/(m2 + 6*mu*beta*dl) ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev      = a*(sig[i] - p_t*id[i]) ;
      double dgsds    = (2*p + pc - ps)*id[i]/3 + 3*beta/m2*dev ;
    
      sig[i]    = p*id[i] + dev ;
      eps_p[i] += dl*dgsds ;
    }
  }
  
  /* Consolidation pressure */
  hardv[0] = log(pc) ;
  
  /* Plastic muliplier */
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}




double (PlasticityCamClayOffset_YF)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the value of the yield function. 
 **/
{
  double m     = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double pc    = exp(hardv[0]) ;
  double ps    = hardv[1] ;
  double m2    = m*m ;
  double p     = (stress[0] + stress[4] + stress[8])/3. ;
  double q     = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
  double yield = q*q/m2 + (p - ps)*(p + pc) ;

  return(yield) ;
}



double* (PlasticityCamClayOffset_FR)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Modified Cam-Clay criterion with offset.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 *  the isotropic tensile strength (ps=hardv[1]),
 * 
 *  Return the direction of the plastic flows based on the flow rules:
 *    - the plastic strain rate (i.e. the potential gradient)
 *    - the rate of log(pre-consolidation pressure) (1/dlambda * d(ln(pc))/dt)
 **/
{
  size_t SizeNeeded = (9+2)*(sizeof(double)) ;
  double* flow   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double pc      = exp(hardv[0]) ;
  double ps      = hardv[1] ;
  double m2      = m*m ;
  double beta    = 1 ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p     = (stress[0] + stress[4] + stress[8])/3. ;
  //double q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
  
  /*
    Potential function: g = beta*q*q/m2 + (p - ps)*(p + pc)
  */
  
  /*
    Plastic strain: deps^p_ij = dl * dg/dstress_ij
    ----------------------------------------------
    dp/dstress_ij = 1/3 delta_ij
    dq/dstress_ij = 3/2 dev_ij/q 
    dg/dstress_ij = 1/3 (dg/dp) delta_ij + 3/2 (dg/dq) dev_ij/q 
    dg/dp      = 2*p + pc - ps
    dg/dq      = beta*2*q/m2
    
    dg/dstress_ij = 1/3 (2*p + pc - ps) delta_ij + beta*(3/m2) dev_ij 
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = stress[i] - p*id[i] ;

      flow[i] = (2*p + pc - ps)*id[i]/3 + 3*beta/m2*dev ;
    }
  }
  
  /*
    The hardening flow: d(ln(pc)) = - (1 + e0)*v*deps_p
   * --------------------------------------------------
   * Using a = ln(pc) as hardening variable.
   * d(a) = - dl * (1 + e0) * v * (dg/dp)
   * So h(p,a) = - (1 + e0) * v * (dg/dp)
   */
  {
    double v = 1./(lambda - kappa) ;
    
    //flow[9] = - (1 + e0)*v*(2*p + pc - ps)*pc ;
    flow[9]  = - (1 + e0)*v*(2*p + pc - ps) ;
    flow[10] = 0 ;
  }
  
  return(flow) ;
}
