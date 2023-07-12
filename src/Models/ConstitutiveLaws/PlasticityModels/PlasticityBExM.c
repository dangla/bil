static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityBExM_CT ;
static Plasticity_ReturnMapping_t                    PlasticityBExM_RM ;
static Plasticity_YieldFunction_t                    PlasticityBExM_YF ;
static Plasticity_FlowRules_t                        PlasticityBExM_FR ;
static Plasticity_SetParameters_t                    PlasticityBExM_SP ;
static Plasticity_SetModelProp_t                     PlasticityBExM_SetModelProp ;


#define Plasticity_GetSlopeSwellingLine(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetSlopeVirginConsolidationLine(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetSlopeCriticalStateLine(PL) \
        Plasticity_GetParameter(PL)[2]
        
#define Plasticity_GetInitialPreconsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[3]
        
#define Plasticity_GetMacroscopicInitialVoidRatio(PL) \
        Plasticity_GetParameter(PL)[4]
        
#define Plasticity_GetMicroscopicInitialVoidRatio(PL) \
        Plasticity_GetParameter(PL)[5]

#define Plasticity_GetSuctionCohesionCoefficient(PL) \
        Plasticity_GetParameter(PL)[6]

#define Plasticity_GetReferenceConsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[7]
        
#define Plasticity_GetAlphaMicrostructure(PL) \
        Plasticity_GetParameter(PL)[8]
        
#define Plasticity_GetInitialSuctionDecreaseHardening(PL) \
        Plasticity_GetParameter(PL)[9]
        
#define Plasticity_GetInitialSuctionIncreaseHardening(PL) \
        Plasticity_GetParameter(PL)[10]
        
#define Plasticity_GetMicroElasticStiffnessEffectivePressureChanges(PL) \
        Plasticity_GetParameter(PL)[11]
        
#define Plasticity_GetLoadingCollapseFactorCurve(PL) \
        Curves_GetCurve(Plasticity_GetCurves(PL))
        
#define Plasticity_GetSuctionIncreasingFactor(PL) \
        (Curves_GetCurve(Plasticity_GetCurves(PL)) + 1)
        
#define Plasticity_GetSuctionDecreasingFactor(PL) \
        (Curves_GetCurve(Plasticity_GetCurves(PL)) + 2)
        
#define Plasticity_GetSaturationDegreeCurve(PL) \
        (Curves_GetCurve(Plasticity_GetCurves(PL)) + 3)
        
        



void PlasticityBExM_SetModelProp(Plasticity_t* plasty)
{
  
  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityBExM_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityBExM_RM ;
    Plasticity_GetYieldFunction(plasty)                 = PlasticityBExM_YF ;
    Plasticity_GetFlowRules(plasty)                     = PlasticityBExM_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityBExM_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 4 ;
    Plasticity_GetNbOfCriteria(plasty)                  = 2 ;
    
  }
  
}


void PlasticityBExM_SP(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetMacroscopicInitialVoidRatio(plasty)     = va_arg(args,double) ;
    Plasticity_GetMicroscopicInitialVoidRatio(plasty)     = va_arg(args,double) ;
    Plasticity_GetSuctionCohesionCoefficient(plasty)      = va_arg(args,double) ;
    Plasticity_GetReferenceConsolidationPressure(plasty)  = va_arg(args,double) ;
    Plasticity_GetAlphaMicrostructure(plasty)             = va_arg(args,double) ;
    Plasticity_GetInitialSuctionDecreaseHardening(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialSuctionIncreaseHardening(plasty) = va_arg(args,double) ;
    Plasticity_GetMicroElasticStiffnessEffectivePressureChanges(plasty) = va_arg(args,double) ;
    Curve_t* lc                                           = va_arg(args,Curve_t*) ;
    Curve_t* fi                                           = va_arg(args,Curve_t*) ;
    Curve_t* fd                                           = va_arg(args,Curve_t*) ;
    Curve_t* sl                                           = va_arg(args,Curve_t*) ;
    int i = Curves_Append(Plasticity_GetCurves(plasty),lc) ;
    
    if(i != 0) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    i = Curves_Append(Plasticity_GetCurves(plasty),fi) ;
    
    if(i != 1) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    i = Curves_Append(Plasticity_GetCurves(plasty),fd) ;
    
    if(i != 2) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    i = Curves_Append(Plasticity_GetCurves(plasty),sl) ;
    
    if(i != 3) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      double si = Plasticity_GetInitialSuctionIncreaseHardening(plasty) ;
      double sd = Plasticity_GetInitialSuctionDecreaseHardening(plasty) ;

      Plasticity_GetHardeningVariable(plasty)[0] = log(pc) ;
      Plasticity_GetHardeningVariable(plasty)[1] = si ;
      Plasticity_GetHardeningVariable(plasty)[2] = sd ;
      
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[1] = 1.e-6*si ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[2] = 1.e-6*sd ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[3] = 1 ;
    }
  }

  va_end(args) ;
}



double* (PlasticityBExM_CT)(Plasticity_t* plasty,const double* sig,const double* hardv,const double* zeta)
/** Barcelona Expansive model.
 **/
{
  double* crit = Plasticity_GenericTangentStiffnessTensor(plasty,sig,hardv,zeta) ;
  
  return(crit) ;
}



double* (PlasticityBExM_RM)(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Barcelona Expansive model.
 **/
{
  double* crit = Plasticity_GenericReturnMapping(plasty,sig,eps_p,hardv) ;
  
  return(crit) ;
}




double* (PlasticityBExM_YF)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the value of the yield functions. 
 **/
{
  size_t SizeNeeded = 2*sizeof(double) ;
  double* yield   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m     = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double k     = Plasticity_GetSuctionCohesionCoefficient(plasty) ;
  double p_r   = Plasticity_GetReferenceConsolidationPressure(plasty) ;
  double alpha = Plasticity_GetAlphaMicrostructure(plasty) ;
  Curve_t* lc  = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
  Curve_t* sl  = Plasticity_GetSaturationDegreeCurve(plasty) ;
  double lnpc0 = hardv[0] ;
  double si    = hardv[1] ;
  double sd    = hardv[2] ;
  double s     = hardv[3] ;
  double ps    = k*s ;
  double lc_s  = Curve_ComputeValue(lc,s) ;
  double sl_s  = Curve_ComputeValue(sl,s) ;
  double lnp_r = log(p_r) ;
  double lnpc  = lnp_r + lc_s * (lnpc0 - lnp_r) ;
  double pc    = exp(lnpc) ;
  double m2    = m*m ;
  double p     = (stress[0] + stress[4] + stress[8])/3. ;
  double q     = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
  double peff  = p - pow(sl_s,alpha)*s ;
  
  yield[0] = q*q/m2 + (p - ps)*(p + pc) ;
  yield[1] = fabs(peff + 0.5*(si + sd)) - 0.5*(si + sd) ;

  return(yield) ;
}



double* (PlasticityBExM_FR)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Barcelona Expansive model.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the slope of the critical state line (M),
 *  the suction cohesion coefficient (k),
 *  the reference consolidation pressure (p_r),
 *  the loading collapse factor curve (lc),
 *  the pre-consolidation pressure at suction=0 (pc0=exp(hardv[0])),
 *  the initial void ratio (e0).
 *  the suction (s=hardv[1]),
 * 
 *  Return the direction of the plastic flows based on the flow rules:
 *    - the plastic strain rate (i.e. the potential gradient)
 *    - the rate of log(pre-consolidation pressure) (1/dlambda * d(ln(pc0))/dt)
 **/
{
  size_t SizeNeeded = 2*(9+4)*(sizeof(double)) ;
  double* flow   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double eM0     = Plasticity_GetMacroscopicInitialVoidRatio(plasty) ;
  double em0     = Plasticity_GetMicroscopicInitialVoidRatio(plasty) ;
  double k       = Plasticity_GetSuctionCohesionCoefficient(plasty) ;
  double p_r     = Plasticity_GetReferenceConsolidationPressure(plasty) ;
  double alpha   = Plasticity_GetAlphaMicrostructure(plasty) ;
  double kappa_m = Plasticity_GetMicroElasticStiffnessEffectivePressureChanges(plasty) ;
  Curve_t* lc    = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
  Curve_t* sl    = Plasticity_GetSaturationDegreeCurve(plasty) ;
  Curve_t* fi    = Plasticity_GetSuctionIncreasingFactor(plasty) ;
  Curve_t* fd    = Plasticity_GetSuctionDecreasingFactor(plasty) ;
  double lnpc0   = hardv[0] ;
  double si      = hardv[1] ;
  double sd      = hardv[2] ;
  double s       = hardv[3] ;
  double ps      = k*s ;
  double lc_s    = Curve_ComputeValue(lc,s) ;
  double sl_s    = Curve_ComputeValue(sl,s) ;
  double lnp_r   = log(p_r) ;
  double lnpc    = lnp_r + lc_s * (lnpc0 - lnp_r) ;
  double pc      = exp(lnpc) ;
  double m2      = m*m ;
  double p       = (stress[0] + stress[4] + stress[8])/3. ;
  double q       = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress)) ;
  double peff    = p - pow(sl_s,alpha)*s ;
  double pstar   = - p - q*q/(m2*(p - ps)) ;
  double fi_p    = Curve_ComputeValue(fi,pstar/pc) ;
  double fd_p    = Curve_ComputeValue(fd,pstar/pc) ;
  double beta    = 1 ;
  double id[9]   = {1,0,0,0,1,0,0,0,1} ;
  
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
    The hardening flow: d(ln(pc0)) = - (1 + eM0)*v*deps_p
   * --------------------------------------------------
   * Using a = ln(pc0) as hardening variable.
   * d(a) = - dl * (1 + eM0) * v * (dg/dp)
   * So h(p,a) = - (1 + eM0) * v * (2*p + pc(a) - ps)
   */
  {
    double v = 1./(lambda - kappa) ;

    flow[9]  = - (1 + eM0)*v*(2*p + pc - ps) ;
    flow[10] = 0 ;
    flow[11] = 0 ;
    flow[12] = 0 ;
  }
  
  {
    double* flow2 = flow + 9 + 4 ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      flow2[i] = (peff + 0.5*(si + sd))*id[i]/3 ;
    }
  }
  
  {
    double* flow2 = flow + 9 + 4 ;
    double v = 1./(lambda - kappa) ;
    double hm = (peff + 0.5*(si + sd)) ;
    double km = (1 + em0)/kappa_m *(-peff) ;
    
    if(km < 0) {
      arret("") ;
    }

    flow2[9]  = - (1 + eM0)*v*hm ;
    flow2[10] = - km/fi_p * hm ;
    flow2[11] = - km/fd_p * hm ;
    flow2[12] = 0 ;
  }
  
  return(flow) ;
}


#undef Plasticity_GetSlopeSwellingLine
#undef Plasticity_GetSlopeVirginConsolidationLine
#undef Plasticity_GetSlopeCriticalStateLine
#undef Plasticity_GetInitialPreconsolidationPressure
#undef Plasticity_GetMacroscopicInitialVoidRatio
#undef Plasticity_GetMicroscopicInitialVoidRatio
#undef Plasticity_GetSuctionCohesionCoefficient
#undef Plasticity_GetReferenceConsolidationPressure
#undef Plasticity_GetAlphaMicrostructure
#undef Plasticity_GetInitialSuctionDecreaseHardening
#undef Plasticity_GetInitialSuctionIncreaseHardening
#undef Plasticity_GetMicroElasticStiffnessEffectivePressureChanges
#undef Plasticity_GetLoadingCollapseFactorCurve
#undef Plasticity_GetSaturationDegreeCurve
#undef Plasticity_GetSuctionIncreasingFactor
#undef Plasticity_GetSuctionDecreasingFactor
