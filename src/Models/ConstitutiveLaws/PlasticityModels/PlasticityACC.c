static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityACC_CT ;
static Plasticity_ReturnMapping_t                    PlasticityACC_RM ;
static Plasticity_YieldFunction_t                    PlasticityACC_YF ;
static Plasticity_FlowRules_t                        PlasticityACC_FR ;
static Plasticity_SetParameters_t                    PlasticityACC_SP ;
static Plasticity_SetModelProp_t                     PlasticityACC_SetModelProp ;


static void PlasticityACC_ScaleStress(Plasticity_t*,double*) ;


#define Elasticity_GetInverseStiffnessTensor \
        Elasticity_GetComplianceTensor


#define Plasticity_GetACC_k(PL) \
        Plasticity_GetParameter(PL)[0]

#define Plasticity_GetACC_M(PL) \
        Plasticity_GetParameter(PL)[1]
        
#define Plasticity_GetACC_N(PL) \
        Plasticity_GetParameter(PL)[2]
        
#define Plasticity_GetInitialIsotropicTensileLimit(PL) \
        Plasticity_GetParameter(PL)[3]
        
#define Plasticity_GetACCInitialPreconsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[4]

#define Plasticity_GetVolumetricStrainHardeningParameter(PL) \
        Plasticity_GetParameter(PL)[5]

#define Plasticity_GetThermalHardeningParameter(PL) \
        Plasticity_GetParameter(PL)[6]

#define Plasticity_GetACC_f_par(PL) \
        Plasticity_GetParameter(PL)[7]

#define Plasticity_GetACC_f_perp(PL) \
        Plasticity_GetParameter(PL)[8]
        
        



void PlasticityACC_SetModelProp(Plasticity_t* plasty)
{
  
  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityACC_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityACC_RM ;
    Plasticity_GetYieldFunction(plasty)                 = PlasticityACC_YF ;
    Plasticity_GetFlowRules(plasty)                     = PlasticityACC_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityACC_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 1 ;
  }
  
}
        

void PlasticityACC_SP(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  {
    Plasticity_GetACC_k(plasty) = va_arg(args,double) ;
    Plasticity_GetACC_M(plasty) = va_arg(args,double) ;
    Plasticity_GetACC_N(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialIsotropicTensileLimit(plasty)       = va_arg(args,double) ;
    Plasticity_GetACCInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetVolumetricStrainHardeningParameter(plasty) = va_arg(args,double) ;
    Plasticity_GetThermalHardeningParameter(plasty)          = va_arg(args,double) ;
    Plasticity_GetACC_f_par(plasty)  = va_arg(args,double) ;
    Plasticity_GetACC_f_perp(plasty) = va_arg(args,double) ;
    
    {
      double pc = Plasticity_GetACCInitialPreconsolidationPressure(plasty) ;
      
      //Plasticity_GetHardeningVariable(plasty)[0] = pc ;
      Plasticity_GetHardeningVariable(plasty)[0] = log(pc) ;
      
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6*log(pc) ;
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
    }
    
  }

  va_end(args) ;
}



double* PlasticityACC_CT(Plasticity_t* plasty,const double* stress,const double* hardv,const double* plambda)
/** Asymmetric Cam-Clay criterion 
 *  first hardening parameter is the preconsolidation pressure
 *  second hardening parameter is the isotropic tensile elastic limit
*/
{
  double* crit   = Plasticity_GetCriterionValue(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;

  {
    double* yield = PlasticityACC_YF(plasty,stress,hardv) ;
  
    crit[0] = yield[0] ;
    
    Plasticity_FreeBufferFrom(plasty,yield) ;
  }

  /*
    Flow directions
  */
  {
    int    i ;
    double* dyield = Plasticity_DerivativeOfYieldFunction(plasty,PlasticityACC_YF,stress,hardv) ;
    double* flow = PlasticityACC_FR(plasty,stress,hardv) ;
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = dyield[i] ;
      dgsds[i] = flow[i] ;
    }
  
    /* The hardening modulus */
    /* H is defined by: df = (df/dsig_ij) dsig_ij - dl H
     * H = - df/da * ha
     */

    hm[0] = - dyield[9] * flow[9] ;
    
    Plasticity_FreeBufferFrom(plasty,dyield) ;
  }
  
  /*
   * Continuum tangent stiffness matrix
   */
   #if 0
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  #endif
    
  /*
   * Consistent tangent stiffness matrix
   */
  {
    Plasticity_GenericTangentStiffnessTensor(plasty,stress,hardv,plambda) ;
  }
  
  return(crit) ;
}



double* PlasticityACC_RM(Plasticity_t* plasty,double* stress,double* strain_p,double* hardv)
{
  double* crit = Plasticity_GenericReturnMapping(plasty,stress,strain_p,hardv) ;
  
  return(crit) ;
}



double* PlasticityACC_YF(Plasticity_t* plasty,const double* stress,const double* hardv)
{
  size_t SizeNeeded = sizeof(double) ;
  double* yield   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m  = Plasticity_GetACC_M(plasty)  ;
  double k  = Plasticity_GetACC_k(plasty)   ;
  double ps = Plasticity_GetInitialIsotropicTensileLimit(plasty) ;
  double pc = exp(hardv[0]) ;

  double stress_tilde[9] ;
  
  {
    int i ;
    
    for (i = 0 ; i < 9; i++) {
      stress_tilde[i] = stress[i] ;
    }
    
    PlasticityACC_ScaleStress(plasty,stress_tilde) ;
  }

  //evaluate p and q for scaled stress_tilde
  {
    double p = (stress_tilde[0] + stress_tilde[4] + stress_tilde[8])/3. ;  
    double q = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress_tilde)) ;
    
    yield[0] = q*q*exp(k*(2*p+pc-ps)/(pc-ps)) + m*m*(p - ps)*(p + pc) ;
  }
  
  return(yield) ;
}



double* (PlasticityACC_FR)(Plasticity_t* plasty,const double* stress,const double* hardv)
{
  size_t SizeNeeded = (9+1)*(sizeof(double)) ;
  double* flow   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m  = Plasticity_GetACC_M(plasty)  ;
  double n  = Plasticity_GetACC_N(plasty)  ;
  double k  = Plasticity_GetACC_k(plasty)   ;
  double ps = Plasticity_GetInitialIsotropicTensileLimit(plasty) ;
  double beta_eps =  Plasticity_GetVolumetricStrainHardeningParameter(plasty) ;
  double pc = exp(hardv[0]) ;

  //grab parameters f_par and f_perp
  double f_par  = Plasticity_GetACC_f_par(plasty)   ;
  double f_perp = Plasticity_GetACC_f_perp(plasty)   ;
  
  //scale stress
  int axis_3  = Elasticity_GetAxis3(Plasticity_GetElasticity(plasty)) ;
  
  double stress_tilde[9] ;
  double f_scale[3] = {f_par,f_par,f_par} ;
  
  f_scale[axis_3] = f_perp;
  
  {
    int i ;
    
    for (i = 0 ; i < 9; i++) {
      stress_tilde[i] = stress[i] ;
    }
    
    PlasticityACC_ScaleStress(plasty,stress_tilde) ;
  }


  /*
    Potential function: 
    g = q*q*exp(k*(-2*p+ps-pc)/(ps-pc)) + n*n*(p - ps)*(p + pc)
  */
  {
    double p    = (stress_tilde[0] + stress_tilde[4] + stress_tilde[8])/3. ;  
    double q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(stress_tilde)) ;
    double id[9] = {1,0,0,0,1,0,0,0,1} ;
    int    i ;
    
    /*
      Plastic strain: deps^p_ij = dl * dg/dstress_ij
      ----------------------------------------------
      dp/dstress_ij = 1/3 delta_ij
      dq/dstress_ij = 3/2 dev_ij/q 
      dg/dstress_ij = 1/3 (dg/dp) delta_ij + 3/2 (dg/dq) dev_ij/q 
      dg/dp         = n*n*(2*p + pc - ps) + q*q*exp(k*(-2*p+ps-pc)/(ps-pc))*(-2k)/(ps-pc)
      dg/dq         = 2*q*exp(k*(-2*p+ps-pc)/(ps-pc))
    
      dg/dstress_ij = n*n*(2*p + pc - ps) * (1/3)* delta_ij
                    - 2k/(ps-pc)*q*q*exp(k*(-2*p+ps-pc)/(ps-pc)) * (1/3) * delta_ij
                    +            2*q*exp(k*(-2*p+ps-pc)/(ps-pc)) * (3/2) * dev_ij/q
    */
    for(i = 0 ; i < 9 ; i++) {
      double dev = stress[i] - p*id[i] ;

      flow[i]  = n*n*(2*p + pc - ps)*id[i]/3 ;
      flow[i] += (2*k/(pc-ps)*q*q*id[i]/3 + 3*dev)*exp(k*(2*p-ps+pc)/(pc-ps)) ;
    }
    
    flow[0] *= f_scale[0];
    flow[4] *= f_scale[1];
    flow[8] *= f_scale[2];
  
    /*
      The hardening flow: d(ln(pc)) = - beta_eps*deps_p
      -------------------------------------------------
      Using a = ln(pc) as hardening variable.
      d(a) = - dl * beta_eps * (dg/dp)
      So h(p,a) = - beta_eps * (dg/dp)
    */
    // condition for plastic dilation/softening    //! remove this to activate softening
    {
      double a  = (pc+ps)/(pc-ps) ;
      double k1 = k*m*m/(n*n) ;
      double k2 = (k1 == 0) ? 1 : 1 + (sqrt(1 + k1*k1*a*a) - 1)/k1 ;
      double p_star = 0.5*(ps-pc)*k2 ;
      
      if(p > p_star) {
        flow[9] = 0 ;
      } else {
        flow[9] = - beta_eps * (flow[0] + flow[4] + flow[8]) / 3 ;
      }
    }
  }

  return(flow) ;
}



void PlasticityACC_ScaleStress(Plasticity_t* plasty,double* stress)
{
  //grab parameters f_par and f_perp
  double f_par  = Plasticity_GetACC_f_par(plasty)  ;
  double f_perp = Plasticity_GetACC_f_perp(plasty) ;
  //scale stress
  int axis_3  = Elasticity_GetAxis3(Plasticity_GetElasticity(plasty)) ;
  double f_scale[3] = {f_par,f_par,f_par} ;
  
  f_scale[axis_3] = f_perp;
  
  stress[0] *= f_scale[0];
  stress[4] *= f_scale[1];
  stress[8] *= f_scale[2];
}




#undef Plasticity_GetACC_k
#undef Plasticity_GetACC_M
#undef Plasticity_GetACC_N
#undef Plasticity_GetInitialIsotropicTensileLimit
#undef Plasticity_GetACCInitialPreconsolidationPressure
#undef Plasticity_GetVolumetricStrainHardeningParameter
#undef Plasticity_GetThermalHardeningParameter
#undef Plasticity_GetACC_f_par
#undef Plasticity_GetACC_f_perp
