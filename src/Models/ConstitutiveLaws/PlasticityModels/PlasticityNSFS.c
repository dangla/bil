/* Modified Cam-Clay model with NFSF theory*/

/* Hao */
static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityNSFS_CT;
static Plasticity_ReturnMapping_t                    PlasticityNSFS_RM;
static Plasticity_YieldFunction_t                    PlasticityNSFS_YF ;
static Plasticity_FlowRules_t                        PlasticityNSFS_FR ;
static Plasticity_SetParameters_t                    PlasticityNSFS_SP ;
static Plasticity_SetModelProp_t                     PlasticityNSFS_SetModelProp ;

static double (lnxgt1_smooth)(double) ;
static double (dlnxgt1_smooth)(double) ;


#define LNXGT1(x) \
        (((x) > 1) ? log(x) : 0)
        
#define DLNXGT1(x) \
        (((x) > 1) ? 1/(x) : 0)


//#define LNPC_STRAINRATE(X,C) \
        ((C)*LNXGT1(X))
        
//#define DLNPC_STRAINRATE(X,C) \
        ((C)*DLNXGT1(X))

#define LNPC_STRAINRATE(X,C) \
        ((C)*lnxgt1_smooth(X))
        
#define DLNPC_STRAINRATE(X,C) \
        ((C)*dlnxgt1_smooth(X))


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

#define Plasticity_GetSuctionCohesionCoefficient(PL) \
        Plasticity_GetParameter(PL)[5]

#define Plasticity_GetReferenceConsolidationPressure(PL) \
        Plasticity_GetParameter(PL)[6]
        
#define Plasticity_GetViscousExponent(PL) \
        Plasticity_GetParameter(PL)[7]

#define Plasticity_GetReferenceStrainRate(PL) \
        Plasticity_GetParameter(PL)[8]
        
#define Plasticity_GetLoadingCollapseFactorCurve(PL) \
        Curves_GetCurve(Plasticity_GetCurves(PL))
        
#define Plasticity_GetSaturationDegreeCurve(PL) \
        (Curves_GetCurve(Plasticity_GetCurves(PL)) + 1)
        
        



void PlasticityNSFS_SetModelProp(Plasticity_t* plasty)
{
  
  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityNSFS_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityNSFS_RM ;
    //Plasticity_GetComputeTangentStiffnessTensor(plasty) = NULL ;
    //Plasticity_GetReturnMapping(plasty)                 = NULL ;
    Plasticity_GetYieldFunction(plasty)                 = PlasticityNSFS_YF ;
    Plasticity_GetFlowRules(plasty)                     = PlasticityNSFS_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityNSFS_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 4 ;
  }
  
}




void PlasticityNSFS_SP(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialVoidRatio(plasty)                = va_arg(args,double) ;
    Plasticity_GetSuctionCohesionCoefficient(plasty)      = va_arg(args,double) ;
    Plasticity_GetReferenceConsolidationPressure(plasty)  = va_arg(args,double) ;
    Plasticity_GetReferenceStrainRate(plasty)             = va_arg(args,double) ;
    Plasticity_GetViscousExponent(plasty)                 = va_arg(args,double) ;
    Curve_t* lc                                           = va_arg(args,Curve_t*) ;
    Curve_t* sl                                           = va_arg(args,Curve_t*) ;
    int i = Curves_Append(Plasticity_GetCurves(plasty),lc) ;
    
    if(i != 0) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    i = Curves_Append(Plasticity_GetCurves(plasty),sl) ;
    
    if(i != 1) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      
      Plasticity_GetHardeningVariable(plasty)[0] = pc ;
      //Plasticity_GetHardeningVariable(plasty)[0] = log(pc) ;
      
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6*pc ;
      //Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6 ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[1] = 1 ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[2] = 1 ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[3] = 1.e-6 ;
    }
  }

  va_end(args) ;
}



double* PlasticityNSFS_CT(Plasticity_t* plasty, const double* sig, const double* hardv, const double* plambda)
/** Modified Cam-Clay criterion */
{
    double* yield  = Plasticity_GetCriterionValue(plasty) ;
    double m      = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa  = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0     = Plasticity_GetInitialVoidRatio(plasty);
    double pc_ref = Plasticity_GetReferenceConsolidationPressure(plasty) ; 
    Curve_t* lc   = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
    Curve_t* sl   = Plasticity_GetSaturationDegreeCurve(plasty) ;
    double CA     = Plasticity_GetViscousExponent(plasty);
    double refstrainrate = Plasticity_GetReferenceStrainRate(plasty);
    
    double pc0    = hardv[0];
    double s      = hardv[1];
    double dt     = hardv[2];
    double epsv_p_rate = hardv[3];
    
    double sl_s   = Curve_ComputeValue(sl,s) ;
    double lc_s   = Curve_ComputeValue(lc,s) ;
    double ps     = sl_s*s ;
    
    double pc_b   = pc_ref * pow(pc0/pc_ref,lc_s) ;
    double pc_p   = pc_b + ps;

    double m2 = m * m;

    double* dfsds = Plasticity_GetYieldFunctionGradient(plasty);
    double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty);
    double* hm = Plasticity_GetHardeningModulus(plasty);
    //double beta    = 1 ;


    double id[9] = { 1,0,0,0,1,0,0,0,1 };
    double p, q, crit;

    /*
       The yield criterion
    */
    p = (sig[0] + sig[4] + sig[8]) / 3.;
    q = sqrt(3 * Math_ComputeSecondDeviatoricStressInvariant(sig));
    crit = q * q / m2 + p * (p + pc_p);

    /*
      Gradients
      ---------
      dp/dsig_ij = 1/3 delta_ij
      dq/dsig_ij = 3/2 dev_ij/q
      df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q
      df/dp      = 2*p + pc
      df/dq      = 2*q/m2

      df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij
    */
    {
        int    i;

        for (i = 0; i < 9; i++) {
            double dev = sig[i] - p * id[i];

            dfsds[i] = (2 * p + pc_p) * id[i] / 3. + 3. / m2 * dev;
            dgsds[i] = dfsds[i];
        }
    }

    /* The hardening modulus */
    /* H is defined by: df = (df/dsig_ij)*dsig_ij - dl*H
     * But df = (df/dsig_ij)*dsig_ij + (df/dpc_p)*dpc_p
     * Hence: H = - (df/dpc_p)*dpc_p / dl
     * On the other hand
     * dpc_p = dpc_b = pc_b * dlnpc_b
     * dlnpc_b = lc_s * dlnpc0
     * dlnpc0 = - v1 d(volstrain) + d(ln(f(volstrainrate)))
     * dlnpc0 = dl * (- v1 + df/f) * dg/dp
     * dp_c_p = dl * pc_b * lc_s * (- v1 + df/f) * dg/dp
     * Hence: H = - (df/dpc_p) * pc_b * lc_s * (- v1 + df/f) * dg/dp
     */
    {
        double v = 1. / (lambda - kappa);
        double v1 = (1 + e0)*v ;
        double df_dpc_p = p ;
        double dg_dp    = (2 * p + pc_p) ;

        /* The derivative of pc_p wrt epsv_p */
        double strainrateratio = -epsv_p_rate/refstrainrate ;
        double dstrainrateratiosdepsv_p = (dt > 0) ? -1/(dt*refstrainrate) : 0 ;
        double dlnpc0sdepsv_p   = - v1 + DLNPC_STRAINRATE(strainrateratio,CA)*dstrainrateratiosdepsv_p ;

        hm[0] = - df_dpc_p * pc_b * lc_s * dlnpc0sdepsv_p * dg_dp ;
    }

    /*
     * Tangent matrix
     */
    {
        double* c = Plasticity_GetTangentStiffnessTensor(plasty);

        Plasticity_CopyElasticTensor(plasty, c);

        Plasticity_UpdateElastoplasticTensor(plasty, c);
    }

  yield[0] = crit ;
  
  return(yield);
}


#define DEBUGPLASTICITYNSFS_RM 0

double* PlasticityNSFS_RM(Plasticity_t* plasty, double* sig, double* eps_p, double* hardv)
/** Modified NFSF Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Wang.
 *
 *  Inputs are:
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 *  the initial viscoplastic strain rate (epsv_p_rate).
 *  the time increment (dt consumed in the calculation of trial stresses).
 *
 *  On outputs, the following values are modified:
 *  the stresses (sig),
 *  the viscoplastic strains (eps_p)
 *  the pre-consolidation pressure (pc=hardv[0]).
 *
 *  Return the value of the yield function.
 **/
{
    double* yield  = Plasticity_GetCriterionValue(plasty) ;
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty);
    double bulk   = Elasticity_GetBulkModulus(elasty);
    double mu     = Elasticity_GetShearModulus(elasty);
    double m      = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa  = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0     = Plasticity_GetInitialVoidRatio(plasty);
    double pc_ref = Plasticity_GetReferenceConsolidationPressure(plasty) ; 
    Curve_t* lc   = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
    Curve_t* sl   = Plasticity_GetSaturationDegreeCurve(plasty) ;
    double CA     = Plasticity_GetViscousExponent(plasty);
    double pc_ini = Plasticity_GetInitialPreconsolidationPressure(plasty);
    double refstrainrate = Plasticity_GetReferenceStrainRate(plasty);
    
    double pc0    = hardv[0];
    double s      = hardv[1];
    double dt     = hardv[2];
    double epsv_p_rate = hardv[3];
    
    double sl_s   = Curve_ComputeValue(sl,s) ;
    double lc_s   = Curve_ComputeValue(lc,s) ;
    double ps     = sl_s*s ;

    double lnpc0    = log(pc0) ;
    double lnpc_ini = log(pc_ini) ;
    double lnpc_ref = log(pc_ref) ;
    
    double lnpc_b = lnpc_ref + lc_s*(lnpc0 - lnpc_ref) ;
    double pc_b   = exp(lnpc_b) ;
    double pc_p   = pc_b + ps;

    double m2 = m * m;
    double p, q, crit;
    double p_t, q_t;
    double dl;
  

    /*
       The criterion
    */
    p = (sig[0] + sig[4] + sig[8]) / 3.;
    q = sqrt(3 * Math_ComputeSecondDeviatoricStressInvariant(sig));
    crit = q * q / m2 + p * (p + pc_p);

    /*
       Closest point projection algorithm.
     * Only one iterative loop is used to solve
                      q*q/m2 + p*(p + pc) = 0
       for p. The other variables (pc,q,dl) are expressed with p.
     */
    dl = 0.;
    p_t = p;
    q_t = q;
  
    if (crit > 0 && dt > 0) {
        double epsv_p = eps_p[0] + eps_p[4] + eps_p[8];
        double epsv_p_n = epsv_p;
        double pc_n = pc_p;
        double fcrit = crit;
        int nf = 0;
        double tol = 1.e-8;
        double v = 1. / (lambda - kappa);
        double v1  = (1+e0)*v ;

        while (fabs(fcrit) > tol * pc_n * pc_n) {
            /*
             * Flow rule
             * ---------
             * g(p,q,pc)  = q*q*beta/m2 + (p - ps)*(p + pc)
             * dg/dp      = 2*p + pc
             * dg/dq      = (2*beta/m2)*q
             * dp/dsig_ij = 1/3*delta_ij
             * dq/dsig_ij = 3/2*dev_ij/q
             * dg/dsig_ij = 1/3*dg/dp*delta_ij + (3*beta/m2)*dev_ij
             * deps_p     = dl*dg/dp
             * deij_p     = dl*(3*beta/m2)*dev_ij
             */
            /* The partial derivative of f */
            double dfsdp = 2 * p + pc_p;
            double dfsdq = 2 * q / m2;
            double dfsdpc_p = p;

            /* The derivative of pc_p wrt p */
            double depsv_psdp  = - 1 / bulk ;
            double strainrateratio = -epsv_p_rate/refstrainrate ;
            double dstrainrateratiosdp = -depsv_psdp/(dt*refstrainrate) ;
            double dlnpc_strainratesdp  = DLNPC_STRAINRATE(strainrateratio,CA)*dstrainrateratiosdp ;
            double dlnpc0sdp_1 = - v1 * depsv_psdp ;
            double dlnpc0sdp_2 = dlnpc_strainratesdp ;
            double dlnpc0sdp   = dlnpc0sdp_1 + dlnpc0sdp_2 ;
            double dlnpc_bsdp  = lc_s * dlnpc0sdp ;
            double dpc_bsdp    = pc_b * dlnpc_bsdp ;
            double dpc_psdp    = dpc_bsdp ;
            
            /* The partial derivative of g */
            double dgsdp = 2 * p + pc_p;
            double d2gsdp2 = 2 + dpc_psdp;

            double ddlsdp = (-1 / bulk - dl * d2gsdp2) / dgsdp;

            double dqsdp = -q * 6 * mu / (m2 + 6 * mu * dl) * ddlsdp;

            double df = dfsdp + dfsdq * dqsdp + dfsdpc_p * dpc_psdp;
            
            /* Under-relaxation */
            double relax = 0.5 ;
            
            double dp = (fabs(df) > 0) ? - relax * fcrit / df : 0 ;

            p += dp;

            /* Variables (pc,dl,q) are explicit in p */
            /* Plastic multiplier (dl):
             * ------------------------
             * deps_ev  = (p - p_n)/bulk ;
             * deps_v   = (p_t - p_n) / bulk ;
             * depsv_p  = (p_t - p) / bulk = dl * dgsdp
             * Hence
             * dl = (p_t - p) / bulk / dfsdp
             * 
             * Pre-consolidation pressure (pc):
             * --------------------------------
             * epsv_p = epsv_p_n + depsv_p
             * eps_pr = depsv_p / dt
             * pc0  = pc_ini * exp(-v1*epsv_p) * (-eps_pr/refstrainrate)^CA
             * pc_b = pc_ref * pow(pc0/pc_ref,lc_s)
             * pc_p = pc_b + ps
             * 
             * Deviatoric behavior (q):
             * ------------------------
             * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
             * Hence
             * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
             * q      = q_t / (1 + 6 mu / m2 dl)
             */
            
            // The volumetric plastic strain 
            epsv_p = epsv_p_n + (p_t - p) / bulk ;
             
            // The volumetric plastic strain rate
            epsv_p_rate = (epsv_p - epsv_p_n) / dt;
            
            // The strain rate ratio
            strainrateratio = -epsv_p_rate/refstrainrate ;
            
            // The consolidation pressure at 0 suction
            {
              double lnpc_strainrate  = LNPC_STRAINRATE(strainrateratio,CA) ;

              lnpc0  = lnpc_ini - v1 * epsv_p + lnpc_strainrate ;
              pc0    = exp(lnpc0) ;
            }

            // The consolidation pressure pc_p
            lnpc_b = lnpc_ref + lc_s*(lnpc0 - lnpc_ref) ;
            pc_b   = exp(lnpc_b) ;
            pc_p   = pc_b + ps;

            // The plastic multiplier dl, q and fcrit
            dl    = -(p - p_t) / bulk / (2 * p + pc_p);
            q     = q_t * m2 / (m2 + 6 * mu * dl);
            fcrit = q * q / m2 + p * (p + pc_p);
            
            #if DEBUGPLASTICITYNSFS_RM
            {
              printf("\n") ;
              printf("=========\n") ;
              printf("iter = %d\n",nf) ;
              printf("=========\n") ;
              printf("p = %e\n",p) ;
              printf("pc0 = %e\n",pc0) ;
              printf("dl = %e\n",dl) ;
              printf("q = %e\n",q) ;
              printf("fcrit = %e\n",fcrit) ;
              printf("df = %e\n",df) ;
              printf("strainrateratio = %e\n",strainrateratio) ;
            }
            #endif

            if (nf++ > 30) {
                Message_FatalError("PlasticityNSFS_RM: no convergence");
            }
        }
    }

    /*
      Stresses and plastic strains
    */

    {
        double a = m2 / (m2 + 6 * mu * dl) ;
        double id[9] = { 1,0,0,0,1,0,0,0,1 };
        int    i;

        for (i = 0; i < 9; i++) {
            double dev = a * (sig[i] - p_t * id[i]);
            double dfsds = (2 * p + pc_p) * id[i] / 3. + 3. / m2 * dev;

            sig[i] = p * id[i] + dev;
            eps_p[i] += dl * dfsds;
        }
    }

    /* Consolidation pressure parameters */
    hardv[0] = pc0;
                                                         

    /* Plastic muliplier */
    Plasticity_GetPlasticMultiplier(plasty)[0] = dl;
    
    yield[0] = crit ;
    
    return(yield);
}


#if 0 // version Hao
double* PlasticityNSFS_RM(Plasticity_t* plasty, double* sig, double* eps_p, double* hardv)
/** Modified NFSF Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Wang.
 *
 *  Inputs are:
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 *  the initial viscoplastic strain rate (epsr_pr).
 *  the time increment (dt consumed in the calculation of trial stresses).
 *
 *  On outputs, the following values are modified:
 *  the stresses (sig),
 *  the viscoplastic strains (eps_pv),
 *  the viscoplastic strain rates (eps_pvr),
 *  the time increment (dt),
 *  the pre-consolidation pressure (pc=hardv[0]).
 *
 *  Return the value of the yield function.
 **/
{
    double* yield  = Plasticity_GetCriterionValue(plasty) ;
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty);
    double young = Elasticity_GetYoungModulus(elasty);
    double poisson = Elasticity_GetPoissonRatio(elasty);
    double mu = 0.5 * young / (1 + poisson);
    double bulk  =  young / (1 - 2*poisson)/3;
    double m = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0 = Plasticity_GetInitialVoidRatio(plasty);
        
    double p_c     = Plasticity_GetReferenceConsolidationPressure(plasty) ; 
    Curve_t* lc    = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
    Curve_t* sl    = Plasticity_GetSaturationDegreeCurve(plasty) ;
    //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
    double pc0 = hardv[0]; /* pc0 */
    double s  = hardv[1];
    double dt = hardv[2];

    /* To call the reference yield stress and  viscoplastic volumetric strain rate.
    (the lower bound point as the reference one)*/
    double CA = Plasticity_GetViscousExponent(plasty); // The viscous parameter.
    double epsr_pr = Plasticity_GetReferenceStrainRate(plasty);
    double p_pr = Plasticity_GetInitialPreconsolidationPressure(plasty); /* The reference Pressure is taken as the initial Preconsolidation one */

    /* To call the the eps_pv qnd epsr_pvr values at the previous time step.*/
    double eps_pv = eps_p[0] + eps_p[4] + eps_p[8]; //viscoplastic volumetric strain                                                  
    double epsr_pvr = - hardv[3]; //viscoplastic volumetric strain.
    
    double sl_s = Curve_ComputeValue(sl,s) ;
    double lc_s = Curve_ComputeValue(lc,s) ;
    double ps   = sl_s*s ;

    double phi0 = e0 / (1 + e0);
    double m2 = m * m;
    double v = 1. / (lambda - kappa);
    //double beta    = 1 ;

    //double pcs = p_c * pow(pc0/p_c*exp(-v / (1 - phi0)*dt*epsr_pvr),lc_s) ;//net yield stress without viscoplastic influence.
    double pcs = p_c * pow(pc0/p_c,lc_s) ;//net yield stress without viscoplastic influence.
    
    if (epsr_pvr < epsr_pr) {
                epsr_pvr = epsr_pr;
    } 
    
    double pc_b = pcs * pow(epsr_pvr/epsr_pr,CA);
    double pc_p = pc_b + ps;
  



    double id[9] = { 1,0,0,0,1,0,0,0,1 };
    double p, q, crit;
    double p_t, q_t;
    double dl;
  

    /*
       The criterion
    */
    p = (sig[0] + sig[4] + sig[8]) / 3.;
    q = sqrt(3 * Math_ComputeSecondDeviatoricStressInvariant(sig));
    crit = q * q / m2 + p * (p + pc_p);

    /*
       Closest point projection algorithm.
     * Only one iterative loop is used to solve
                      q*q/m2 + p*(p + pc) = 0
       for p. The other variables (pc,q,dl) are expressed with p.
     */
    dl = 0.;
    p_t = p;
    q_t = q;
  


    if (crit > 0. && dt > 0) {
        
        double eps_pv_n = eps_pv;
        double pc_n = pc_p;
                                                            
        double fcrit = crit;
        int nf = 0;
        double tol = 1.e-8;


        while (fabs(fcrit) > tol * pc_n * pc_n) {
            /*
             * Flow rule
             * ---------
             * df/dp      = 2*p + pc
             * df/dq      = (2/m2) q
             * dp/dsig_ij = 1/3 delta_ij
             * dq/dsig_ij = 3/2 dev_ij/q
             * df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij
             * deps_p     = dl (2*p + pc)
             * deij_p     = dl (3/m2) dev_ij
             */
            double dfsdp = 2 * p + pc_p;
            double dfsdq = 2 * q / m2;
            double dfsdpc_p = p;

            /* For the NFSF modified BBM model,
             * the increase of pc is induced by the increase of viscoplastic volumetric strain and
             * viscoplastic volumetric strain rate. The increase of viscoplastic volumetric strain
             * should be equal to the decrease of elastic volumetric strain induced by the decrease
             * of effective mean p. During the iteration, the influence of viscoplastic
             * volumetric strain rate is contained.
             */     
      
      // i.e., dpc_r-dp        
            double dpc_psdp = pc_b * ( lc_s * v / (1 - phi0)  + CA / Math_Max(dt*epsr_pr, dt*epsr_pvr ) )/ bulk; 
      

            /*
             * The basis is d(dp)/K - (dl*d(dg/dp) + dg/dp*d(dl)) = 0.
             * It means that the total volumetric strain keeps unchanged.
             */
            double ddlsdp = (-1 / bulk - dl * (2 + dpc_psdp)) / dfsdp;

            /*
             * The basis is d(dq)/(3G) - (dl*d(dg/dq) + dg/dq*d(dl)) = 0.
             * It means that the total deviatoric strain keeps unchanged.
             */
            double dqsdp = -q * 6 * mu / (m2 + 6 * mu * dl) * ddlsdp;

            /*
             * To seek the total differential of f(yield function) as a function of p,
             * including the differential of p, pc and q to p.
             */
            double df = dfsdp + dfsdq * dqsdp + dfsdpc_p * dpc_psdp;

            /*
             * To correct the drift of p.
             */
            p -= fcrit / df;


            /* Variables (pc,dl,q) are explicit in p */
            /* Plastic multiplier (dl):
             * ------------------------
             * deps_ev = (p-p_n)/bulk ;
             * deps_v   = (p_t - p_n) / bulk ;
             * deps_pv = (p_t - p) / bulk = dl (2*p + pc)
             * Hence
             * dl = (p_t - p) / bulk / (2*p + pc)
             * 
             * Pre-consolidation pressure (pc):
             * --------------------------------
             * deps_pv = - (1-phi0) (lambda - kappa)*ln(pc0/pc0_n)
             * eps_pv = depsv_p + eps_pv_n
             * eps_pr = - deps_p / dt
             * Hence
             * pc0 = p_pr * exp(-v/(1 - phi0)*eps_pv)
             * pc_b = p_c * pow(pc0/p_c,lc_s)* (eps_pr/epsr_pr)^CA ;
       * pc_p = pc_b + ps;
             * Deviatoric behavior (q):
             * ------------------------
             * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
             * Hence
             * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
             * q      = q_t / (1 + 6 mu / m2 dl)
             */
            
            // To calculate the viscoplastic volumetric strain using p and pt 
            eps_pv = eps_pv_n + (p_t - p) / bulk; // the compression strain (total, elastic and plastic) is negative
             
            // To judge the viscoplastic volumetric strain rate.
            epsr_pvr = -(p_t - p) / bulk / dt;                                                       
            if (epsr_pvr < epsr_pr) {
                epsr_pvr = epsr_pr;
            } 
            

            // To calculate pc0 with the complete equation for pc_p
            pc0 = p_pr  * exp(-v / (1 - phi0) * eps_pv);
      pc_b = p_c * pow(pc0/p_c,lc_s) * pow(epsr_pvr / epsr_pr, CA);
      pc_p = pc_b + ps;

            dl = -(p - p_t) / bulk / (2 * p + pc_p);
            q = q_t * m2 / (m2 + 6 * mu * dl);
            fcrit = q * q / m2 + p * (p + pc_p);

            if (nf++ > 500) {
                Message_FatalError("ViscoPlasticity_ReturnMappingCamClay: no convergence");
            }
        }
    }

    /*
      Stresses and plastic strains
    */

    {
        double a = 1. / (1 + 6 * mu / m2 * dl);
        int    i;

        for (i = 0; i < 9; i++) {
            double dev = a * (sig[i] - p_t * id[i]);
            double dfsds = (2 * p + pc_p) * id[i] / 3. + 3. / m2 * dev;

            sig[i] = p * id[i] + dev;
            eps_p[i] += dl * dfsds;

        }
    }

    /* Consolidation pressure parameters */
    hardv[0] = pc0;
                                                         

    /* Plastic muliplier */
    Plasticity_GetPlasticMultiplier(plasty)[0] = dl;
    
    yield[0] = crit ;
    
    return(yield);
}
#endif




double* (PlasticityNSFS_YF)(Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the value of the yield functions. 
 **/
{
  size_t SizeNeeded = sizeof(double) ;
  double* yield   = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  double m        = Plasticity_GetSlopeCriticalStateLine(plasty);
  double pc_ref   = Plasticity_GetReferenceConsolidationPressure(plasty) ; 
  Curve_t* lc     = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
  Curve_t* sl     = Plasticity_GetSaturationDegreeCurve(plasty) ;
  double pc0      = hardv[0];
  double s        = hardv[1];
  double sl_s     = Curve_ComputeValue(sl,s) ;
  double lc_s     = Curve_ComputeValue(lc,s) ;
  double ps       = sl_s*s ;
  double pc_b     = pc_ref * pow(pc0/pc_ref,lc_s) ;
  double pc_p     = pc_b + ps;
  double p        = (stress[0] + stress[4] + stress[8])/3.;
  double q2       = 3*Math_ComputeSecondDeviatoricStressInvariant(stress) ;
  double m2       = m*m ;
  
  yield[0] = q2/m2 + p*(p + pc_p);

  return(yield) ;
}



double* (PlasticityNSFS_FR)(Plasticity_t* plasty,const double* sig,const double* hardv)
/** 
 **/
{
  size_t SizeNeeded = (9+4)*(sizeof(double)) ;
  double* flow    = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double bulk     = Elasticity_GetBulkModulus(elasty) ;
  double m        = Plasticity_GetSlopeCriticalStateLine(plasty);
  double kappa    = Plasticity_GetSlopeSwellingLine(plasty);
  double lambda   = Plasticity_GetSlopeVirginConsolidationLine(plasty);
  double e0       = Plasticity_GetInitialVoidRatio(plasty);
  double pc_ref   = Plasticity_GetReferenceConsolidationPressure(plasty) ;
  Curve_t* lc     = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
  Curve_t* sl     = Plasticity_GetSaturationDegreeCurve(plasty) ;
  double CA       = Plasticity_GetViscousExponent(plasty);
  double pc_ini   = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double refstrainrate  = Plasticity_GetReferenceStrainRate(plasty);
  double pc0      = hardv[0];
  double s        = hardv[1];
  double dt       = hardv[2];
  double epsv_p_rate = hardv[3];
  double sl_s     = Curve_ComputeValue(sl,s) ;
  double lc_s     = Curve_ComputeValue(lc,s) ;
  double ps       = sl_s*s ;
  double pc_b     = pc_ref * pow(pc0/pc_ref,lc_s) ;
  double pc_p     = pc_b + ps;
  double m2       = m*m;
  
  /*
   * Potential function: g(p,q,pc) = q*q/m2 + p*(p + pc)
   */
  
  /*
   * Plastic strain flow:
   * --------------------
   * deps^p_ij     = dl * dg/dstress_ij
   * dp/dstress_ij = 1/3 delta_ij
   * dq/dstress_ij = 3/2 dev_ij/q 
   * dg/dstress_ij = 1/3 (dg/dp) delta_ij + 3/2 (dg/dq) dev_ij/q 
   * dg/dp         = 2*p + pc
   * dg/dq         = 2*q/m2
   * 
   * dg/dstress_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij 
   */

  {
    double N        = 4/(pc_ini*pc_ini*bulk) ;
    double p        = (sig[0] + sig[4] + sig[8]) / 3. ;
    double id[9]    = { 1,0,0,0,1,0,0,0,1 };
    int    i;

    for (i = 0; i < 9; i++) {
      double dev = sig[i] - p * id[i];

      flow[i] = N*((2 * p + pc_p) * id[i] / 3. + 3. / m2 * dev) ;
    }
  }
  
  /*
   * The hardening flows:
   * --------------------
   * d(ln(pc0)) = - (1 + e0)*v*d(volstrain) + d(ln(f(volstrainrate)))
   * Using a = ln(pc0) as hardening variable.
   * d(a) = - dl * (1 + e0) * v * (dg/dp) + dl * df/f * 1/dt * (dg/dp)
   * So h(p,a) = (- (1 + e0) * v + df/f * 1/dt) * (dg/dp)
   * d(volstrainrate) = 1/dt d(volstrain - volstrain_n) = 1/dt d(volstrain)
   */
  {
    double v = 1./(lambda - kappa) ;
    double v1  = (1+e0)*v ;
    double strainrateratio = -epsv_p_rate/refstrainrate ;
    /* The derivative of lnpc0 wrt epsv_p */
    double dstrainrateratiosdepsv_p = (dt > 0) ? -1/(dt*refstrainrate) : 0 ;
    double dlnpc0sdepsv_p   = - v1 + DLNPC_STRAINRATE(strainrateratio,CA)*dstrainrateratiosdepsv_p ;
    double h = flow[0] + flow[4] + flow[8] ;

    flow[9]  = pc0*dlnpc0sdepsv_p*h ;
    flow[10] = 0 ;
    flow[11] = 0 ;
    flow[12] = (dt > 0) ? h/dt : 0 ; // Not correct -> doesn't work at present
  }
  
  return(flow) ;
}


#define A 0.2

double (lnxgt1_smooth)(double x)
{
  double a = A ;
  double y = 1 + exp((x-1)/a) ;
  
  return(log(1 + a*log(y))) ;
}
double (dlnxgt1_smooth)(double x)
{
  double a = A ;
  double y = 1 + exp((x-1)/a) ;
  
  return((1 - 1/y)/(1 + a*log(y))) ;
}
#undef A



#undef Plasticity_GetSlopeSwellingLine
#undef Plasticity_GetSlopeVirginConsolidationLine
#undef Plasticity_GetSlopeCriticalStateLine
#undef Plasticity_GetInitialPreconsolidationPressure
#undef Plasticity_GetInitialVoidRatio
#undef Plasticity_GetSuctionCohesionCoefficient
#undef Plasticity_GetReferenceConsolidationPressure
#undef Plasticity_GetViscousExponent
#undef Plasticity_GetReferenceStrainRate
#undef Plasticity_GetLoadingCollapseFactorCurve
#undef Plasticity_GetSaturationDegreeCurve
