/* Modified Cam-Clay model with NFSF theory*/

/* Hao */
static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityNSFS_CT;
static Plasticity_ReturnMapping_t                    PlasticityNSFS_RM;
static Plasticity_SetParameters_t                    PlasticityNSFS_SP ;
static Plasticity_SetModelProp_t                     PlasticityNSFS_SetModelProp ;


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
    //Plasticity_GetYieldFunction(plasty)                 = PlasticityNSFS_YF ;
    //Plasticity_GetFlowRules(plasty)                     = PlasticityNSFS_FR ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityNSFS_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 3 ;
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
      
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
    }
  }

  va_end(args) ;
}



double* PlasticityNSFS_CT(Plasticity_t* plasty, const double* sig, const double* hardv, const double* plambda)
/** Modified Cam-Clay criterion */
{
  double* yield  = Plasticity_GetCriterionValue(plasty) ;
    double m = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0 = Plasticity_GetInitialVoidRatio(plasty);
    //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
    double* dfsds = Plasticity_GetYieldFunctionGradient(plasty);
    double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty);
    double* hm = Plasticity_GetHardeningModulus(plasty);
  
  double p_c     = Plasticity_GetReferenceConsolidationPressure(plasty) ;
  Curve_t* lc    = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
    Curve_t* sl    = Plasticity_GetSaturationDegreeCurve(plasty) ;
  
    double pc0 = hardv[0]; /* pc0 */
  double s  = hardv[1];
  double sl_s = Curve_ComputeValue(sl,s) ;
    double lc_s = Curve_ComputeValue(lc,s) ;
  double ps   = sl_s*s ;
  double pc_b = p_c * pow(pc0/p_c,lc_s);
  double pc_p = pc_b + ps; 
    double m2 = m*m;
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
    /* H is defined by: df = (df/dsig_ij) dsig_ij - dl H
     * But df = (df/dsig_ij) dsig_ij + (df/dpc) dpc
     * Hence: H = - (df/dpc) dpc / dl
     * On the other hand
     * dpc/pc = - v de_p = - (1 + e0) v deps_p = - (1 + e0) v dl (dg/dp)
     * i.e. dpc = dl hpc with hpc = - (1 + e0)*v*pc*(2*p + pc)
     * Hence: H = (1 + e0) v (df/dpc) (dg/dp) pc
     */
    {
        double v = 1. / (lambda - kappa);
        double v1 = (1 + e0)*v ;
        double h = - v1*(2 * p + pc_p) ;
        double dpcda    = pc_b*lc_s ;
        double dfda     = p*dpcda ;
        hm[0] = - dfda * h ; // It needs to be further checked.
    }

    /*
     * Tangent matrix
     */
    {
        double* c = Plasticity_GetTangentStiffnessTensor(plasty);

        Plasticity_CopyElasticTensor(plasty, c);


        /*
        if (dlambda > 0) {
            double v = 1. / (lambda - kappa);
            double v1 = (1 + e0) * v;
            double hpc = -v1 * (2 * p + pc) * pc;
            double dhpcdpc = -v1 * 2 * (p + pc);
            double dhpcdp = -v1 * 2 * pc;
            double dfdpc = p;
            double ddgdpdpc = 1;
            double dlambda1 = dlambda / (1 - dlambda * dhpcdpc);
            Elasticity_t* elasty = Plasticity_GetElasticity(plasty);
            double bulk = Elasticity_GetBulkModulus(elasty);
            double shear = Elasticity_GetShearModulus(elasty);
            double g0 = shear;
            double k0 = bulk;
            double g1 = g0 / (6 * g0 * dlambda / m2 + 1);
            double k1 = k0 / (2 * k0 * dlambda * (1 - dlambda1 * v1 * pc) + 1);

            /*
            {
                double mu = g1;
                double lame = k1 - 2. / 3. * mu;
                int    i;

                for (i = 0; i < 81; i++) c[i] = 0.;

                for (i = 0; i < 3; i++) {
                    int j;

                    #define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
                    for (j = 0; j < 3; j++) {
                        C(i, i, j, j) += lame;
                        C(i, j, i, j) += mu;
                        C(i, j, j, i) += mu;
                    }
                    #undef C
                }
            }


            {
                int    i;

                for (i = 0; i < 9; i++) {
                    dfsds[i] += dlambda1 * dfdpc * dhpcdp * id[i] / 3;
                    dgsds[i] += dlambda1 * hpc * ddgdpdpc * id[i] / 3;
                }
            }

            hm[0] /= (1 - dlambda * dhpcdpc);
        }
        */

        Plasticity_UpdateElastoplasticTensor(plasty, c);
    }

  yield[0] = crit ;
  
  return(yield);
}



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
    double sl_s = Curve_ComputeValue(sl,s) ;
    double lc_s = Curve_ComputeValue(lc,s) ;
  double ps   = sl_s*s ;
  double pc_b = p_c * pow(pc0/p_c,lc_s);
  double pc_p = pc_b + ps;  
  
    double phi0 = e0 / (1 + e0);
    double m2 = m * m;
    double v = 1. / (lambda - kappa);
    //double beta    = 1 ;


    /* To call the reference yield stress and  viscoplastic volumetric strain rate.
    (the lower bound point as the reference one)*/
    double CA = Plasticity_GetViscousExponent(plasty); // These parameters need to be identified from experiments.
    double epsr_pr = Plasticity_GetReferenceStrainRate(plasty);
    double p_pr = Plasticity_GetInitialPreconsolidationPressure(plasty); /* The initial Preconsolidation Pressure  
  is taken as the reference one */

    double id[9] = { 1,0,0,0,1,0,0,0,1 };
    double p, q, crit;
    double p_t, q_t;
    double dl;
  
    /* To call the the eps_pv value at the previous time step.*/
    double eps_pv = eps_p[0] + eps_p[4] + eps_p[8]; //viscoplastic volumetric strain
    double eps_pv_t = eps_pv;
    double epsr_pvr = epsr_pr; //Take the reference strain rate as the trial one.

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
        double pc_n = pc_p; // pc_p = pc in Camclay
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
             
            double dpc_psdp = pc_b * lc_s * ( - v / (1 - phi0)  + CA / Math_Max(dt*epsr_pvr, eps_pv - eps_pv_t) )/ bulk; // i.e., dpc_r-dp


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
             * deps_e = - (1-phi0) kappa ln(p/p_n) ;
             * deps   = - (1-phi0) kappa ln(p_t/p_n) ;
             * deps_p = - (1-phi0) kappa ln(p_t/p) = dl (2*p + pc)
             * Hence
             * dl = (1-phi0) kappa ln(p/p_t) / (2*p + pc)
             * pc =
             * Pre-consolidation pressure (pc):
             * --------------------------------
             * deps_p = - (1-phi0) kappa ln(p_t/p) = - (1-phi0) (lambda - kappa) ln(pc/pc_n)
             * eps_p = deps_p + eps_p_n
             * eps_pr = - deps_p / dt
             * Hence
             * pc = p_pr * (eps_pr/epsr_pr)^CA * exp(-v/(1 - phi0)*eps_p)
             * Deviatoric behavior (q):
             * ------------------------
             * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
             * Hence
             * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
             * q      = q_t / (1 + 6 mu / m2 dl)
             */


             // To calculate the viscoplastic volumetric strain using p and pt 
             eps_pv = eps_pv_t + (p - p_t) / bulk; // the compression strain (total, elastic and plastic) is negative
             
            // To judge the viscoplastic volumetric strain rate.
            epsr_pvr = (p - p_t) / bulk / dt;
            if (epsr_pvr < epsr_pr) {
                epsr_pvr = epsr_pr;
            } 
      
            // To calculate pc0 with the complete format
            pc0 = p_pr * pow(epsr_pvr / epsr_pr, CA) * exp(-v / (1 - phi0) * eps_pv);
          pc_b = p_c * pow(pc0/p_c,lc_s);
          pc_p = pc_b + ps;

            dl = (p - p_t) / bulk / (2 * p + pc_p);
            q = q_t * m2 / (m2 + 6 * mu * dl);
            fcrit = q * q / m2 + p * (p + pc_p);

            if (nf++ > 20) {
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
