/* Modified Cam-Clay model with NFSF theory*/

/* Hao */
static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_CTNSFSHao;
static Plasticity_ReturnMapping_t                    Plasticity_RMNSFSHao;
static Plasticity_SetParameters_t                    Plasticity_SPNSFSHao ;


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




void Plasticity_SPNSFSHao(Plasticity_t* plasty,...)
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
    int i = Curves_Append(Plasticity_GetCurves(plasty),lc) ;
    
    if(i != 0) {
      Message_RuntimeError("Plasticity_SetParameters: illegal curve") ;
    }
    //Plasticity_GetLoadingCollapseFactorCurve(plasty)[0] = lc[0] ;
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      
      Plasticity_GetHardeningVariable(plasty)[0] = pc ;
      //Plasticity_GetHardeningVariable(plasty)[0] = log(pc) ;
      
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[0] = 1.e-6 ;
      Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty)[1] = 1.e-6*pc ;
      Plasticity_GetTypicalSmallIncrementOfStress(plasty) = 1.e-6*pc ;
    }
  }

  va_end(args) ;
}



double Plasticity_CTNSFSHao(Plasticity_t* plasty, const double* sig, const double* hardv, const double dlambda)
/** Modified Cam-Clay criterion */
{
    double m = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0 = Plasticity_GetInitialVoidRatio(plasty);
    //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
    double* dfsds = Plasticity_GetYieldFunctionGradient(plasty);
    double* dgsds = Plasticity_GetPotentialFunctionGradient(plasty);
    double* hm = Plasticity_GetHardeningModulus(plasty);
    double pc = hardv[0];
    double m2 = m * m;

    double id[9] = { 1,0,0,0,1,0,0,0,1 };
    double p, q, crit;

    /*
       The yield criterion
    */
    p = (sig[0] + sig[4] + sig[8]) / 3.;
    q = sqrt(3 * Math_ComputeSecondDeviatoricStressInvariant(sig));
    crit = q * q / m2 + p * (p + pc);

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

            dfsds[i] = (2 * p + pc) * id[i] / 3. + 3. / m2 * dev;
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

        hm[0] = (1 + e0) * v * p * (2 * p + pc) * pc;
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

    return(crit);
}



double Plasticity_RMNSFSHao(Plasticity_t* plasty, double* sig, double* eps_p, double* hardv)
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
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty);
    double young = Elasticity_GetYoungModulus(elasty);
    double poisson = Elasticity_GetPoissonRatio(elasty);
    double dmu = young / (1 + poisson);
    double mu = 0.5 * dmu;
    double m = Plasticity_GetSlopeCriticalStateLine(plasty);
    double kappa = Plasticity_GetSlopeSwellingLine(plasty);
    double lambda = Plasticity_GetSlopeVirginConsolidationLine(plasty);
    double e0 = Plasticity_GetInitialVoidRatio(plasty);
    //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
    double pc = hardv[0];
    double dt = hardv[1];
    double phi0 = e0 / (1 + e0);
    double m2 = m * m;
    double v = 1. / (lambda - kappa);


    /* To call the reference yield stress and  viscoplastic volumetric strain rate.
    (the lower bound point as the reference one)*/
    double CA = Plasticity_GetViscousExponent(plasty); // These parameters need to be identified from experiments.
    double epsr_pr = Plasticity_GetReferenceStrainRate(plasty);
    double p_pr = Plasticity_GetInitialPreconsolidationPressure(plasty); // Is it proper to take this one as the reference???

    double id[9] = { 1,0,0,0,1,0,0,0,1 };
    double p, q, crit;
    double p_t, q_t, Fpc_eps_pvr; // Fpc_eps_pvr is the ratio of the strain rate to the reference one.
    double dl;
    /* To call the the eps_pv value at the previous time step.*/
    double eps_pv = eps_p[0] + eps_p[4] + eps_p[8]; //viscoplastic volumetric strain
    double eps_pv_t = eps_pv;
    double eps_pvr = epsr_pr; //Take the reference one as the trial viscoplastic volumetric strain rate???

    /*
       The criterion
    */
    p = (sig[0] + sig[4] + sig[8]) / 3.;
    q = sqrt(3 * Math_ComputeSecondDeviatoricStressInvariant(sig));
    crit = q * q / m2 + p * (p + pc);

    /*
       Closest point projection algorithm.
     * Only one iterative loop is used to solve
                      q*q/m2 + p*(p + pc) = 0
       for p. The other variables (pc,q,dl) are expressed with p.
     */
    dl = 0.;
    p_t = p;
    q_t = q;
    


    if (crit > 0.) {
        double pc_n = pc;
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
            double dfsdp = 2 * p + pc;
            double dfsdq = 2 * q / m2;
            double dfsdpc = p;

            /* For the NFSF modified CamClay model,
             * the increase of pc is induced by the increase of viscoplastic volumetric strain and
             * viscoplastic volumetric strain rate. The increase of viscoplastic volumetric strain
             * should be equal to the decrease of elastic volumetric strain induced by the decrease
             * of effective mean p. During the iteration, the influence of viscoplastic
             * volumetric strain rate is contained.
             */
            double dpcsdp = -(v / (1 - phi0) + CA / (eps_pvr * dt)) * kappa * (1 - phi0) * pc / p;

            /*
             * The basis is d(dp)/p*kappa/(1+e0) - (dl*d(dg/dp) + dg/dp*d(dl)) = 0.
             * It means that the total volumetric strain keeps unchanged.
             */
            double ddlsdp = ((1 - phi0) * kappa / p - dl * (2 + dpcsdp)) / dfsdp;

            /*
             * The basis is d(dq)/(3G) - (dl*d(dg/dq) + dg/dq*d(dl)) = 0.
             * It means that the total deviatoric strain keeps unchanged.
             */
            double dqsdp = -q * 6 * mu / (m2 + 6 * mu * dl) * ddlsdp;

            /*
             * To seek the total differential of f(yield function) as a function of p,
             * including the differential of p, pc and q to p.
             */
            double df = dfsdp + dfsdq * dqsdp + dfsdpc * dpcsdp;

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


             // To calculate the viscoplastic volumetric strain rate-using p and pt to calculate viscoplastic volumetric strain
             eps_pv = eps_pv_t - (1 - phi0) * kappa * log(p_t / p); // the compression strain (total, elastic and plastic) is negative???
             
            // To judge the viscoplastic volumetric strain rate.
      Fpc_eps_pvr = (1 - phi0) * kappa * log(p_t / p) / dt / epsr_pr;
            if (eps_pvr < epsr_pr) {
                Fpc_eps_pvr = 1;
            } 
      
            // To calculate pc with the complete format
            pc = p_pr * pow(Fpc_eps_pvr, CA) * exp(-v / (1 - phi0) * eps_pv);

            dl = (1 - phi0) * kappa * log(p / p_t) / (2 * p + pc);
            q = q_t * m2 / (m2 + 6 * mu * dl);
            fcrit = q * q / m2 + p * (p + pc);

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
            double dfsds = (2 * p + pc) * id[i] / 3. + 3. / m2 * dev;

            sig[i] = p * id[i] + dev;
            eps_p[i] += dl * dfsds;

        }
    }

    /* Consolidation pressure parameters */
    hardv[0] = pc;
    hardv[1] = dt;

    /* Plastic muliplier */
    Plasticity_GetPlasticMultiplier(plasty) = dl;

    return(crit);
}


