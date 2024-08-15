#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include "Message.h"
#include "Math_.h"
#include "Plasticity.h"
#include "autodiff.h"

static Plasticity_ComputeTangentStiffnessTensor_t    PlasticityACCPierre_CT ;
static Plasticity_ReturnMapping_t                    PlasticityACCPierre_RM ;
static Plasticity_SetParameters_t                    PlasticityACCPierre_SP ;
extern Plasticity_SetModelProp_t                     PlasticityACCPierre_SetModelProp ;


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
        
#define Plasticity_GetHydrationHardeningParameter(PL) \
        Plasticity_GetParameter(PL)[7]
        



void PlasticityACCPierre_SetModelProp(Plasticity_t* plasty)
{
  {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = PlasticityACCPierre_CT ;
    Plasticity_GetReturnMapping(plasty)                 = PlasticityACCPierre_RM ;
    Plasticity_GetSetParameters(plasty)                 = PlasticityACCPierre_SP ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 2 ;
  }
}




void PlasticityACCPierre_SP(Plasticity_t* plasty,...)
{
  va_list args ;

  va_start(args,plasty) ;


  {
    Plasticity_GetACC_k(plasty)                           = va_arg(args,double) ;
    Plasticity_GetACC_M(plasty)                           = va_arg(args,double) ;
    Plasticity_GetACC_N(plasty)                           = va_arg(args,double) ;
    Plasticity_GetInitialIsotropicTensileLimit(plasty)    = va_arg(args,double) ;
    Plasticity_GetACCInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetVolumetricStrainHardeningParameter(plasty) = va_arg(args,double) ;
    Plasticity_GetThermalHardeningParameter(plasty)          = va_arg(args,double) ;
    Plasticity_GetHydrationHardeningParameter(plasty)        = va_arg(args,double) ;
    
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



double* PlasticityACCPierre_CT(Plasticity_t* plasty,const double* sig,const double* hardv,const double* plambda)
/** Assymmetric Cam-Clay criterion
 * first hardening parameter is the preconsolidation pressure
 * second hardening parameter is the isotropic tensile elastic limit
*/
{
  double alpha   = hardv[1] ;
  double m       = alpha * Plasticity_GetACC_M(plasty)  ;
  double n       = alpha * Plasticity_GetACC_N(plasty)  ;
  double k       = Plasticity_GetACC_k(plasty)   ;
  //double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double pc0     = Plasticity_GetACCInitialPreconsolidationPressure(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double pc      = hardv[0] ;
  // double ps      = Plasticity_GetInitialIsotropicTensileLimit(plasty) ;
  double ps = pc/10.;
  double beta_eps=  alpha * Plasticity_GetVolumetricStrainHardeningParameter(plasty) ;


  double eps_v_pl = 1/beta_eps * log(pc/pc0) ; // we need here the current absolute plastic vol strain

  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;

  /*
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q*exp(k*(2*p+ps-pc)/(ps+pc))+m*m*(p+ps)*(p-pc);  //Insert here the yield function
  //crit = q*q*exp(k*(2*p+ps-pc)/(pc+ps))+m*m*(p+ps)*(p-pc);


  /*
    Gradients
    ---------
    dp/dsig_ij = 1/3 delta_ij = id[i]/3.
    dq/dsig_ij = 3/2 dev_ij/q

    further equations are taken from the jupyter notebook

    OLD:
    //df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q
    //df/dp      = 2*p + pc - ps
    //df/dq      = 2*q/m2

    //df/dsig_ij = 1/3 (2*p + pc - ps) delta_ij + (3/m2) dev_ij
  */
  {
    int    i ;

    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
      double dpsdsig = id[i]/3. ;
      //double dqsdsig = 3/2*dev/q ; //taken from existing camclay
      double dqsdsig = 3/2*dev ; //multiply already with q and remove q in the later equation

      //paste function calculated in jupyter
      dfsds[i] = (-2*dpsdsig*k*pow(q, 2) + dpsdsig*pow(m, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)) + 2*dqsdsig*(pc - ps))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps);

      //the potential gradient should be the same, only replace m by n
      dgsds[i] = (-2*dpsdsig*k*pow(q, 2) + dpsdsig*pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)) + 2*dqsdsig*(pc - ps))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps);

      //OLD:
      ////dfsds[i] = (2*p + pc - ps)*id[i]/3. + 3./m2*dev ;
      ////dgsds[i] = dfsds[i] ;
    }
  }

  /* The hardening modulus */
  /* H is defined by: df = (df/dsig_ij) dsig_ij - dl H
  for further caluclations see jupyter
   */
  {
    ////double v = 1./(lambda - kappa) ;
    double dfsdpc = -pow(m, 2)*(p + ps) + pow(q, 2)*(-k/(-pc + ps) + k*(2*p - pc + ps)/pow(-pc + ps, 2))*exp(k*(2*p - pc + ps)/(-pc + ps));
    // double dpcsdeps = pc0 * beta_eps * exp(beta_eps*eps_v_pl); //! check this !!!
    double dpcsdeps = beta_eps ; // Linear Hardening
    double dgsdp = (-2*k*pow(q, 2) + pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps);

    double pc_t = pc ;
    double p_star = (1.0/2.0)*((k - 1)*(pc_t - ps) - sqrt(pow(k, 2)*pow(pc_t, 2) + 2*pow(k, 2)*pc_t*ps + pow(k, 2)*pow(ps, 2) + pow(pc_t, 2) - 2*pc_t*ps + pow(ps, 2)))/k;
    //double softening_off = (p_t-p_star)*softening_off_switch;

    *hm = - dfsdpc * dpcsdeps * dgsdp ; // this is negative
    //*hm = MIN(- dfsdpc * dpcsdeps * dgsdp, 0.) ; //if hm is zero, no stress increments are generated, we stay at one point

    // if(p-p_star>0) {    // condition for plastic dilation/softening    //! remove this loop to activate softening
    //   *hm = 0. ; // constant yield surface
    // }

    //*hm = beta_eps*pc0*(-2*k*pow(q, 2) + pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*(-2*k*p*pow(q, 2) + pow(m, 2)*(p + ps)*pow(pc - ps, 2)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp((beta_eps*eps_v_pl*(pc - ps) - 2*k*(2*p - pc + ps))/(pc - ps))/pow(pc - ps, 3);
  }
  
  {
    double* pcrit = &crit ;
    return(pcrit) ;
  }
}



double* PlasticityACCPierre_RM(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Assymmetric Cam-Clay return mapping. Inputs are:
 *  the elastic properties K and G
 *  the plasticity parmeters,
 *  the trial preconsolidation pressure pc=hardv[0]
 *  -> this is new ! for thermal hardening, temperature effects are evaluated outside of _plasticity
 *  On outputs, the following values are modified:
 *  the stresses (sig),
 *  the plastic strains (eps_p),
 *  the pre-consolidation pressure (pc=hardv[0]).
 *
 * Braun: For the beginning we only implement the mapping with isotropic elasticity.
 * Note that this does not affect the previous function (ComputeFunctionGradients), which remains general.
 * */
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double K     = young / 3 / (1 - 2 * poisson) ;       //ADDED drained bulk modulus
  double G     = young / (2*(1+poisson)) ;             //ADDED shear modulus

  //! update anisotropic elasticity

  double alpha   = hardv[1] ;
  double m       = alpha * Plasticity_GetACC_M(plasty)  ;
  double n       = alpha * Plasticity_GetACC_N(plasty)  ;
  double k       = Plasticity_GetACC_k(plasty)   ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double pc      = hardv[0] ;
  double pc0     = Plasticity_GetACCInitialPreconsolidationPressure(plasty) ;
  // double ps      = Plasticity_GetInitialIsotropicTensileLimit(plasty) ;
  double ps = pc/10.;
  double beta_eps= alpha * Plasticity_GetVolumetricStrainHardeningParameter(plasty) ;
  //double beta_T0 =  Plasticity_GetThermalHardeningParameter(plasty) ;
  //double beta_alpha = Plasticity_GetHydrationHardeningParameter(plasty) ;


  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;

  //calculate the vol plastic strain from the previous step
  double eps_p_v_n = 0. ;
  int i ;
  for(i = 0 ; i < 9 ; i++) eps_p_v_n += eps_p[i] * id[i] ;


  /*
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;

  //this part has been moved out of _plasticity to the model:
  //double pc_t = pc0*(beta_T0*theta + 1)*exp(beta_eps*eps_p_v_n);   //the check for free thermal hardening is done outside of _plast
  double pc_t = hardv[0] ;

  pc = pc_t ;
  //check if the trial state (p,q) is within the elastic domain (crit<=0) or if we plastify (crit>0)
  //the criterion is evaluated using pc from the previous step

  // crit = q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc);  //Insert here the yield function
  crit = q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc);


  //Message_Direct("probe = %4.2e\n",p/1.e6) ;


  //! until here just copy from above

  double p_t,q_t ;
  double dl ;



  /*
    Newton-raphson return mapping
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;

  if(crit > 0.) {
    //double pc_n  = pc ;

    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.; //! set tolerance
    double R_rm[3] = {0,0,0} ;
    double x_rm[3] = {p_t,q_t,dl} ;  //give initial values: trial stresses and dlam = zero
    double J_rm[9] = {0,0,0,0,0,0,0,0,0} ;
    double dx_rm[3] = {0,0,0} ;
    double n_j[3] = {0,0,0} ; //Line Search Direction
    double theta_j = 1 ; //Step length
    double R_rmj[3] = {0,0,0} ;
    double ls_g = 0 ;
    double g_min = 0;


    while(fcrit > tol) {
      /*
      here we need the residual vector Res, with
      R[0] = f = 0
      R[1] = -p_t + p + K*dlam*dgsdp = 0
      R[2] = -q_t + q + 3*G*dlam*dgsdq = 0

      the unknown vector x contains
      x[0] = p
      x[1] = q
      x[2] = dlam

      Then we require the Jacobian J of R with respect to x

      the increment of unknowns dx is calculated by
      dx = J^-1 * R

      In the following, the Jacobian is gven in an analytical form, evaluated numerically and then inverted

      the subscript _rm is added to variable for identification
      */

      //assign x
      p = x_rm[0] ;
      q = x_rm[1] ;
      dl = x_rm[2] ;

      //calculate necessary variables and gradients
      //pc = pc_n*exp(beta_eps*(-p + p_t)/K);  //isothermal version
      //double dpcsdp = -beta_eps*pc_n*exp(beta_eps*(-p + p_t)/K)/K;//isothermal version
      //double dpcsdq = 0. ; //this is zero //isothermal version



      // pc = pc_t*exp(beta_eps*(-p + p_t)/K);
      // printf("p_t-p / K : %e ",(p_t-p)/K) ;
      // printf("pc_t : %e ",pc_t) ;
      pc = pc_t + beta_eps*(-p + p_t)/K ; //Linear Hardening
      // printf("pc : %e \n",pc) ;

      // double dpcsdp = -beta_eps*pc_t*exp(beta_eps*(-p + p_t)/K)/K;
      double dpcsdp = - beta_eps/K ; //Linear Hardening
      double dpcsdq = 0;

      //! deactivate softening by setting = 1. to activate, set = 0
      int softening_off_switch = 0 ;

      //check if softening
      double p_star = (1.0/2.0)*((k - 1)*(pc_t - ps) - sqrt(pow(k, 2)*pow(pc_t, 2) + 2*pow(k, 2)*pc_t*ps + pow(k, 2)*pow(ps, 2) + pow(pc_t, 2) - 2*pc_t*ps + pow(ps, 2)))/k;
      double softening_off = (p_t-p_star)*softening_off_switch;


      //double dgsdp_n = (-2*k*q*q + n*n*(pc_t - ps)*(2*p - pc_t + ps)*exp(k*(2*p - pc_t + ps)/(pc_t - ps)))*exp(-k*(2*p - pc_t + ps)/(pc_t - ps))/(pc_t - ps) ;



      //double dgsdp = (-2*k*q*q + n*n*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
      //double softening_off = dgsdp*softening_off_switch;

      //double softening_off = (p_t-p)*softening_off_switch;
      //double softening_off = dgsdp_n*softening_off_switch;


      if(softening_off>0) {    // condition for plastic dilation/softening and forced softening_off
         pc = pc_t ; // no yield surface shrinking
         //dpcsdp = 0. ; //force zero //! ACTIVATE ???
      }

      double dgsdp = (-2*k*q*q + n*n*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
      double dgsdq = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps)) ;
      double dgsdpc = -pow(n, 2)*(p + ps) + pow(q, 2)*(-k/(-pc + ps) + k*(2*p - pc + ps)/pow(-pc + ps, 2))*exp(k*(2*p - pc + ps)/(-pc + ps));

      //  dg² /dp /dp
      double ddgsdpdp = 4*pow(k, 2)*pow(q, 2)*exp(-2*k*p/(pc - ps) + k*pc/(pc - ps) - k*ps/(pc - ps))/pow(pc - ps, 2) + 2*pow(n, 2)*pow(pc, 2)/pow(pc - ps, 2) - 4*pow(n, 2)*pc*ps/pow(pc - ps, 2) + 2*pow(n, 2)*pow(ps, 2)/pow(pc - ps, 2);
      //  dg² /dq /dq
      double ddgsdqdq = 2*exp(-k*(2*p - pc + ps)/(pc - ps));
      //  dg² /dp /dq
      double ddgsdpdq  = -4*k*q*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
      //  dg² /dp /dpc
      double ddgsdpdpc = (-2*k*pow(q, 2) + pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*(k/(pc - ps) + k*(2*p - pc + ps)/pow(pc - ps, 2))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) + (pow(n, 2)*(pc - ps)*(-k/(pc - ps) - k*(2*p - pc + ps)/pow(pc - ps, 2))*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)) - pow(n, 2)*(pc - ps)*exp(k*(2*p - pc + ps)/(pc - ps)) + pow(n, 2)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) - (-2*k*pow(q, 2) + pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/pow(pc - ps, 2);
      //  dg² /dq /dpc
      double ddgsdqdpc = 2*q*(k/(pc - ps) + k*(2*p - pc + ps)/pow(pc - ps, 2))*exp(-k*(2*p - pc + ps)/(pc - ps));

      double dfsdp  = (-2*k*q*q + m*m*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
      double dfsdq  = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps))   ;
      double dfsdpc = -m*m*(p + ps) + q*q*(-k/(-pc + ps) + k*(2*p - pc + ps)/(-pc + ps)/(-pc + ps))*exp(k*(2*p - pc + ps)/(-pc + ps)) ;

      //calculate R
      R_rm[0] = q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc);  //copy here the yield function from above
      R_rm[1] = -p_t + p + K*dl*dgsdp ;
      R_rm[2] = -q_t + q + 3*G*dl*dgsdq ;


      //calulate J
      J_rm[0] = dfsdpc * dpcsdp + dfsdp ;//diag
      J_rm[1] = dfsdpc * dpcsdq + dfsdq ;
      J_rm[2] = 0. ;
      J_rm[3] = K * dl * (ddgsdpdpc*dpcsdp + ddgsdpdp) + 1 ;
      J_rm[4] = K * dl* (ddgsdpdpc*dpcsdq + ddgsdpdq) ; //diag
      J_rm[5] = K * dgsdp ;
      J_rm[6] = 3 * G * dl* (ddgsdqdpc * dpcsdp + ddgsdpdq) ;
      J_rm[7] = 3 * G * dl* (ddgsdqdpc * dpcsdq + ddgsdqdq) + 1 ;
      J_rm[8] = 3 * G * dgsdq ;  //diag

      // if((p_t-p)*softening_off>0) {    // condition for plastic dilation/softening and forced softening_off
      //   R_rm[1] = 0. ;
      //   J_rm[3] = 1. ;
      //   J_rm[4] = 0. ;
      //   J_rm[5] = 0. ;
      // }


      //assign R to dx, to be overwritten after the the solution dx is found by the following function
      int i ;
      int j ;
      for(i = 0 ; i < 3 ; i++) dx_rm[i] = - R_rm[i] ;

      // solves Ax=b where A is J, b is R and x is dx; the function takes J and R as arguments
      //the solution dx is written on b, therefore we have to declare first dx = R
      //Math_SolveByGaussElimination(J_rm, dx_rm, 3); //!
      //update x=dx+x
      //for(i = 0 ; i < 3 ; i++) x_rm[i] += dx_rm[i] ;  //!

      //double J_inv[9] = {0,0,0,0,0,0,0,0,0} ;
      double* J_inv = Math_Inverse3x3Matrix(J_rm) ;
      // a[(i)*2+(j)]
      // double test_2 = -J_inv[0] * R_rm[0] ;
      // double test_dq_invers = -J_inv[3] * R_rm[0] ;

      //
      // Line Search Part
      //

      for(i = 0 ; i < 3 ; i++) { // Compute search direction n_j
          n_j[i] = 0 ;
        for(j = 0 ; j < 3 ; j++) {
          n_j[i] += -J_inv[(i)*3+(j)] * R_rm[j] ;
        }
      }

      // printf("n_j = %e %e %e\n", n_j[0], n_j[1], n_j[2]);

      theta_j = 1 ; // Reset step length

      p = x_rm[0] + n_j[0] ;
      q = x_rm[1] + n_j[1] ;
      dl = x_rm[2] + n_j[2] ;
      dgsdp = (-2*k*q*q + n*n*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
      dgsdq = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps)) ;

      R_rmj[0] = ( q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc) );
      R_rmj[1] =  ( -p_t + p + K*dl*dgsdp ) ;
      R_rmj[2] = ( -q_t + q + 3*G*dl*dgsdq ) ;



      // Armijo Goldstein conditions

      double cond = 0.5*MAX(fabs(R_rmj[0]),MAX(R_rmj[1]*R_rmj[1],R_rmj[2]*R_rmj[2])) - ( 1 - 0.0001*theta_j )* 0.5*MAX(fabs(R_rm[0]),MAX(R_rm[1]*R_rm[1],R_rm[2]*R_rm[2])) ;


      if( cond > 0 ){ // Unacceptable step length, line search required
        // g_min = 0 ;
        // for(i = 0 ; i < 3 ; i++) {
        //   g_min += R_rmj[i]*R_rmj[i] ;
        // }
        g_min = MAX(fabs(R_rmj[0]),MAX(R_rmj[1]*R_rmj[1],R_rmj[2]*R_rmj[2])) ;
        // if(nf > 40){
        //   printf("step n°%d, Initial error : %e   ",nf,g_min);
        // }
        for(i = 1 ; i < 10 ; i++) {
          p = x_rm[0] + 0.1*i*n_j[0] ;
          q = x_rm[1] + 0.1*i*n_j[1] ;
          dl = x_rm[2] + 0.1*i*n_j[2] ;
          pc = pc_t + beta_eps*(-p + p_t)/K ;
          dgsdp = (-2*k*q*q + n*n*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;
          dgsdq = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps)) ;

          R_rmj[0] = ( q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc) ) ;
          R_rmj[1] = ( -p_t + p + K*dl*dgsdp ) ;
          R_rmj[2] = ( -q_t + q + 3*G*dl*dgsdq ) ;


          // ls_g = 0 ;
          // for(j = 0 ; j < 3 ; j++) {
          //   ls_g += R_rmj[j]*R_rmj[j] ;
          // }
          ls_g = MAX(fabs(R_rmj[0]),MAX(R_rmj[1]*R_rmj[1],R_rmj[2]*R_rmj[2])) ;
          // if(nf > 40){
          //   printf("theta_j : %e, Initial error : %e   ",0.1*i,ls_g);
          // }
          if( ls_g < g_min ) {
            theta_j = 0.1*i ;
            g_min = ls_g ;
          }
        }
      }

      // if(nf > 40){
      //   printf("best theta : %e, g_min : %e .\n  ",theta_j,g_min);
      //   getchar();
      // }

      //
      // End of Line Search
      //



      for(i = 0 ; i < 3 ; i++) {
        dx_rm[i] = 0 ;

        for(j = 0 ; j < 3 ; j++) {
          dx_rm[i] += -J_inv[(i)*3+(j)] * R_rm[j] ;

        }
        x_rm[i] += theta_j*dx_rm[i] ;
      }



      //update stress state
      double zero = 0.0;
      p = x_rm[0] ;
      q = x_rm[1] ;
      dl = x_rm[2] ;

      //pc = pc_n*exp(beta_eps*(-p + p_t)/K); //isothermal version
      // pc = pc_t*exp(beta_eps*(-p + p_t)/K);
      pc = pc_t + beta_eps*(-p + p_t)/K ; //Linear hardening


      if(softening_off>0) {    // condition for plastic dilation/softening and forced softening_off
         pc = pc_t ; // no yield surface shrinking
      }

      //STOP CRITERION
      //tolerances have to be without unit
      //calculate relative (dimensionless) residuals and use the maximum one as criterion
      //calculate R: //!copy the following from above
      dgsdp = (-2*k*q*q + n*n*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps) ;

      // if((p_t-p)*softening_off>0) {    // condition for plastic dilation/softening and forced softening_off
      //   dgsdp = 0. ;
      // }

      dgsdq = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps)) ;
      R_rm[0] = q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc);
      R_rm[1] = -p_t + p + K*dl*dgsdp ;
      R_rm[2] = -q_t + q + 3*G*dl*dgsdq ;
      //here we have the following dimensions:
      // R_rm[0] dimensionless
      // R_rm[1] stress
      // R_rm[2] stress

      //hence we divide the stress residuals by the trial stress (p or q), to obtain relative dimensionless values
      R_rm[0] = fabs(R_rm[0]) ;
      R_rm[1] = (R_rm[1] )*(R_rm[1] ) ; //1.e6)//!
      R_rm[2] = (R_rm[2] )*(R_rm[2] ) ; //1.e6) //!

      //relative errors, check maximum
      fcrit = MAX(R_rm[0],MAX(R_rm[1],R_rm[2])) ;

      // if(fcrit < tol) {
      //   Message_Direct("converged after %d iters\n",nf+1) ;
      // }

      // if(fcrit <= tol) {
      //   printf("Number of iterations: %i ",nf) ;
      //   printf("f = %e ",R_rm[0]) ;
      //   printf("\n") ;
      // }

      // printf(" R_f = %e ",R_rm[0]) ;
      // printf("R_p = %e ",R_rm[1]) ;
      // printf("R_q = %e ",R_rm[2]) ;
      // printf("pc = %e MPa ",pc/1.e6) ;
      // printf("p_t = %e MPa ",p_t/1.e6) ;
      // printf("q_t = %e MPa ",q_t/1.e6) ;
      // printf("hm = %e MPa \n",*hm/1.e6) ;
      //stop after too many iterations
      if(nf++ > 500) {
        printf(" R_f = %e ",R_rm[0]) ;
        printf("R_p = %e ",R_rm[1]) ;
        printf("R_q = %e ",R_rm[2]) ;
        printf("p_t = %e MPa ",p_t/1.e6) ;
        printf("q_t = %e MPa ",q_t/1.e6) ;
        printf("hm = %e MPa ",*hm/1.e6) ;
        Message_FatalError("PlasticityACCPierre_RM: no convergence") ;
      }
    }
    // printf("p = %e MPa ",p/1.e6);
    // printf("q = %e MPa \n",q/1.e6);
  }
  // if(dl < 0) {
  //   printf("dl = %e",dl) ;
  //   printf("\n") ;
  // }




  /*
    Stresses and plastic strains and yield criterion update
  */

  {
    //dgsdq = 2*q*exp(-k*(2*p - pc + ps)/(pc - ps)) ; copy from above
    // a = 1./(1 + 3*G*dl*dgsdq/q)
    //hence a = 1./(1 + 3*G*dl*(2*exp(-k*(2*p - pc + ps)/(pc - ps))))            //necessary to avoid division by 0

    double a = 1./(1 + 3*G*dl*(2*exp(-k*(2*p - pc + ps)/(pc - ps))))   ; // this is a scalar, independent on the entry of the deviator tensor
    ////double a = 1./(1 + 6*mu/m2*dl) ;
    int    i ;

    for(i = 0 ; i < 9 ; i++) {
      //decompose p and q into the stress tensor
      //go over all tensor entries
      /*the final stress state is obtained by:
      sig[i]    = p*id[i] + dev
      The value of p has been already obtained above. We need to compute dev, given by:
      dev_ij = dev_ij_t * a
      a is given above
      dev_ij_t = sig_t - p_t * I
      */

      double dev      = a*(sig[i] - p_t*id[i]) ; //!

      double dpsdsig =  id[i]/3. ;  //taken from existing camclay
      // dqsdsig = 3/2*dev/q ; //taken from existing camclay
      double dqsdsig_ =  3/2*dev ; //= dqsdsig *q
      double dgsds = (-2*dpsdsig*k*pow(q, 2) + dpsdsig*pow(n, 2)*(pc - ps)*(2*p - pc + ps)*exp(k*(2*p - pc + ps)/(pc - ps)) + 2*dqsdsig_*(pc - ps))*exp(-k*(2*p - pc + ps)/(pc - ps))/(pc - ps);
      ////double dfsds    = (2*p + pc - ps)*id[i]/3. + 3./m2*dev ;

      sig[i]    = p*id[i] + dev ;
      eps_p[i] += dl*dgsds ;    //changed to non associative
      crit = q*q*exp(k*(2*p+ps-pc)/(ps-pc))+m*m*(p+ps)*(p-pc);  //Insert here the yield function //!!!!!!!!
      double test = 0;
    }
  }

  /* Consolidation pressure */
  hardv[0] = pc ;
  
  {
    double* pcrit = &crit ;
    return(pcrit) ;
  }
}

#undef Plasticity_GetACC_k
#undef Plasticity_GetACC_M
#undef Plasticity_GetACC_N
#undef Plasticity_GetInitialIsotropicTensileLimit
#undef Plasticity_GetACCInitialPreconsolidationPressure
#undef Plasticity_GetVolumetricStrainHardeningParameter
#undef Plasticity_GetThermalHardeningParameter
#undef Plasticity_GetHydrationHardeningParameter
