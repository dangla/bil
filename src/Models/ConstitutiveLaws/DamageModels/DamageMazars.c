static Damage_ComputeTangentStiffnessTensor_t    DamageMazars_CT ;
static Damage_ReturnMapping_t                    DamageMazars_RM ;
static Damage_SetParameters_t                    DamageMazars_SP ;



        
#define Damage_GetStrainAtUniaxialTensileStrength(D) \
        Damage_GetParameter(D)[0]
        
#define Damage_GetA_t(D) \
        Damage_GetParameter(D)[1]
        
#define Damage_GetB_t(D) \
        Damage_GetParameter(D)[2]
        
#define Damage_GetA_c(D) \
        Damage_GetParameter(D)[3]
        
#define Damage_GetB_c(D) \
        Damage_GetParameter(D)[4]



void DamageMazars_SetModelProp(Damage_t* damage)
{
  
  {
    Damage_GetComputeTangentStiffnessTensor(damage) = DamageMazars_CT ;
    Damage_GetReturnMapping(damage)                 = DamageMazars_RM ;
    Damage_GetSetParameters(damage)                 = DamageMazars_SP ;
      
  }
  
}
        


void DamageMazars_SP(Damage_t* damage,...)
{
  va_list args ;
  
  va_start(args,damage) ;
  
  {
    double strain0  = va_arg(args,double) ;
    double A_c  = va_arg(args,double) ;
    double A_t  = va_arg(args,double) ;
    double B_c  = va_arg(args,double) ;
    double B_t  = va_arg(args,double) ;
    
    Damage_GetStrainAtUniaxialTensileStrength(damage) = strain0 ;
    Damage_GetA_c(damage)  = A_c ;
    Damage_GetA_t(damage)  = A_t ;
    Damage_GetB_c(damage)  = B_c ;
    Damage_GetB_t(damage)  = B_t ;
    
    Damage_GetHardeningVariable(damage)[0] = strain0 ;
    
  }

  va_end(args) ;
}



double DamageMazars_CT(Damage_t* damage,const double* strain,const double* d,const double* hardv)
/** Mazars criterion.
 *  J. Mazars, A description of micro- and macroscale damage of concrete structures,
 *  Engineering Fracture Mechanics (1986), 25(5/6): 729-737.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest equivalent strain (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the Poisson's ratio (poisson),
 *  the limit elastic strain at the peak of uniaxial tensile stress test,
 *  the parameters of Mazars' model: A_t,B_t,A_c,B_c.
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  dgsde = derivative of the potential function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double* dfsde  = Damage_GetYieldFunctionGradient(damage) ;
  double* dgsde  = Damage_GetPotentialFunctionGradient(damage) ;
  double* hm     = Damage_GetHardeningModulus(damage) ;
  double kappa0  = Damage_GetStrainAtUniaxialTensileStrength(damage) ;
  double A_t     = Damage_GetA_t(damage) ;
  double B_t     = Damage_GetB_t(damage) ;
  double A_c     = Damage_GetA_c(damage) ;
  double B_c     = Damage_GetB_c(damage) ;

  double* strain_val ; //= Math_ComputePrincipalStresses(strain) ;
  double strain_dir[9] ;
  double strain_pos[3] ;
  double strain_eq ;
  double crit ;
  
  /* Principal values and principal directions */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      strain_dir[i] = strain[i] ;
    }
    
    strain_val = Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix(strain_dir,'r') ;
  }
  
  /*
    Positive principal strains
  */
  {
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_pos[i] = Math_Max(strain_val[i],0) ;
    }
  }
  
  /*
    Equivalent tensile strain
  */
  {
    int    i ;
    
    strain_eq = 0 ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_eq += strain_pos[i]*strain_pos[i] ;
    }
    
    strain_eq = sqrt(strain_eq) ;
  }
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    
    crit = strain_eq - kappa ;
  }
  
  /*
    Function gradients
    For an isotropic function of the strain tensor both
    the strain and the function gradient have the same eigenvectors.
    Denote with (e1,e2,e3) the eigenvalues of the strain tensor
    and consider the function f(e1,e2,e3): the eigenvalues
    of the function gradient are df/de1, df/de2, df/de3.
  */
  {
    int    i ;

    #define DF(i,j)   dfsde[3*(i) + (j)]
    for(i = 0 ; i < 3 ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) {
        int k ;
        
        DF(i,j) = 0 ;
        
        for(k = 0 ; k < 3 ; k++) {
          double dij = strain_dir[3*(k) + i] * strain_dir[3*(k) + j] ;
          
          DF(i,j) += strain_pos[k] / strain_eq * dij ;
        }
      }
    }
    #undef DF
    
    for(i = 0 ; i < 9 ; i++) {
      dgsde[i] = dfsde[i] ;
    }
  }
  
  /* Hardening modulus: it is the derivative of kappa w.r.t. damage */
  {
    double kappa = hardv[0] ;
    double kappa2 = kappa * kappa ;
    double dd_t = kappa0 / kappa2 * (1 - A_t) + A_t * B_t * exp(-B_t * (kappa - kappa0)) ;
    double dd_c = kappa0 / kappa2 * (1 - A_c) + A_c * B_c * exp(-B_c * (kappa - kappa0)) ;
    double strain_c[3] ;
    double strain_t[3] ;
    double strain_v = strain[0] + strain[4] + strain[8] ;
    double strain_tv = 0 ;
    double strain_cv = 0 ;
    double alpha_t = 0 ;
    double alpha_c = 0 ;
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      double stress_val = poisson / (1 - 2*poisson) * strain_v + strain_val[i] ;
      
      strain_t[i] = (stress_val > 0) ? stress_val : 0 ;
      strain_c[i] = (stress_val < 0) ? stress_val : 0 ;
      strain_tv += strain_t[i] ;
      strain_cv += strain_c[i] ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      strain_t[i] -= poisson / (1 + poisson) * strain_tv ;
      strain_c[i] -= poisson / (1 + poisson) * strain_cv ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      alpha_t += strain_t[i] * strain_pos[i] ;
      alpha_c += strain_c[i] * strain_pos[i] ;
    }
    
    alpha_t /= (strain_eq * strain_eq) ;
    alpha_c /= (strain_eq * strain_eq) ;
    
    {
      double dd = alpha_t * dd_t + alpha_c * dd_c ;
      double dd1 = (fabs(dd) > 0) ? 1./dd : 0 ;
    
      hm[0] = dd1 ;
    }
  }
  
  Damage_GetCriterionValue(damage) = crit ;

  /* 
   * Tangent matrix
   */
  {
    double* ct = Damage_GetTangentStiffnessTensor(damage) ;

    {
      double* cel = Elasticity_GetStiffnessTensor(elasty) ;
      int i ;
    
      for(i = 0 ; i < 81 ; i++) {
        ct[i] = (1 - d[0]) * cel[i] ;
      }
    }
    
    Damage_UpdateTangentStiffnessTensor(damage,ct) ;
  }
  
  return(crit) ;
}



double DamageMazars_RM(Damage_t* damage,double* strain,double* d,double* hardv)
/** Mazars return mapping.
 *  Inputs are:
 *  the strains (strain), 
 *  the largest equivalent strain (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the Poisson's ratio (poisson),
 *  the limit elastic strain at the peak of uniaxial tensile stress test,
 *  the parameters of Mazars' model: A_t,B_t,A_c,B_c.
 * 
 *  On outputs, the following values are modified:
 *  the damage (d), 
 *  the largest equivalent strain reached in material history (kappa),
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double kappa0  = Damage_GetStrainAtUniaxialTensileStrength(damage) ;
  double A_t     = Damage_GetA_t(damage) ;
  double B_t     = Damage_GetB_t(damage) ;
  double A_c     = Damage_GetA_c(damage) ;
  double B_c     = Damage_GetB_c(damage) ;

  double* strain_val ; //= Math_ComputePrincipalStresses(strain) ;
  double strain_dir[9] ;
  double strain_pos[3] ;
  double strain_eq ;
  double crit ;
  
  /* Principal values and principal directions */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      strain_dir[i] = strain[i] ;
    }
    
    strain_val = Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix(strain_dir,'r') ;
  }
  
  /*
    Positive principal strains
  */
  {
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_pos[i] = Math_Max(strain_val[i],0) ;
    }
  }
  
  /*
    Equivalent tensile strain
  */
  {
    int    i ;
    
    strain_eq = 0 ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_eq += strain_pos[i]*strain_pos[i] ;
    }
    
    strain_eq = sqrt(strain_eq) ;
  }
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    
    crit = strain_eq - kappa ;
  }
  
  
  /*
    Return mapping: update damage and hardening variable
  */
  if(crit > 0.) {
    double kappa = strain_eq ;
    double d_t = 1 - kappa0 / kappa * (1 - A_t) - A_t * exp(-B_t * (kappa - kappa0)) ;
    double d_c = 1 - kappa0 / kappa * (1 - A_c) - A_c * exp(-B_c * (kappa - kappa0)) ;
    double strain_c[3] ;
    double strain_t[3] ;
    double strain_v = strain[0] + strain[4] + strain[8] ;
    double strain_tv = 0 ;
    double strain_cv = 0 ;
    double alpha_t = 0 ;
    double alpha_c = 0 ;
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      double stress_val = poisson / (1 - 2*poisson) * strain_v + strain_val[i] ;
      
      strain_t[i] = (stress_val > 0) ? stress_val : 0 ;
      strain_c[i] = (stress_val < 0) ? stress_val : 0 ;
      strain_tv += strain_t[i] ;
      strain_cv += strain_c[i] ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      strain_t[i] -= poisson / (1 + poisson) * strain_tv ;
      strain_c[i] -= poisson / (1 + poisson) * strain_cv ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      alpha_t += strain_t[i] * strain_pos[i] ;
      alpha_c += strain_c[i] * strain_pos[i] ;
    }
    
    alpha_t /= (strain_eq * strain_eq) ;
    alpha_c /= (strain_eq * strain_eq) ;


    d[0] = alpha_t * d_t + alpha_c * d_c ;
    hardv[0] = kappa ;
  }
  
  /* Current stiffness matrix */
  {
    double* c  = Elasticity_GetStiffnessTensor(elasty) ;
    double* cd = Damage_GetDamagedStiffnessTensor(damage) ;
    int i ;
    
    for(i = 0 ; i < 81 ; i++) {
      cd[i] = (1 - d[0]) * c[i] ;
    }
  }
  
  return(crit) ;
}
