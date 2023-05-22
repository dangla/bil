static Damage_ComputeTangentStiffnessTensor_t   DamageMarigoJirasek_CT ;
static Damage_ReturnMapping_t                   DamageMarigoJirasek_RM ;
static Damage_SetParameters_t                   DamageMarigoJirasek_SP ;


        
#define Damage_GetCriticalEnergyReleaseRate(D) \
        Damage_GetParameter(D)[0]
        
#define Damage_GetMaximumEnergyReleaseRate(D) \
        Damage_GetParameter(D)[1]
        
#define Damage_GetUniaxialTensileStrength(D) \
        Damage_GetParameter(D)[2]
        
#define Damage_GetFractureEnergy(D) \
        Damage_GetParameter(D)[3]
        
#define Damage_GetCrackBandWidth(D) \
        Damage_GetParameter(D)[4]



void DamageMarigoJirasek_SetModelProp(Damage_t* damage)
{
  
  {
    Damage_GetComputeTangentStiffnessTensor(damage) = DamageMarigoJirasek_CT ;
    Damage_GetReturnMapping(damage)                 = DamageMarigoJirasek_RM ;
    Damage_GetSetParameters(damage)                 = DamageMarigoJirasek_SP ;
  }
  
}





void DamageMarigoJirasek_SP(Damage_t* damage,...)
{
  va_list args ;
  
  va_start(args,damage) ;
  
  {
    double ft  = va_arg(args,double) ;
    double Gf  = va_arg(args,double) ;
    double w   = va_arg(args,double) ;
    Elasticity_t* elasty = Damage_GetElasticity(damage) ;
    double E       = Elasticity_GetYoungModulus(elasty) ;
    double strain0 = ft/E ;
    double g0      = 0.5*E*strain0*strain0 ;
    
    Damage_GetUniaxialTensileStrength(damage)   = ft ;
    Damage_GetFractureEnergy(damage)            = Gf ;
    Damage_GetCrackBandWidth(damage)            = w ;
    
    Damage_GetHardeningVariable(damage)[0]      = g0 ;
    
  }

  va_end(args) ;
}



double DamageMarigoJirasek_CT(Damage_t* damage,const double* strain,const double* d,const double* hardv)
/** Energy release rate criterion
 *  Refs:
 *  [1] J.J. Marigo, Modelling of brittle and fatigue damage for elastic material 
 *  by growth of microvoids, Engineering Fracture Mechanics, 21(4):861-874, 1985.
 *  [2] M. Jirasek, B. Patzak, Consistent tangent stiffness for nonlocal damage models,
 *  Computers and Structures 80 (2002) 1279-1293.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest energy release rate (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the uniaxial tensile strength
 *  the fracture energy
 *  the width of the crack band
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  dgsde = derivative of the potential function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double* dfsde  = Damage_GetYieldFunctionGradient(damage) ;
  double* dgsde  = Damage_GetPotentialFunctionGradient(damage) ;
  double* hm     = Damage_GetHardeningModulus(damage) ;
  
  double ft      = Damage_GetUniaxialTensileStrength(damage) ;
  double Gf      = Damage_GetFractureEnergy(damage) ;
  double w       = Damage_GetCrackBandWidth(damage) ;
  double E       = Elasticity_GetYoungModulus(elasty) ;
  
  double eps0    = ft/E ;
  double gf      = Gf/w ;
  double epsf    = 0.5 * eps0 + gf/ft ;
  double kappa0  = 0.5 * E * eps0 * eps0 ;
  double kappaf  = 0.5 * E * epsf * epsf ;

  double crit ;
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    double energy = Elasticity_ComputeElasticEnergy(elasty,strain) ;
    /* G = energy release rate = - d(Free Energy)/d(d), i.e.
     * the opposite to the derivative of free energy wrt 
     * the damage variable at constant strains. 
     * Here G = 1/2 e:C0:e independant of the damage */
    double G = energy ;
    
    crit = G - kappa ;
  }
  
  /*
    Gradients
  */
  {
    int i ;
    
    Elasticity_ComputeStressTensor(elasty,strain,dfsde) ;
    
    for(i = 0 ; i < 9 ; i++) {
      dgsde[i] = dfsde[i] ;
    }
  }
  
  /* Hardening modulus: H = - d(crit)/d(d) = d(kappa)/d(d) */
  {
    double kappa = hardv[0] ;
    /* 
     * 1 - d = sqrt(kappa0/kappa) * exp(-(sqrt(kappa)-sqrt(kappa0)) / (sqrt(kappaf)-sqrt(kappa0)))
     */
    double mh = 0.5 * (1 - d[0]) / kappa * (1 + sqrt(kappa)/(sqrt(kappaf)-sqrt(kappa0))) ;
    
    hm[0] = 1/mh ;
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





double DamageMarigoJirasek_RM(Damage_t* damage,double* strain,double* d,double* hardv)
/** Energy release rate criterion
 *  Refs:
 *  [1] J.J. Marigo, Modelling of brittle and fatigue damage for elastic material 
 *  by growth of microvoids, Engineering Fracture Mechanics, 21(4):861-874, 1985.
 *  [2] M. Jirasek, B. Patzak, Consistent tangent stiffness for nonlocal damage models,
 *  Computers and Structures 80 (2002) 1279-1293.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest energy release rate (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the uniaxial tensile strength
 *  the fracture energy
 *  the width of the crack band
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  
  double ft      = Damage_GetUniaxialTensileStrength(damage) ;
  double Gf      = Damage_GetFractureEnergy(damage) ;
  double w       = Damage_GetCrackBandWidth(damage) ;
  double E       = Elasticity_GetYoungModulus(elasty) ;
  
  double eps0    = ft/E ;
  double gf      = Gf/w ;
  double epsf    = 0.5 * eps0 + gf/ft ;
  double kappa0  = 0.5 * E * eps0 * eps0 ;
  double kappaf  = 0.5 * E * epsf * epsf ;

  double G ;
  double crit ;
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    double energy = Elasticity_ComputeElasticEnergy(elasty,strain) ;
    /* G = energy release rate, i.e. here
     * the opposite to the derivative of free energy wrt 
     * the damage variable at constant strains. */
    
    G = energy ;
    
    crit = G - kappa ;
  }
  
  
  /*
    Return mapping: update damage, hardening variable
   */
  if(crit > 0.) {
    double kappa = G ;
    double d1 = 1 - sqrt(kappa0 / kappa) * exp(-(sqrt(kappa)-sqrt(kappa0)) / (sqrt(kappaf)-sqrt(kappa0))) ;

    d[0] = d1 ;
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
