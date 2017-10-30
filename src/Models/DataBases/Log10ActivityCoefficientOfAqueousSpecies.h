#ifndef LOG10ACTIVITYCOEFFICIENTOFAQUEOUSSPECIES_H
#define LOG10ACTIVITYCOEFFICIENTOFAQUEOUSSPECIES_H


#define LOG10ACTIVITYCOEFFICIENTOFAQUEOUSSPECIES(M,Z,I) (LOG10ACTIVITYCOEFFICIENTOFAQUEOUSSPECIES_##M(Z,I))

/*
 * M   = Model (Davies,...)
 * Z   = electric charge of ion
 * I   = ionic strength = 0.5 * Sum{ m_k*Z_k*Z_k }
 * m_k = molal concentration
 */

/* Davies equation:
 * 
 * Log(g) = Z*Z*( - A*sqrt(I) / (1 + B*sqrt(I)) + C*I)
 * 
 * A = 0.5
 * B = 1
 * C = 0.15
 * 
 * Attention: I given in mol/kg (non homogeneous formula!)
 * Validity for I < 0.5 mol/kg
 * 
 * Source:
 * Davies, C.W. (1962). Ion Association. London: Butterworths. pp. 37-53
 */

#define LOG10ACTIVITYCOEFFICIENTOFAQUEOUSSPECIES_DAVIES(Z,I) ((Z)*(Z)*(-0.5*sqrt(I)/(1 + sqrt(I)) + 0.15*(I)))

#endif
