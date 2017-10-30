#ifndef LOG10ACTIVITYOFWATERINBRINE_H
#define LOG10ACTIVITYOFWATERINBRINE_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Log10ActivityOfWaterInBrine_s     ; 
typedef struct Log10ActivityOfWaterInBrine_s     Log10ActivityOfWaterInBrine_t ;


extern double Log10ActivityOfWaterInBrine(const char*,double) ;


/* 
 * Dependance on the temperature
 * log(a_w(T)) = log(a_w(T0)) + h_w/R (1/T - 1/T0)
 */




#define Log10ActivityOfWaterInBrine_GetCurves(log10aw)   ((log10aw)->curves)
#define Log10ActivityOfWaterInBrine_GetNbOfSalt(log10aw) ((log10aw)->nbofsalt)


struct Log10ActivityOfWaterInBrine_s {
  int nbofsalt ;
  Curves_t* curves ;
} ;


#endif
