#ifndef OBVAL_H
#define OBVAL_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct ObVal_s        ; typedef struct ObVal_s        ObVal_t ;



#define ObVal_MaxLengthOfKeyWord        (30)


#define ObVal_GetType(OV)             ((OV)->type)
#define ObVal_GetNameOfUnknown(OV)    ((OV)->inc)
#define ObVal_GetValue(OV)            ((OV)->val)
#define ObVal_GetRelaxationFactor(OV) ((OV)->relaxfactor)



#define ObVal_IsRelativeValue(OV) \
        (ObVal_GetType(OV) == 'r')

#define ObVal_IsAbsoluteValue(OV) \
        (ObVal_GetType(OV) == 'a')




struct ObVal_s {              /* Objective variation */
  char    type ;              /* Type = a(bsolute) or r(elative) */
  char*   inc ;               /* Name of the unknown */
  double  val ;               /* Objective variation */
  double  relaxfactor ;       /* Relaxation factor */
} ;

#endif
