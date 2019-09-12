#ifndef PIECEWISELINEARTEMPERATUREDEPENDENCE_H
#define PIECEWISELINEARTEMPERATUREDEPENDENCE_H


//#include "Message.h"
#include "LinearTemperatureDependence.h"

#define PiecewiseLinearTemperatureDependence(T,...) \
        PiecewiseLinearTemperatureDependence_InC((T)-273.15,__VA_ARGS__)


#define PiecewiseLinearTemperatureDependence_InC(...) \
        PiecewiseLinearTemperatureDependence_Lt0(__VA_ARGS__)



/* Implementation */
#define PiecewiseLinearTemperatureDependence_Lt0(T,L1,...) \
        (((T) < 0  ) ? (L1) : \
        PiecewiseLinearTemperatureDependence_Lt25(T,L1,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt25(T,L1,L2,...) \
        (((T) < 25 ) ? LinearTemperatureDependence(T,0,L1,25,L2) : \
        PiecewiseLinearTemperatureDependence_Lt60(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt60(T,L2,L3,...) \
        (((T) < 60 ) ? LinearTemperatureDependence(T,25,L2,60,L3) : \
        PiecewiseLinearTemperatureDependence_Lt100(T,L3,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt100(T,L3,L4,...) \
        (((T) < 100) ? LinearTemperatureDependence(T,60,L3,100,L4) : \
        PiecewiseLinearTemperatureDependence_Lt150(T,L4,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt150(T,L4,L5,...) \
        (((T) < 150) ? LinearTemperatureDependence(T,100,L4,150,L5) : \
        PiecewiseLinearTemperatureDependence_Lt200(T,L5,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt200(T,L5,L6,...) \
        (((T) < 200) ? LinearTemperatureDependence(T,150,L5,200,L6) : \
        PiecewiseLinearTemperatureDependence_Lt250(T,L6,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt250(T,L6,L7,...) \
        (((T) < 250) ? LinearTemperatureDependence(T,200,L6,250,L7) : \
        PiecewiseLinearTemperatureDependence_Lt300(T,L7,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt300(T,L7,L8) \
        (((T) < 300) ? LinearTemperatureDependence(T,250,L7,300,L8) : \
        (L8))



#if 0
#define PiecewiseLinearTemperatureDependence_CAT_NARG(U,...) \
        Utils_CAT_NARG(Utils_CAT(PiecewiseLinearTemperatureDependence,U),__VA_ARGS__)

#define PiecewiseLinearTemperatureDependence_Gt0_1(T,L1) \
        (L1)

#define PiecewiseLinearTemperatureDependence_Gt0_2(T,L1,L2,...) \
        (((T) < 25 ) ? LinearTemperatureDependence(T,0,L1,25,L2) : \
        PiecewiseLinearTemperatureDependence_Gt25(T,L2,__VA_ARGS__))
        
        

#define PiecewiseLinearTemperatureDependence_InC(T,...) \
        Utils_CAT_NARG(PiecewiseLinearTemperatureDependence_,__VA_ARGS__)

#define PiecewiseLinearTemperatureDependence_3(T,TA,LA,TB,LB,...) \
        (((T) < TB ) ? LinearTemperatureDependence(T,TA,LA,TB,LB,__VA_ARGS__) : \
        PiecewiseLinearTemperatureDependence_2(T,TB,LB,...))
        
#define PiecewiseLinearTemperatureDependence_2(T,TA,LA,TB,LB,...) \
        (((T) < TB ) ? LinearTemperatureDependence(T,TA,LA,TB,LB,__VA_ARGS__) : \
        PiecewiseLinearTemperatureDependence_3(T,TB,LB,...))
#endif


#endif
