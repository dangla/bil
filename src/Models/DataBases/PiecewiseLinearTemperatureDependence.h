#ifndef PIECEWISELINEARTEMPERATUREDEPENDENCE_H
#define PIECEWISELINEARTEMPERATUREDEPENDENCE_H


//#include "Message.h"
#include "LinearTemperatureDependence.h"

#define PiecewiseLinearTemperatureDependence(T,...) \
        PiecewiseLinearTemperatureDependence_InC((T)-273.15,__VA_ARGS__)


#define PiecewiseLinearTemperatureDependence_InC(...) \
        PiecewiseLinearTemperatureDependence_Lt0(__VA_ARGS__)



#define PiecewiseLinearTemperatureDependence_Lt0(T,L1,...) \
        (((T) < 0  ) ? (L1) : \
        PiecewiseLinearTemperatureDependence_Lt25(T,L1,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt25(T,L1,L2,...) \
        (((T) < 25 ) ? LinearTemperatureDependence(T,0,L1,25,L2) : \
        PiecewiseLinearTemperatureDependence_Lt60(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt60(T,L1,L2,...) \
        (((T) < 60 ) ? LinearTemperatureDependence(T,25,L1,60,L2) : \
        PiecewiseLinearTemperatureDependence_Lt100(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt100(T,L1,L2,...) \
        (((T) < 100) ? LinearTemperatureDependence(T,60,L1,100,L2) : \
        PiecewiseLinearTemperatureDependence_Lt150(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt150(T,L1,L2,...) \
        (((T) < 150) ? LinearTemperatureDependence(T,100,L1,150,L2) : \
        PiecewiseLinearTemperatureDependence_Lt200(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt200(T,L1,L2,...) \
        (((T) < 200) ? LinearTemperatureDependence(T,150,L1,200,L2) : \
        PiecewiseLinearTemperatureDependence_Lt250(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt250(T,L1,L2,...) \
        (((T) < 250) ? LinearTemperatureDependence(T,200,L1,250,L2) : \
        PiecewiseLinearTemperatureDependence_Lt300(T,L2,__VA_ARGS__))

#define PiecewiseLinearTemperatureDependence_Lt300(T,L1,L2) \
        (((T) < 300) ? LinearTemperatureDependence(T,250,L1,300,L2) : \
        (L2))



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
