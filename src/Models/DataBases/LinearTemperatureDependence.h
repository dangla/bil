#ifndef LINEARTEMPERATUREDEPENDENCE_H
#define LINEARTEMPERATUREDEPENDENCE_H


#define LinearTemperatureDependence(T,T1,L1,T2,L2,...) \
        ((((T) - (T1)) * (L2) + ((T2) - (T)) * (L1)) / ((T2)-(T1)))


#endif
