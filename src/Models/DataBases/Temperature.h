#ifndef TEMPERATURE_H
#define TEMPERATURE_H

/* class-like structures "Temperature_t" */

/* vacuous declarations and typedef names */
struct Temperature_s ; typedef struct Temperature_s Temperature_t ;



extern Temperature_t* (Temperature_Create)(void) ;



#define Temperature_DefaultValue        (293.)


#define Temperature_GetRoomValue(TEMP)     ((TEMP)->roomvalue)


#define Temperature_SetRoomTemperature(TEMP,T) \
        do {Temperature_GetRoomValue(TEMP) = (T) ;} while(0)


struct Temperature_s {
  double roomvalue ;          /* Room temperature */
} ;

#endif
