#ifndef TEMPERATURE_H
#define TEMPERATURE_H

/* class-like structures "Temperature_t" */

/* vacuous declarations and typedef names */
struct Temperature_s ; typedef struct Temperature_s Temperature_t ;



extern Temperature_t* (Temperature_GetInstance)(void) ;



#define Temperature_DefaultValue        (293.)


#define Temperature_GetRoomValue(temperature)     ((temperature)->roomvalue)

#define Temperature_SetRoomTemperature(T) \
        do {Temperature_RoomValue = (T) ;} while(0)
        
#define Temperature_RoomValue \
        Temperature_GetRoomValue(Temperature_GetInstance())


struct Temperature_s {
  double roomvalue ;          /* Room temperature */
} ;

#endif
