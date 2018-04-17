#ifndef INTERNATIONALSYSTEMOFUNITS_H
#define INTERNATIONALSYSTEMOFUNITS_H

/*
 * Seven Base Units:
 * -----------------
 * Base quantity              base unit (symbol)
 * ---------------------------------------------
 * Length:                    meter     (m)
 * Mass:                      kilogram  (kg)
 * Time:                      second    (s)
 * Electric Current:          ampere    (A)
 * Thermodynamic Temperature: kelvin    (K)
 * Luminous Intensity:        candela   (cd)
 * Amount of Substance:       mole      (mol)
 */
 
 /*
  * Derived Units:
  * --------------
  * Angle:                    radian         (rad)     ; (m/m)
  * Solid Angle:              steradian      (sr)      ; (m2/m2)
  * Frequency:                hertz          (Hz)      ; (1/s)
  * Force:                    newton         (N)       ; (m.kg/s2)
  * Pressure:                 pascal         (Pa)      ; (N/m2)
  * Energy:                   joule          (J)       ; (N.m)
  * Power:                    watt           (W)       ; (J/s)
  * Electric Charge:          coulomb        (C)       ; (s.A)
  * Electrical Potential:     volt           (V)       ; (W/A)
  * Electric Capacitance:     farad          (F)       ; (C/V)
  * Electric Resistance:      ohm            (\Omega)  ; (V/A)
  * Electric Conductance:     siemens        (S)       ; (A/V)
  * Magnetic Flux:            weber          (Wb)      ; (V.s)
  * Magnetic Field:           tesla          (T)       ; (Wb/m2)
  * Inductance:               henry          (H)       ; (Wb/A)
  * Celsius Temperature:      degree Celsius (°C)      ; (K)
  * Luminous Flux:            lumen          (lm)      ; (cd.sr)
  * Illuminance:              lux            (lx)      ; (lm/m2)
  * Activity of radionuclide: becquerel      (Bq)      ; (1/s)
  * Aborbed Dose:             gray           (Gy)      ; (J/kg)
  * Equivalent Dose:          sievert        (Sv)      : (J/kg)
  * Catalytic Activity:       katal          (kat)     ; (mol/s)
  * 
  * 
  * (\Omega as a capital greek letter)
  */
  
  


/* vacuous declarations and typedef names */

/* class-like structure */
struct InternationalSystemOfUnits_s ; 
typedef struct InternationalSystemOfUnits_s     InternationalSystemOfUnits_t ;


extern InternationalSystemOfUnits_t* (InternationalSystemOfUnits_GetInstance)(void) ;
extern void (InternationalSystemOfUnits_Delete)(void*) ;
extern void (InternationalSystemOfUnits_UseAsLength)(const char*) ;
extern void (InternationalSystemOfUnits_UseAsTime)(const char*) ;
extern void (InternationalSystemOfUnits_UseAsMass)(const char*) ;




#define InternationalSystemOfUnits_GetMeter(si)     ((si)->meter)
#define InternationalSystemOfUnits_GetKilogram(si)  ((si)->kilogram)
#define InternationalSystemOfUnits_GetSecond(si)    ((si)->second)
#define InternationalSystemOfUnits_GetAmpere(si)    ((si)->ampere)
#define InternationalSystemOfUnits_GetKelvin(si)    ((si)->kelvin)
#define InternationalSystemOfUnits_GetCandela(si)   ((si)->candela)
#define InternationalSystemOfUnits_GetMole(si)      ((si)->mole)
#define InternationalSystemOfUnits_GetHertz(si)     ((si)->hertz)
#define InternationalSystemOfUnits_GetNewton(si)    ((si)->newton)
#define InternationalSystemOfUnits_GetPascal(si)    ((si)->pascal)
#define InternationalSystemOfUnits_GetJoule(si)     ((si)->joule)
#define InternationalSystemOfUnits_GetWatt(si)      ((si)->watt)
#define InternationalSystemOfUnits_GetCoulomb(si)   ((si)->coulomb)
#define InternationalSystemOfUnits_GetVolt(si)      ((si)->volt)
#define InternationalSystemOfUnits_GetFarad(si)     ((si)->farad)
#define InternationalSystemOfUnits_GetOhm(si)       ((si)->ohm)
#define InternationalSystemOfUnits_GetSiemens(si)   ((si)->siemens)
#define InternationalSystemOfUnits_GetWeber(si)     ((si)->weber)
#define InternationalSystemOfUnits_GetTesla(si)     ((si)->tesla)
#define InternationalSystemOfUnits_GetHenry(si)     ((si)->henry)
#define InternationalSystemOfUnits_GetDegree(si)    ((si)->degree)
#define InternationalSystemOfUnits_GetLumen(si)     ((si)->lumen)
#define InternationalSystemOfUnits_GetLux(si)       ((si)->lux)
#define InternationalSystemOfUnits_GetBecquerel(si) ((si)->becquerel)
#define InternationalSystemOfUnits_GetGray(si)      ((si)->gray)
#define InternationalSystemOfUnits_GetSievert(si)   ((si)->sievert)
#define InternationalSystemOfUnits_GetKatal(si)     ((si)->katal)
#define InternationalSystemOfUnits_GetDelete(si)    ((si)->Delete)



#define InternationalSystemOfUnits_OneMeter \
        InternationalSystemOfUnits_GetMeter(InternationalSystemOfUnits_GetInstance())
        
#define InternationalSystemOfUnits_OneKilogram \
        InternationalSystemOfUnits_GetKilogram(InternationalSystemOfUnits_GetInstance())
        
#define InternationalSystemOfUnits_OneSecond \
        InternationalSystemOfUnits_GetSecond(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneAmpere \
        InternationalSystemOfUnits_GetAmpere(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneKelvin \
        InternationalSystemOfUnits_GetKelvin(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneCandela \
        InternationalSystemOfUnits_GetCandela(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneMole \
        InternationalSystemOfUnits_GetMole(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneHertz \
        InternationalSystemOfUnits_GetHertz(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneNewton \
        InternationalSystemOfUnits_GetNewton(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OnePascal \
        InternationalSystemOfUnits_GetPascal(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneJoule \
        InternationalSystemOfUnits_GetJoule(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneWatt \
        InternationalSystemOfUnits_GetWatt(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneCoulomb \
        InternationalSystemOfUnits_GetCoulomb(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneVolt \
        InternationalSystemOfUnits_GetVolt(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneFarad \
        InternationalSystemOfUnits_GetFarad(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneOhm \
        InternationalSystemOfUnits_GetOhm(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneSiemens \
        InternationalSystemOfUnits_GetSiemens(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneWeber \
        InternationalSystemOfUnits_GetWeber(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneTesla \
        InternationalSystemOfUnits_GetTesla(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneHenry \
        InternationalSystemOfUnits_GetHenry(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneDegree \
        InternationalSystemOfUnits_GetDegree(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneLumen \
        InternationalSystemOfUnits_GetLumen(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneLux \
        InternationalSystemOfUnits_GetLux(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneBecquerel \
        InternationalSystemOfUnits_GetBecquerel(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneGray \
        InternationalSystemOfUnits_GetGray(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneSievert \
        InternationalSystemOfUnits_GetSievert(InternationalSystemOfUnits_GetInstance())

#define InternationalSystemOfUnits_OneKatal \
        InternationalSystemOfUnits_GetKatal(InternationalSystemOfUnits_GetInstance())


#include <GenericObject.h>

struct InternationalSystemOfUnits_s {
  /* Base units */
  double meter ;
  double kilogram ;
  double second ;
  double ampere ;
  double kelvin ;
  double candela ;
  double mole ;
 /* Derived Units */
  double hertz     ; /* (1/s)     */
  double newton    ; /* (m.kg/s2) */
  double pascal    ; /* (N/m2)    */
  double joule     ; /* (N.m)     */
  double watt      ; /* (J/s)     */
  double coulomb   ; /* (s.A)     */
  double volt      ; /* (W/A)     */
  double farad     ; /* (C/V)     */
  double ohm       ; /* (V/A)     */
  double siemens   ; /* (A/V)     */
  double weber     ; /* (V.s)     */
  double tesla     ; /* (Wb/m2)   */
  double henry     ; /* (Wb/A)    */
  //double degree    ; /* (K)       */
  double lumen     ; /* (cd.sr)   */
  double lux       ; /* (lm/m2)   */
  double becquerel ; /* (1/s)     */
  double gray      ; /* (J/kg)    */
  double sievert   ; /* (J/kg)    */
  double katal     ; /* (mol/s)   */
  GenericObject_Delete_t* Delete ;
} ;

#endif
