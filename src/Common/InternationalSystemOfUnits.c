#include "Message.h"
#include "InternationalSystemOfUnits.h"
#include "Session.h"
#include "GenericData.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


//static InternationalSystemOfUnits_t* instancesi = NULL ;

static InternationalSystemOfUnits_t* InternationalSystemOfUnits_Create(void) ;
static void InternationalSystemOfUnits_UpdateDerivedUnits(InternationalSystemOfUnits_t*) ;



InternationalSystemOfUnits_t* (InternationalSystemOfUnits_Create)(void)
{
  InternationalSystemOfUnits_t* si = (InternationalSystemOfUnits_t*) malloc(sizeof(InternationalSystemOfUnits_t)) ;

  assert(si) ;
  
  
  /* Initialization */
  InternationalSystemOfUnits_GetMeter(si)    = 1 ;
  InternationalSystemOfUnits_GetKilogram(si) = 1 ;
  InternationalSystemOfUnits_GetSecond(si)   = 1 ;
  InternationalSystemOfUnits_GetAmpere(si)   = 1 ;
  InternationalSystemOfUnits_GetKelvin(si)   = 1 ;
  InternationalSystemOfUnits_GetCandela(si)  = 1 ;
  InternationalSystemOfUnits_GetMole(si)     = 1 ;
  
  InternationalSystemOfUnits_UpdateDerivedUnits(si) ;
  
  //InternationalSystemOfUnits_GetDelete(si) = InternationalSystemOfUnits_Delete ;

  return(si) ;
}


void (InternationalSystemOfUnits_Delete)(void* self)
{
  InternationalSystemOfUnits_t** pis = (InternationalSystemOfUnits_t**) self ;

  free(*pis) ;
}


#if 0
//InternationalSystemOfUnits_t* (InternationalSystemOfUnits_GetInstance0)(void)
{
  if(!instancesi) {
    instancesi = InternationalSystemOfUnits_Create() ;
  }
  
  return(instancesi) ;
}
#endif



InternationalSystemOfUnits_t*  (InternationalSystemOfUnits_GetInstance)(void)
{
  GenericData_t* gdat = Session_FindGenericData(InternationalSystemOfUnits_t,"InternationalSystemOfUnits") ;
  
  if(!gdat) {
    InternationalSystemOfUnits_t* isu = InternationalSystemOfUnits_Create() ;
    
    gdat = GenericData_Create(1,isu,InternationalSystemOfUnits_t,"InternationalSystemOfUnits") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(InternationalSystemOfUnits_t,"InternationalSystemOfUnits")) ;
  }
  
  return((InternationalSystemOfUnits_t*) GenericData_GetData(gdat)) ;
}



void InternationalSystemOfUnits_UseAsLength(const char* unit)
{
  InternationalSystemOfUnits_t* si = InternationalSystemOfUnits_GetInstance() ;
  
  if(!strcmp(unit,"nanometer")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e9 ;
  } else if(!strcmp(unit,"micrometer")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e6 ;
  } else if(!strcmp(unit,"millimeter")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e3 ;
  } else if(!strcmp(unit,"centimeter")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e2 ;
  } else if(!strcmp(unit,"decimeter")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e1 ;
  } else if(!strcmp(unit,"meter")) {
    InternationalSystemOfUnits_GetMeter(si) = 1 ;
  } else if(!strcmp(unit,"decameter")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e-1 ;
  } else if(!strcmp(unit,"hectometer")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e-2 ;
  } else if(!strcmp(unit,"kilometer")) {
    InternationalSystemOfUnits_GetMeter(si) = 1.e-3 ;
  } else if(!strcmp(unit,"inch")) {
    InternationalSystemOfUnits_GetMeter(si) = 1/2.54e-2 ;
  } else {
    arret("InternationalSystemOfUnits_UseAsLength") ;
  }
  
  InternationalSystemOfUnits_UpdateDerivedUnits(si) ;
}



void InternationalSystemOfUnits_UseAsTime(const char* unit)
{
  InternationalSystemOfUnits_t* si = InternationalSystemOfUnits_GetInstance() ;
  
  if(!strcmp(unit,"year")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/3.1536e7 ;
  } else if(!strcmp(unit,"month")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/2.592e6 ;
  } else if(!strcmp(unit,"week")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/6.048e5 ;
  }else if(!strcmp(unit,"day")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/86400. ;
  } else if(!strcmp(unit,"hour")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/3600. ;
  } else if(!strcmp(unit,"minute")) {
    InternationalSystemOfUnits_GetSecond(si) = 1/60. ;
  } else if(!strcmp(unit,"second")) {
    InternationalSystemOfUnits_GetSecond(si) = 1 ;
  } else {
    arret("InternationalSystemOfUnits_UseAsTime") ;
  }
  
  InternationalSystemOfUnits_UpdateDerivedUnits(si) ;
}



void InternationalSystemOfUnits_UseAsMass(const char* unit)
{
  InternationalSystemOfUnits_t* si = InternationalSystemOfUnits_GetInstance() ;
  
  if(!strcmp(unit,"gram")) {
    InternationalSystemOfUnits_GetKilogram(si) = 1.e3 ;
  } else if(!strcmp(unit,"decagram")) {
    InternationalSystemOfUnits_GetKilogram(si) = 1.e2 ;
  } else if(!strcmp(unit,"hectogram")) {
    InternationalSystemOfUnits_GetKilogram(si) = 1.e1 ;
  } else if(!strcmp(unit,"kilogram")) {
    InternationalSystemOfUnits_GetKilogram(si) = 1 ;
  } else if(!strcmp(unit,"tonne")) {
    InternationalSystemOfUnits_GetKilogram(si) = 1.e-3 ;
  } else {
    arret("InternationalSystemOfUnits_UseAsMass") ;
  }
  
  InternationalSystemOfUnits_UpdateDerivedUnits(si) ;
}



void InternationalSystemOfUnits_UpdateDerivedUnits(InternationalSystemOfUnits_t* si)
{
  double m   = InternationalSystemOfUnits_GetMeter(si) ;
  double kg  = InternationalSystemOfUnits_GetKilogram(si) ;
  double s   = InternationalSystemOfUnits_GetSecond(si)   ;
  double A   = InternationalSystemOfUnits_GetAmpere(si)  ;
  //double K   = InternationalSystemOfUnits_GetKelvin(si)  ;
  double cd  = InternationalSystemOfUnits_GetCandela(si) ;
  double mol = InternationalSystemOfUnits_GetMole(si)    ;
  
  double Hz  = 1/s ;
  double N   = m*kg/(s*s) ;
  double Pa  = N/(m*m) ;
  double J   = N*m ;
  double W   = J/s ;
  double C   = s*A ;
  double V   = W/A ;
  double F   = C/V ;
  double Ohm = V/A ;
  double S   = A/V ;
  double Wb  = V*s ;
  double T   = Wb/(m*m) ;
  double H   = Wb/A ;
  double lm  = cd ;
  double lx  = cd/(m*m) ;
  double Bq  = 1/s ;
  double Gy  = J/kg ;
  double Sv  = J/kg ;
  double kat = mol/s ;
  
  InternationalSystemOfUnits_GetHertz(si)     = Hz ;
  InternationalSystemOfUnits_GetNewton(si)    = N ;
  InternationalSystemOfUnits_GetPascal(si)    = Pa ;
  InternationalSystemOfUnits_GetJoule(si)     = J ;
  InternationalSystemOfUnits_GetWatt(si)      = W ;
  InternationalSystemOfUnits_GetCoulomb(si)   = C ;
  InternationalSystemOfUnits_GetVolt(si)      = V ;
  InternationalSystemOfUnits_GetFarad(si)     = F ;
  InternationalSystemOfUnits_GetOhm(si)       = Ohm ;
  InternationalSystemOfUnits_GetSiemens(si)   = S ;
  InternationalSystemOfUnits_GetWeber(si)     = Wb ;
  InternationalSystemOfUnits_GetTesla(si)     = T ;
  InternationalSystemOfUnits_GetHenry(si)     = H ;
  InternationalSystemOfUnits_GetLumen(si)     = lm ;
  InternationalSystemOfUnits_GetLux(si)       = lx ;
  InternationalSystemOfUnits_GetBecquerel(si) = Bq ;
  InternationalSystemOfUnits_GetGray(si)      = Gy ;
  InternationalSystemOfUnits_GetSievert(si)   = Sv ;
  InternationalSystemOfUnits_GetKatal(si)     = kat ;
}
