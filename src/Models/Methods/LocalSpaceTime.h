#ifndef LOCALSPACETIME_H
#define LOCALSPACETIME_H


#include "Element.h"

struct LocalSpaceTime_t {
  private:
  const Element_t* _el;
  const double& _t;
  const double& _dt;
  
  public:
  /* Constructors */
  LocalSpaceTime_t(const Element_t* el,const double& t,const double& dt):
  _el(el),_t(t),_dt(dt) {}
  //LocalSpaceTime_t(void) {}
  
  /* Destructor */
  ~LocalSpaceTime_t(void) {}
  
  /* Function call operator */
  //double* operator()(const int) ;
  
  /* Accessors */
  Element_t* GetElement(void) const {return _el;}
  double     GetCurrentTime(void) const {return _t;}
  double     GetTimeIncrement(void) const {return _dt;}
} ;


#endif
