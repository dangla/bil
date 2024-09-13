#ifndef LOCALSPACETIME_H
#define LOCALSPACETIME_H


#include "Element.h"

struct LocalSpaceTime_t {
  private:
  Element_t const* _el;
  double const& _t;
  double const& _dt;
  
  public:
  /* Constructors */
  LocalSpaceTime_t(Element_t const* el,double const& t,double const& dt):
  _el(el),_t(t),_dt(dt) {}
  //LocalSpaceTime_t(void) {}
  
  /* Destructor */
  ~LocalSpaceTime_t(void) {}
  
  /* Function call operator */
  //double* operator()(int const) ;
  
  /* Accessors */
  Element_t const* GetElement(void) const {return _el;}
  double const&     GetCurrentTime(void) const {return _t;}
  double const&     GetTimeIncrement(void) const {return _dt;}
} ;


#endif
