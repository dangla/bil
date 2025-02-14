#ifndef LOCALSPACETIME_H
#define LOCALSPACETIME_H


#include "Element.h"

struct LocalSpaceTime_t {
  private:
  Element_t* _el;
  double _t;
  double _dt;
  
  public:
  /* Constructors */
  LocalSpaceTime_t(void) {}
  LocalSpaceTime_t(Element_t* el,double const& t,double const& dt):
  _el(el),_t(t),_dt(dt) {}
  LocalSpaceTime_t(LocalSpaceTime_t const& a) {
    _el = a.GetElement();
    _t  = a.GetCurrentTime();
    _dt = a.GetTimeIncrement();
  }
  
  /* Destructor */
  ~LocalSpaceTime_t(void) {}
  
  /* Assignement operator */
  inline LocalSpaceTime_t& operator=(LocalSpaceTime_t const& a) {
    if(this != &a) {
      _el = a.GetElement();
      _t  = a.GetCurrentTime();
      _dt = a.GetTimeIncrement();
    }
    return(*this);
  }
  
  /* Accessors */
  Element_t*        GetElement(void) const {return _el;}
  double const&     GetCurrentTime(void) const {return _t;}
  double const&     GetTimeIncrement(void) const {return _dt;}
} ;


#endif
