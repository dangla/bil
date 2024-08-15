#include "BilConfig.h"

#ifdef HAVE_AUTODIFF
  #include <autodiff/forward/real/eigen.hpp>
  #include <autodiff/forward/real.hpp>

  using namespace autodiff;
#else
  using real = void;
#endif
