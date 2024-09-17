#ifndef CUSTOMVALUES_H
#define CUSTOMVALUES_H

#include "Message.h"


template<class... C>
struct CustomValues_t: C... {
  #include "CustomValues_MemberOperations.in"
};


/* Partial template specialization */
template<class IM,class EX,class CO,class... OT>
struct CustomValues_t<IM,EX,CO,OT...>: IM,EX,CO,OT... {
  using ImplicitValues_type = IM;
  using ExplicitValues_type = EX;
  using ConstantValues_type = CO;
  
  #include "CustomValues_MemberOperations.in"
};


/* Convert a pointer V to any type in pointer to type T */
#define CustomValues_Convert(V,T) (((T)*) (V))

#define CustomValues_Index(CV,V,T)  ((int) (((T*) (&(CV)->V)) - ((T*) (CV))))

#define CustomValues_TypeOfImplicitValues(CV)  typename CV::ImplicitValues_type
#define CustomValues_TypeOfExplicitValues(CV)  typename CV::ExplicitValues_type
#define CustomValues_TypeOfConstantValues(CV)  typename CV::ConstantValues_type


#define CustomValues_NbOfMembers(CV,T)  ((int) (sizeof(CV)/sizeof(T)))
#define CustomValues_Size(CV,T)  ((int) (sizeof(CV)/sizeof(T)))

#define CustomValues_SizeOfImplicitValues(CV,T) \
        CustomValues_Size(CustomValues_TypeOfImplicitValues(CV),T)
        
#define CustomValues_SizeOfExplicitValues(CV,T) \
        CustomValues_Size(CustomValues_TypeOfExplicitValues(CV),T)
        
#define CustomValues_SizeOfConstantValues(CV,T) \
        CustomValues_Size(CustomValues_TypeOfConstantValues(CV),T)





//----------------------------------------------------------------------
// Math Operations as non-member functions
//----------------------------------------------------------------------
#define CLASSDEF  class... C
#define CLASSLIST C...
#include "CustomValues_Non-MemberOperations.in"

#undef CLASSDEF
#undef CLASSLIST

#define CLASSDEF  class IM,class EX,class CO,class... OT
#define CLASSLIST IM,EX,CO,OT...
#include "CustomValues_Non-MemberOperations.in"

#undef CLASSDEF
#undef CLASSLIST


#endif

