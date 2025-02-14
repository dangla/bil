#ifndef CUSTOMVALUES_H
#define CUSTOMVALUES_H

#include <type_traits>
#include "Message.h"


/* Primary template */
template<typename T,template<typename> class... C>
struct alignas(T) CustomValues_t: C<T>... {
  using Value_type = T;
  
  #include "CustomValues_MemberOperations.in"
};


/* Partial template specialization */
template<typename T,template<typename> class IM,template<typename> class EX,template<typename> class CO,template<typename> class... OT>
struct alignas(T) CustomValues_t<T,IM,EX,CO,OT...>: IM<T>,EX<T>,CO<T>,OT<T>... {
  using ImplicitValues_type = IM<T>;
  using ExplicitValues_type = EX<T>;
  using ConstantValues_type = CO<T>;
  using Value_type = T;
  
  #include "CustomValues_MemberOperations.in"
};



/* Below CV (=TCV<T>) stands for a class. */
#include "Utils.h"

#define CustomValues_Index(...) \
        Utils_CAT_NARG(CustomValues_Index,__VA_ARGS__)(__VA_ARGS__)
        
#define CustomValues_Index2(TCV,V) \
        CustomValues_Index3(TCV<char>,V,char)
        
#define CustomValues_Index3(CV,V,T)  ((int) ((T*)&((CV*)(0))->V - (T*)0))


#define CustomValues_TypeOfImplicitValues(CV)   typename CV::ImplicitValues_type
#define CustomValues_TypeOfExplicitValues(CV)   typename CV::ExplicitValues_type
#define CustomValues_TypeOfConstantValues(CV)   typename CV::ConstantValues_type
#define CustomValues_TypeOfValue(CV)            typename CV::Value_type

#define CustomValues_IsValueType(CV,T) \
        std::is_same_v<CustomValues_TypeOfValue(CV),T>



#define CustomValues_NbOfImplicitValues(CV) \
        ((int) (sizeof(CustomValues_TypeOfImplicitValues(CV))/sizeof(CustomValues_TypeOfValue(CV))))
        
#define CustomValues_NbOfExplicitValues(CV) \
        ((int) (sizeof(CustomValues_TypeOfExplicitValues(CV))/sizeof(CustomValues_TypeOfValue(CV))))
        
#define CustomValues_NbOfConstantValues(CV) \
        ((int) (sizeof(CustomValues_TypeOfConstantValues(CV))/sizeof(CustomValues_TypeOfValue(CV))))





//----------------------------------------------------------------------
// Math Operations as non-member functions
//----------------------------------------------------------------------
#define CLASSDEF  typename U,template<typename> class... A
#define CLASSLIST U,A...
#include "CustomValues_Non-MemberOperations.in"
#undef CLASSDEF
#undef CLASSLIST

#define CLASSDEF  typename U,template<typename> class A,template<typename> class B,template<typename> class C,template<typename> class... D
#define CLASSLIST U,A,B,C,D...
#include "CustomValues_Non-MemberOperations.in"
#undef CLASSDEF
#undef CLASSLIST


#endif

