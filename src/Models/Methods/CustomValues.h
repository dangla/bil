#ifndef CUSTOMVALUES_H
#define CUSTOMVALUES_H

template<class I,class E,class C>
struct CustomValues_t: I,E,C {
  using ImplicitValues_type = I;
  using ExplicitValues_type = E;
  using ConstantValues_type = C;
};


/* Convert a pointer V to any type in pointer to type T */
#define CustomValues_Convert(V,T) (((T)*) (V))

#define CustomValues_Index(CV,V,T)  (&(CV)->V - ((T*) (CV)))

#define CustomValues_TypeOfImplicitValues(CV)  typename CV::ImplicitValues_type
#define CustomValues_TypeOfExplicitValues(CV)  typename CV::ExplicitValues_type
#define CustomValues_TypeOfConstantValues(CV)  typename CV::ConstantValues_type

#endif

