// File: Spin.H  A 3 state variable indicating the spin polarization.
#ifndef _Spin_H_
#define _Spin_H_

#include "oml/imp/index_t.h"

enum class Spin {Down,None,Up};

inline bool   IsPolarized  (Spin s) {return !(s==Spin::None);}
inline int    GetDegeneracy(Spin s) {return IsPolarized(s) ? 1 : 2;}
inline size_t SequenceIndex(Spin s) {return static_cast<size_t>(s);}

#endif //_Spin_H_
