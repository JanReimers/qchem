module;
#include <cstddef> 

export module qchem.Symmetry.Spin;

export
{
    enum class Spin {Down,None,Up};
    inline bool   IsPolarized  (Spin s) {return !(s==Spin::None);}
    inline int    GetDegeneracy(Spin s) {return IsPolarized(s) ? 1 : 2;}
    inline size_t SequenceIndex(Spin s) {return static_cast<size_t>(s);}
}