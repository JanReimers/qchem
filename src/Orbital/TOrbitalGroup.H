// File: TOrbitalGroup.H  Group of orbitals.
#ifndef _TOrbitalGroup_H_
#define _TOrbitalGroup_H_



#include "Orbital/OrbitalGroup.H"
#include "Functions/VectorFunction.H"
#include "Misc/ptr_vector1.h"

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbitalGroup
    : public virtual OrbitalGroup
    , public virtual VectorFunction<T>
{
public:
    typedef dynamic_cast_iterator<Orbital*,TOrbital<T>* > const_iterator;
    virtual index_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    const_iterator beginT() const {return const_iterator(begin());};

    static  TOrbitalGroup* Factory(std::istream&)        ;
};

#endif //_TOrbitalGroup_H_
