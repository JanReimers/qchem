// File: TOrbitalGroup.H  Group of orbitals.
#ifndef _TOrbitalGroup_H_
#define _TOrbitalGroup_H_



#include "Orbital/OrbitalGroup.H"
#include "Functions/VectorFunction.H"

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
    virtual index_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }

    static  TOrbitalGroup* Factory(std::istream&)        ;
};

#endif //_TOrbitalGroup_H_
