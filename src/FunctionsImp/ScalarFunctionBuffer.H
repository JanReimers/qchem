// File: ScalarFunctionBuffer.H  Provides buffering implementation for one mesh.
#ifndef _ScalarFunctionBuffer_H_
#define _ScalarFunctionBuffer_H_



#include "Functions/ScalarFunction.H"
#include "Misc/UniqueID/UniqueID.H"
#include "oml/imp/stream.h"
#include "oml/vector.h"
#include "oml/vector.h"

//-------------------------------------------------------------------------
//
//  This class implements a buff which stores this ScalarFunction
//  Evaluated at all points in the supplied Mesh.  Derived classes
//  are encouraged to provide a more efficient version ScalarMeshEval.
//
template <class T> class ScalarFunctionBuffer
    : public virtual ScalarFunction<T>
    , public virtual StreamableObject
{
protected:
    typedef ScalarFunction<T> Base;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename Base::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::Vec     Vec;  //Vector of scalars.
    typedef typename Base::Vec3    Vec3;   //3 vector (possibly complex).
    typedef typename Base::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename Base::RVec3   RVec3;  //Real space vector.
    typedef typename Base::RVec    RVec;

public:
    ScalarFunctionBuffer(                            );
    ScalarFunctionBuffer(bool Pickle, bool PickleGrad);

    virtual T       operator()(const RVec3&) const=0;
    virtual Vec     operator()(const Mesh& ) const  ;

    virtual Vec3    Gradient  (const RVec3&) const=0;
    virtual Vec3Vec Gradient  (const Mesh& ) const  ;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void  MakeBufferDirty() const
    {
        LastMeshID=0;
        LastGradMeshID=0;
    }

private:
    mutable UniqueID::IDtype  LastMeshID;
    bool              PickleBuffer;
    mutable Vec               Buffer;

    mutable UniqueID::IDtype  LastGradMeshID;
    bool              PickleGradBuffer;
    mutable Vec3Vec           GradBuffer;
};


#endif // _ScalarFunctionBuffer_H_
