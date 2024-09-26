// File: VectorFunctionBuffer.H  Provides buffering for one mesh.
#ifndef _VectorFunctionBuffer_H_
#define _VectorFunctionBuffer_H_



#include "Functions/VectorFunction.H"
#include "Misc/UniqueID/UniqueID.H"
#include "oml/imp/stream.h"
#include "oml/matrix.h"

//-------------------------------------------------------------------------
//
//  This class implements a buff which stores this Function
//  Evaluated at all points in the supplied Mesh.
//  For matricies the first index is the vector size, the second index is the mesh size.
//
template <class T> class VectorFunctionBuffer
    : public virtual VectorFunction<T>
    , public virtual StreamableObject
{
protected:
    typedef VectorFunction<T> Base;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename Base::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::Vec     Vec;  //Vector of scalars.
    typedef typename Base::Vec3    Vec3;   //3 vector (possibly complex).
    typedef typename Base::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename Base::RVec3   RVec3;  //Real space vector.
    typedef typename Base::RVec    RVec;
    typedef typename Base::Vec3Mat Vec3Mat;//matrix of 3 space vectors.

public:
    VectorFunctionBuffer(                            );
    VectorFunctionBuffer(bool Pickle, bool PickleGrad);

    virtual Vec     operator()(const RVec3&) const=0;
    virtual Mat     operator()(const Mesh& ) const  ;

    virtual Vec3Vec Gradient  (const RVec3&) const=0;
    virtual Vec3Mat Gradient  (const Mesh& ) const  ;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void  MakeBufferDirty() const
    {
        LastMeshID=0;
        LastGradMeshID=0;
    }

private:
    using Base::GetVectorSize;

    mutable UniqueID::IDtype  LastMeshID;
    bool              PickleBuffer;
    mutable Mat               Buffer;

    mutable UniqueID::IDtype  LastGradMeshID;
    bool              PickleGradBuffer;
    mutable Vec3Mat           GradBuffer;
};

#endif // _VectorFunctionBuffer_H_
