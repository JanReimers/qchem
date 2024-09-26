// File: VectorFunctionBuffer.C  Buffer a real space vector function evaluated over a mesh.



#include "FunctionsImp/VectorFunctionBuffer.H"
#include "Mesh/Mesh.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <complex>

template <class T> VectorFunctionBuffer<T>::VectorFunctionBuffer()
    : LastMeshID      (0    )
    , PickleBuffer    (false)
    , Buffer          (     )
    , LastGradMeshID  (0    )
    , PickleGradBuffer(false)
    , GradBuffer      (     )
{};

template <class T> VectorFunctionBuffer<T>::VectorFunctionBuffer(bool Pickle, bool PickleGrad)
    : LastMeshID      (0         )
    , PickleBuffer    (Pickle    )
    , Buffer          (          )
    , LastGradMeshID  (0         )
    , PickleGradBuffer(PickleGrad)
    , GradBuffer      (          )
{};

template <class T> typename VectorFunction<T>::Mat VectorFunctionBuffer<T>::operator() (const Mesh& m) const
{
    if (LastMeshID!=m.GetID())
    {
    std::cout << "VectorFunctionBuffer:: Allocating " << GetVectorSize() << "X" << m.GetNumPoints() << " matrix, "
      << GetVectorSize()*m.GetNumPoints()*8 << " bytes" << std::endl;
        LastMeshID=m.GetID();
        Buffer.SetLimits(MatLimits(GetVectorSize(),m.GetNumPoints()));
        Fill(Buffer,T(0.0));
        VectorFunction<T>::Eval(m,Buffer);
    }

    return Buffer;
}

template <class T> typename VectorFunction<T>::Vec3Mat VectorFunctionBuffer<T>::Gradient(const Mesh& m) const
{
    if (LastGradMeshID!=m.GetID())
    {
//    cout << "VectorFunctionGradBuffer:: Allocating " << GetVectorSize() << "X" << m.GetNumPoints() << " matrix, "
//      << GetVectorSize()*m.GetNumPoints()*8*3 << " bytes" << std::endl;
        LastGradMeshID=m.GetID();
        GradBuffer.SetLimits(MatLimits(GetVectorSize(),m.GetNumPoints()));
        Fill(GradBuffer,Vec3(0,0,0));
        VectorFunction<T>::EvalGrad(m,GradBuffer);
    }

    return GradBuffer;
}



template <class T> std::ostream& VectorFunctionBuffer<T>::Write(std::ostream& os) const
{
    if (PickleBuffer)
    {
        if (Ascii()) os << LastMeshID << " " << Buffer << " ";
        if (Binary())
        {
            BinaryWrite(LastMeshID,os);
            os << Buffer;
        }
    }
    else
    {
        if (Ascii()) os << -1 << " ";
        if (Binary())BinaryWrite(-1,os);
    }

    if (PickleGradBuffer)
    {
        if (Ascii()) os << LastGradMeshID << " " << GradBuffer << " ";
        if (Binary())
        {
            BinaryWrite(LastGradMeshID,os);
            os << GradBuffer;
        }
    }
    else
    {
        if (Ascii()) os << -1 << " ";
        if (Binary())BinaryWrite(-1,os);
    }
    return os;
}

template <class T> std::istream& VectorFunctionBuffer<T>::Read (std::istream& is)
{
    int tempID;
    PickleBuffer=false;
    if (Ascii()) is >> tempID >> std::ws;
    if (Binary()) BinaryRead(tempID,is);
    if (tempID >= 0)
    {
        LastMeshID=tempID;
        is >> Buffer;
        PickleBuffer=true;
    }
    PickleGradBuffer=false;
    if (Ascii()) is >> tempID >> std::ws;
    if (Binary()) BinaryRead(tempID,is);
    if (tempID >= 0)
    {
        LastGradMeshID=tempID;
        is >> GradBuffer;
        PickleGradBuffer=true;
    }
    return is;
}

template class VectorFunctionBuffer<double>;
template class VectorFunctionBuffer<std::complex<double> >;
