// File: ScalarFunctionBuffer.C  Buffer a scalar function evaluated over a mesh.



#include "FunctionsImp/ScalarFunctionBuffer.H"
#include "Mesh/Mesh.H"
#include "oml/imp/binio.h"
#include <iostream>
#include "oml/io3d.h"


template <class T> ScalarFunctionBuffer<T>::ScalarFunctionBuffer()
    : LastMeshID      (0    )
    , PickleBuffer    (false)
    , Buffer          (     )
    , LastGradMeshID  (0    )
    , PickleGradBuffer(false)
    , GradBuffer      (     )
{};

template <class T> ScalarFunctionBuffer<T>::ScalarFunctionBuffer(bool Pickle, bool PickleGrad)
    : LastMeshID      (0         )
    , PickleBuffer    (Pickle    )
    , Buffer          (          )
    , LastGradMeshID  (0         )
    , PickleGradBuffer(PickleGrad)
    , GradBuffer      (          )
{};

template <class T> typename ScalarFunction<T>::Vec ScalarFunctionBuffer<T>::operator() (const Mesh& m) const
{
    if (LastMeshID!=m.GetID())
    {
//    cout << "ScalarFunctionBuffer(" << typeid(*this).name() << "):: Allocating " << m.GetNumPoints() << " vector, "
//      << m.GetNumPoints()*8 << " bytes" << std::endl;
        LastMeshID=m.GetID();
        Buffer.SetLimits(VecLimits(m.GetNumPoints()));
        Fill(Buffer,0.0);
        ScalarFunction<T>::Eval(m, Buffer);
    }

    return Buffer;
}

template <class T> typename ScalarFunction<T>::Vec3Vec ScalarFunctionBuffer<T>::Gradient(const Mesh& m) const
{
    if (LastGradMeshID!=m.GetID())
    {
//    cout << "ScalarFunctionGradBuffer:: Allocating " << m.GetNumPoints() << " vector, "
//      << m.GetNumPoints()*8*3 << " bytes" << std::endl;
        LastGradMeshID=m.GetID();
        GradBuffer.SetLimits(VecLimits(m.GetNumPoints()));
        Fill(GradBuffer,Vec3(0,0,0));
        ScalarFunction<T>::EvalGrad(m, GradBuffer);
    }

    return GradBuffer;
}

template <class T> std::ostream& ScalarFunctionBuffer<T>::Write(std::ostream& os) const
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

template <class T> std::istream& ScalarFunctionBuffer<T>::Read (std::istream& is)
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

template class ScalarFunctionBuffer<double>;
