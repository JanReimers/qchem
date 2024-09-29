// File: SphericalGaussianBS.C  Spherical gaussian basis set.



#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "BasisSetImplementation/NumericalIE.H"
#include "BasisSet/IntegralDataBase.H"
#include "oml/vector.h"
#include <iostream>
#include <cassert>

template <class T> inline void FillPower(Vector<T>& arr,T start, T stop)
{
  double del=(std::log(stop/start))/(double)(arr.size()-1);
  typename Vector<T>::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i=T(start*std::exp(n*del));
}

//#######################################################################
//
//  Concrete  gaussian basis set.
//
SphericalGaussianBS::SphericalGaussianBS()
    :  BasisSetImplementation        ()
    , TBasisSetImplementation<double>()
{};

SphericalGaussianBS::SphericalGaussianBS(IntegralDataBase<double>* theDB,
        index_t size,
        double minexp,
        double maxexp,
        int theL,
        Mesh* theMesh)
    : BasisSetImplementation(new SphericalSymmetryQN(theL))
    , TBasisSetImplementation<double>(theDB)
{
    Vector<double> exp(size);
    FillPower(exp,minexp,maxexp);

    Vector<double>::const_iterator b(exp.begin());
    for(; b!=exp.end(); b++) BasisSetImplementation::Insert(new SphericalGaussianBF(*b, theL) );

    if (theMesh)
    {
        assert(theL==0);
        TBasisSetImplementation<double>::Insert(new NumericalIE<double>(theMesh));
    }
    else
        TBasisSetImplementation<double>::Insert(new SphericalGaussianIE );

};

std::ostream&  SphericalGaussianBS::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        BasisSetImplementation::Write(os);
        TBasisSetImplementation<double>::Write(os);
    }
    else
    {
        const SphericalGaussianBF* sgbf=0;
        os << "Spherical Gaussian L=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        // No UT coverage
        auto bs=begin();
        sgbf=dynamic_cast<const SphericalGaussianBF*>(*bs);
        assert(sgbf);
        os << sgbf << "... ";
        for (unsigned int i=0; i<GetNumFunctions()-1; i++) bs++;
        sgbf=dynamic_cast<const SphericalGaussianBF*>(*bs);
        assert(sgbf);
        os << sgbf << "}" << std::endl;
    }
    return os;
}

std::istream&  SphericalGaussianBS::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    BasisSetImplementation::Read(is);
    TBasisSetImplementation<double>::Read(is);
    return is;
}

BasisSet* SphericalGaussianBS::Clone() const
{
    return new SphericalGaussianBS(*this);
}

BasisSet* SphericalGaussianBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return Clone();
}


