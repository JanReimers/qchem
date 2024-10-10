// File: SphericalGaussianBS.C  Spherical gaussian basis set.



#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <iostream>
#include <cassert>

//#######################################################################
//
//  Concrete  gaussian basis set.
//
SphericalGaussianBS::SphericalGaussianBS()
    :  BasisSetImplementation        ()
    , TBasisSetImplementation<double>()
{};

//
//  We need three constructors type here.  They all need DB, size, exponents,L
//    1) For HF orbitals also need LinearAlgebraParams for secular eq. solving
//    2) For DFT orbitals need 1+LinearAlgebraParams for overlap inversion.
//    3) For DFT Vxc, and ro fitting we need defulat + mesh.
//

SphericalGaussianBS::SphericalGaussianBS(
        const LinearAlgebraParams& lap,
        IntegralDataBase<double>* theDB,
        size_t size,
        double minexp,
        double maxexp,
        size_t L)
    : BasisSetImplementation(new SphericalSymmetryQN(L))
    , TBasisSetImplementation<double>(lap,theDB)
    , SphericalGaussianIEClient(size)
{
    SphericalGaussianIEClient::Init(minexp,maxexp,L);
    TBasisSetImplementation<double>::Insert(new SphericalGaussianIE1());  
    size_t i=1;
    for (auto e:es) 
        BasisSetImplementation::Insert(new SphericalGaussianBF(e,L,ns(i++))); //ns from SphericalGaussianIEClient

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
        os << "Spherical Gaussian L=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
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

IrrepBasisSet* SphericalGaussianBS::Clone() const
{
    return new SphericalGaussianBS(*this);
}

IrrepBasisSet* SphericalGaussianBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return Clone();
}


