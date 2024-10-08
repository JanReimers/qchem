// File: SphericalGaussianBS.C  Spherical gaussian basis set.



#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <iostream>
#include <cassert>

namespace SphericalGaussian
{
  
//#######################################################################
//
//  Concrete  gaussian basis set.
//
IrrepBasisSet::IrrepBasisSet()
    :  BasisSetImplementation        ()
    , TBasisSetImplementation<double>()
{};

//
//  We need three constructors type here.  They all need DB, size, exponents,L
//    1) For HF orbitals also need LinearAlgebraParams for secular eq. solving
//    2) For DFT orbitals need 1+LinearAlgebraParams for overlap inversion.
//    3) For DFT Vxc, and ro fitting we need defulat + mesh.
//

IrrepBasisSet::IrrepBasisSet(
        const LinearAlgebraParams& lap,
        IntegralDataBase<double>* theDB,
        size_t size,
        double minexp,
        double maxexp,
        size_t L)
    : BasisSetImplementation(new SphericalSymmetryQN(L))
    , TBasisSetImplementation<double>(lap,theDB)
    , IrrepIEClient(size)
{
    IrrepIEClient::Init(minexp,maxexp,L);
    TBasisSetImplementation<double>::Insert(new IntegralEngine());  
    size_t i=1;
    for (auto e:es) 
        BasisSetImplementation::Insert(new BasisFunction(e,L,ns(i++))); //ns from SphericalGaussianIEClient

};

std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
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

std::istream&  IrrepBasisSet::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    BasisSetImplementation::Read(is);
    TBasisSetImplementation<double>::Read(is);
    return is;
}

::IrrepBasisSet* IrrepBasisSet::Clone() const
{
    return new IrrepBasisSet(*this);
}

::IrrepBasisSet* IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return Clone();
}


} //namespace
