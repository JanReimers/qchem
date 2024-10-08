// File: SphericalGaussianBS.C  Spherical gaussian basis set.



#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE1.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
#include "BasisSetImplementation/NumericalIEImp.H"
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
        size_t L,
        Mesh* theMesh)
    : BasisSetImplementation(new SphericalSymmetryQN(L))
    , TBasisSetImplementation<double>(lap,theDB)
    , SphericalGaussianIEClient(size)
{
    SphericalGaussianIEClient::Init(minexp,maxexp,L);
//    Vector<double> exp(size);
//    FillPower(exp,minexp,maxexp);
    for (auto e:es) BasisSetImplementation::Insert(new SphericalGaussianBF(e,L));
    TBasisSetImplementation<double>::Insert(new SphericalGaussianIE1(L,es));  
    if (theMesh)
    {
        assert(L==0); //Why???
        TBasisSetImplementation<double>::Insert(new NumericalIEImp<double>(theMesh));
    }
          

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


