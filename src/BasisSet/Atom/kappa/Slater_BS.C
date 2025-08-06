// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;

export namespace Atom_kappa
{
namespace Slater
{

class IE_Primatives_slkappa : public ::Slater::IE_Primatives 
{
    virtual double Inv_r1   (double ea, double eb,size_t l_total) const;
};

//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r) as:
//
//                                  -e*r^l*exp(-e*r), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)g = {
//                                  ((2l+1)/r-e)*r^lexp(-e*r), kappa>0
//


// Irrep basis set

// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Slater function
// specific.
class Orbital_RKB_IBS
    : public IrrepBasisSet_Common<double>
    , public Orbital_RKB_IBS_Common<double> 
    , public IE_Primatives_slkappa
{
public:
    Orbital_RKB_IBS(const DB_cache<double>* db,const ::IE_Primatives* pie,const Vector<double>& exponents, int kappa);
    virtual size_t size() const {return Orbital_RKB_IBS_Common<double>::size();}
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    friend class ::DiracIntegralTests;

};

template <class T> class Orbital_RKBL_IBS
    : public ::Slater::IrrepBasisSet //Use NR slater basis
    , public Atom::Orbital_RKBL_IBS<T>
{
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const ::IE_Primatives* pie, const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return size();}
};

template <class T> class Orbital_RKBS_IBS
    : public IrrepBasisSet_Common<T> 
    , public Atom::Orbital_RKBS_IBS<T> 
    , public AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    using Orbital_RKBS_IBS_Common<T>::large;
public:
    Orbital_RKBS_IBS(const DB_cache<double>*, const ::IE_Primatives* pie,const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return size();}

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;

    virtual const SMatrix<T>& Overlap() const {return DB_Kinetic<T>::Kinetic();}

private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBS_IBS_Common<T>::kappa;
};

// Full basis set
class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
    , public ::Slater::IE_Primatives
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lMax);
    
};

}} //namespace Atom_kappa::Slater

