// File: BasisSet/Atom/BSpline/NR/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.BSpline.NR.BS;

import BasisSet.Atom.BSpline.NR.IBS_Evaluator;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace BSpline
{
template <size_t K> class Orbital_IBS
    : public BSpline_IBS<K>
    , public Atom::IrrepBasisSet
    , public Atom::Orbital_HF_IBS <double>
    , public Atom::Orbital_IBS    <double>
    , public Atom::Orbital_DFT_IBS<double>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L,  const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;

    virtual size_t  size           () const {return BSpline_IBS<K>::size();}

};

template <size_t K> class Fit_IBS 
: public BSpline_IBS<K>
, public Atom::IrrepBasisSet
, public Atom::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L);
};

    // Full basis set.

template <size_t K> class BasisSet 
    : public BSpline_BS<K> 
    , public ::BS_Common
    , public AtomIE_BS_2E<double> //HF support
{
public:
    BasisSet(size_t N, double rmin, double rmax, size_t Lmax); 
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
    using BSpline_BS<K>::BuildCache;
private:
    void Insert(Orbital_IBS<K>*);
};

}} //namespace Atoml::BSpline

