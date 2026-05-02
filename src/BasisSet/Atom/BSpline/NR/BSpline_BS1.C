// File: BasisSet/Atom/BSpline/NR/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.BSpline.NR.BS1;

import BasisSet.Atom.BSpline.NR.IBS_Evaluator;
import BasisSet.Atom.BSpline.NR.IBS_Evaluator_r;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;
import qchem.BasisSet.Atom.IBS1;
import qchem.BasisSet.Internal.Common1;
import qchem.BasisSet.Atom.IE1;
import qchem.BasisSet.DB_Cache1;
import qchem.BasisSet1;
import qchem.Orbital_1E_IBS1;

export namespace AtomBS
{
namespace BSpline
{
template <size_t K> class Orbital_IBS1
    : public virtual ::Orbital_IBS1<double>
    , public BSpline_IBS<K>
    , public AtomBS::IrrepBasisSet1
    , public AtomIE_Overlap1<double>
    , public AtomIE_Kinetic1<double>
    , public AtomIE_Nuclear1<double>
{
public:
    Orbital_IBS1(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t&);
    using AtomIE_Overlap1<double>::MakeOverlap;
    // virtual smat_t<double> MakeOverlap() const {return AtomIE_Overlap1<double>::MakeOverlap();}
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet1*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet1*,const Cluster*) const;

    virtual size_t  size           () const {return BSpline_IBS<K>::size();}

};
template <size_t K> class Orbital_IBS_r1
    : public virtual ::Orbital_IBS1<double>
    , public BSpline_r_IBS<K>
    , public AtomBS::IrrepBasisSet1
    , public AtomIE_Overlap1<double>
    , public AtomIE_Kinetic1<double>
    , public AtomIE_Nuclear1<double>
{
public:
    Orbital_IBS_r1(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t&);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet1*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet1*,const Cluster*) const;

    virtual size_t  size           () const {return BSpline_r_IBS<K>::size();}

};

template <size_t K> class Fit_IBS1 
: public BSpline_IBS<K>
, public AtomBS::IrrepBasisSet1
// , public AtomBS::Fit_IBS
{
public:
    Fit_IBS1(size_t N, double rmin, double rmax, size_t L);
};

    // Full basis set.

template <size_t K> class BasisSet1 
    : public BSpline_BS<K> 
    , public ::BS_Common1
    // , public AtomIE_BS_HF<double> //HF support
{
public:
    BasisSet1(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
    using BSpline_BS<K>::BuildCache;
private:
    void Insert(Orbital_IBS1<K>*);
};

template <size_t K> class BasisSet_r1
    : public BSpline_r_BS<K> 
    , public ::BS_Common1
    // , public AtomIE_BS_HF<double> //HF support
{
public:
    BasisSet_r1(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
    using BSpline_r_BS<K>::BuildCache;
private:
    void Insert(Orbital_IBS_r1<K>*);
};

}} //namespace AtomBSl::BSpline

