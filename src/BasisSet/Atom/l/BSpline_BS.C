// File: BasisSet/Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <bspline/Core.h>
#include <iosfwd>
export module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
import BasisSet.Atom.IBS_Evaluator;
import BasisSet.Atom.BSpline_IBS;
import qchem.Basisset.Atom.radial.BSpline.BS_Common;
import qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet;
import qchem.Fit_IBS;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace BSpline
{
template <size_t K> class Orbital_IBS
    : public         ::BSpline::IrrepBasisSet<K>
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public AtomIE_Overlap<double>
    , public AtomIE_Kinetic<double>
    , public ::BSpline::IE_Inv_r1<double,K>
    , public ::BSpline::IE_DFT<double,K>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const ::BSpline::IE_Primatives<K>*,const IBS_Evaluator* eval,size_t N, double rmin, double rmax, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const ::BSpline::IE_Primatives<K>*,const IBS_Evaluator* eval,size_t N, double rmin, double rmax, size_t L,  const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

template <size_t K> class Fit_IBS 
: public ::BSpline::IrrepBasisSet<K>
, public BSpline_IBS<K>
, public Fit_IBS_Common
, public ::BSpline::IE_Fit<K>
, public ::BSpline::IE_Overlap<double,K>
{
public:
    Fit_IBS(const DB_cache<double>* db,const ::BSpline::IE_Primatives<K>* pie,size_t N, double rmin, double rmax, size_t L);
};

    // Full basis set.

template <size_t K> class BasisSet 
    : public ::BSpline::BS_Common<K>
    , public ::BSpline::IE_Primatives_Imp<K>
{
public:
    BasisSet(size_t N, double rmin, double rmax, size_t Lmax); 
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
   
};

}} //namespace Atoml::BSpline

