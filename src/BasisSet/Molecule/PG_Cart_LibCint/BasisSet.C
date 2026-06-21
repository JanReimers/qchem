// File: BasisSet/Molecule/PG_Cart_LibCint/BasisSet.C
//
// Production orbital-basis tree for the Cartesian PG basis integrated by libcint -- the IBS parallel to
// PG_Cart, but built on the matrix-delivery PG_Cart_LibCint evaluator instead of PG_Cart_MnD.  The data
// (radials, blocks, PGData component layout) is read EXACTLY as in PG_Cart -- only the integral engine
// differs -- so a full SCF run here is the end-to-end cross-check of libcint against M&D (HF total energy,
// which is basis-ordering invariant; the per-element match is already in tests/M_LibCint.C).
//
// HF only (1E + 4-centre).  The basis is a single C1 IrrepBasisSet (no SALC), enough for the energy
// cross-check.  A DFT 3-centre fit (a libcint EFit_IBS) and a native-spherical sibling are later increments.
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PG_Cart_LibCint;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Reader;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_LibCint;   // NR_Evaluator: the IBS IS-A evaluator base
import qchem.BasisSet.Molecule.IBS;                          // Molecule::Orbital_{1E,HF}_IBS<E> mixins

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::PG_Cart_LibCint
{
namespace Cart = ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;
namespace LC   = ::BasisSet::Molecule::Evaluators::PG_Cart_LibCint;

class IrrepBasisSet
        : public virtual Real_IBS,
          public IrrepBasisSetImp<double>,
          public LC::NR_Evaluator                  // IS-A libcint evaluator, which IS-A PGData
    {
    public:
        typedef std::vector<std::unique_ptr<Cart::Block>> bv_t;

        IrrepBasisSet(Reader*, const Cluster*, bool spherical=false);
        IrrepBasisSet(const rvec_t& exponents, size_t L, const Cluster*, bool spherical=false);
        virtual ~IrrepBasisSet();

        // The number of functions is the evaluator's delivered count (Cartesian (l+1)(l+2)/2 or spherical
        // 2l+1 per shell), NOT the raw Cartesian PGData count -- in spherical mode they differ.
        virtual size_t  GetNumFunctions() const {return NR_Evaluator::size();}
        virtual size_t  size() const {return NR_Evaluator::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        // Geometry-aware cache key with a DISTINCT prefix per engine+mode: the global Jac/Kab cache keys on
        // BasisSetID, and PGData's is " PG {.." == the M&D PG_Cart IBS, so without a distinct prefix the
        // libcint run would serve / be served M&D's cached matrices in the same process (as PG1-vs-PG).
        virtual std::string BasisSetID() const
        { return (itsSpherical ? "LibCintSph " : "LibCint ") + PGData::BasisSetID(); }
        virtual std::string Name     () const {return "Pol. Gaussian (libcint) ";}
        virtual std::ostream &Write(std::ostream&) const;

    private:
        bv_t itsBlocks;
        bool itsSpherical=false;
    };

// 1E + 4-centre HF inherited from the Molecule-generic, evaluator-templated mixins instantiated with the
// libcint evaluator (the IBS IS-A that evaluator).  Because the evaluator is isM_*, the mixins FORWARD to
// its assembled matrices instead of running the per-element loop.
class Orbital_IBS
    : public Molecule::Orbital_1E_IBS<LC::NR_Evaluator>
    , public Molecule::Orbital_HF_IBS<LC::NR_Evaluator>
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader*, const Cluster*, bool spherical=false);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Cluster*, bool spherical=false);
};

class BasisSet
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
{
public:
    BasisSet() {};
    BasisSet(Reader*, const Cluster*, bool spherical=false);
    virtual void Insert(bs_t* bs);
};

} //namespace
