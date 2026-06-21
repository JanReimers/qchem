// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Evaluator.C
//
// First concrete molecular evaluator (Goal B): the NON-RELATIVISTIC, CARTESIAN polarized-Gaussian basis
// integrated by McMurchie-Davidson (hence PG_Cart_MnD::NR_Evaluator).  IS-A PolarizedGaussian PGData --
// the flattened (radial x polarization) component set -- and is meant to be a base subobject of the PG
// IrrepBasisSet (the atom evaluators relate to their IBS the same way).  The inline 1E kernels are thin
// wrappers over GaussianRF::Integrate, the very primitives that M_PG_Oracle validates to machine
// precision, with the per-component normalization n_i*n_j folded in as PolarizedGaussian's MakeIntegrals.
//
// Satisfies is1E_Evaluator, isDFT_Evaluator and isHF_Evaluator (checked below).  The VectorFunction
// grid-eval (operator()(r), needed before isOpr can fold into isDFT) comes in a later increment.
//
// Sibling evaluators anticipated: PG_Spherical_MnD, PG_Cart_PRISM, ... each a new Evaluator over (mostly)
// the same shared PGData -- which is why PGData is kept separate, not absorbed here.
module;
#include <cassert>
#include <string>
#include <ostream>
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;      // PGData
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;  // GaussianRF, IType
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;// Polarization
import qchem.BasisSet.Internal.IntegralEnums;                          // qchem::IType3C
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{
namespace PG = ::BasisSet::Molecule::PolarizedGaussian;

class NR_Evaluator : public virtual Evaluator, public PG::PGData
{
public:
    // The evaluator OWNS its data (IS-A PGData) rather than viewing one: it is a base subobject of the PG
    // IrrepBasisSet, exactly as the atom evaluators are bases of their IBS.  The default constructor builds
    // an empty PGData; the host IBS fills it via PGData::Init.  The cluster is NOT evaluator state -- it is
    // passed per call to the Nuclear kernel.
    NR_Evaluator() = default;

    // --- cold-path Evaluator interface ---
    virtual size_t        size    () const {return PGData::size();}
    virtual rvec_t        Norm    () const {return ns;}
    virtual std::string   RadialID() const {return PGData::RadialID();}
    virtual std::string   Name    () const {return "PolarizedGaussian";}
    virtual std::ostream& Write   (std::ostream& os) const {return os << "PG_Cart_MnD::NR_Evaluator[" << size() << "]";}

    // --- 1E inline kernels (hot loops) ---
    // Grad2 is the FULL Cartesian \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block: NO 1/2 (that is
    // the Hamiltonian's), NO centrifugal term (atom-only).  Nuclear is the multi-centre attraction.
    double Norm   (size_t i)          const {return ns[i];}
    double Overlap(size_t i,size_t j) const {return Integrate(PG::Overlap2C, i, j);}
    double Grad2  (size_t i,size_t j) const {return Integrate(PG::Grad2,     i, j);}
    double Nuclear(size_t i,size_t j,const Cluster* cl=0) const {return Integrate(PG::Nuclear, i, j, cl);}

    // --- 3-centre (DFT) and 4-centre (HF) kernels ---------------------------------------------------
    // General multi-evaluator elements: each (evaluator, index) pair names one basis component, so the
    // same kernel serves both Coulomb/Exchange (the caller just maps its loop indices into the slots).
    // `this` is the A slot.  Normalizations of every slot are folded in, matching MakeOverlap3C /
    // MakeDirect / MakeExchange.  (A member may read another same-type evaluator's private data.)

    // <ab|c> 3-centre integral of the given IType3C (M&D Integrate3C, called on the C radial).
    double ThreeC(qchem::IType3C t, size_t iA,
                  const NR_Evaluator& eB, size_t iB,
                  const NR_Evaluator& eC, size_t iC) const
    {
        return eC.radials[iC]->Integrate(t, *radials[iA], *eB.radials[iB],
                                         pols[iA], eB.pols[iB], eC.pols[iC])
             * ns[iA] * eB.ns[iB] * eC.ns[iC];
    }

    // (ab|cd) 4-centre electron-repulsion integral (M&D Integrate4C, called on the D radial).
    double FourC(size_t iA,
                 const NR_Evaluator& eB, size_t iB,
                 const NR_Evaluator& eC, size_t iC,
                 const NR_Evaluator& eD, size_t iD) const
    {
        return eD.radials[iD]->Integrate(*radials[iA], *eB.radials[iB],
                                         *eC.radials[iC],
                                         pols[iA], eB.pols[iB],
                                         eC.pols[iC], eD.pols[iD])
             * ns[iA] * eB.ns[iB] * eC.ns[iC] * eD.ns[iD];
    }

private:
    double Integrate(PG::IType t, size_t i, size_t j, const Cluster* cl=0) const
    {
        return radials[i]->Integrate(t, *radials[j], pols[i], pols[j], cl)
             * ns[i] * ns[j];
    }
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace
