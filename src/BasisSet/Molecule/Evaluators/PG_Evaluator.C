// File: BasisSet/Molecule/Evaluators/PG_Evaluator.C
//
// First concrete molecular evaluator (Goal B): wraps a PolarizedGaussian PGData -- the flattened
// (radial x polarization) component set -- plus the cluster (for the multi-centre Nuclear integral).
// The inline 1E kernels are thin wrappers over GaussianRF::Integrate, the very primitives that
// M_PG_Oracle validates to machine precision, with the per-component normalization n_i*n_j folded in
// exactly as PolarizedGaussian's MakeIntegrals does.
//
// Satisfies is1E_Evaluator (checked below).  DFT (3-centre), HF (4-centre), and the VectorFunction
// grid-eval (operator()(r)) come in later increments.
module;
#include <cassert>
#include <string>
#include <ostream>
export module qchem.BasisSet.Molecule.Evaluators.PG;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;      // PGData
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;  // GaussianRF, IType
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;// Polarization
import qchem.BasisSet.Internal.IntegralEnums;                          // qchem::IType3C
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators
{
namespace PG = ::BasisSet::Molecule::PolarizedGaussian;

class PG_Evaluator : public virtual Evaluator
{
public:
    // The PGData (and cluster) must outlive the evaluator -- this is a lightweight view, like the atom
    // evaluators that reference their host IBS's data.
    PG_Evaluator(const PG::PGData& pg, const Cluster* cl) : itsPG(pg), itsCluster(cl) {}

    // --- cold-path Evaluator interface ---
    virtual size_t        size    () const {return itsPG.size();}
    virtual rvec_t        Norm    () const {return itsPG.ns;}
    virtual std::string   RadialID() const {return itsPG.RadialID();}
    virtual std::string   Name    () const {return "PolarizedGaussian";}
    virtual std::ostream& Write   (std::ostream& os) const {return os << "PG_Evaluator[" << size() << "]";}

    // --- 1E inline kernels (hot loops) ---
    // Grad2 is the FULL Cartesian \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block: NO 1/2 (that is
    // the Hamiltonian's), NO centrifugal term (atom-only).  Nuclear is the multi-centre attraction.
    double Norm   (size_t i)          const {return itsPG.ns[i];}
    double Overlap(size_t i,size_t j) const {return Integrate(PG::Overlap2C, i, j);}
    double Grad2  (size_t i,size_t j) const {return Integrate(PG::Grad2,     i, j);}
    double Nuclear(size_t i,size_t j) const {return Integrate(PG::Nuclear,   i, j, itsCluster);}

    // --- 3-centre (DFT) and 4-centre (HF) kernels ---------------------------------------------------
    // General multi-evaluator elements: each (evaluator, index) pair names one basis component, so the
    // same kernel serves both Coulomb/Exchange (the caller just maps its loop indices into the slots).
    // `this` is the A slot.  Normalizations of every slot are folded in, matching MakeOverlap3C /
    // MakeDirect / MakeExchange.  (A member may read another same-type evaluator's private data.)

    // <ab|c> 3-centre integral of the given IType3C (M&D Integrate3C, called on the C radial).
    double ThreeC(qchem::IType3C t, size_t iA,
                  const PG_Evaluator& eB, size_t iB,
                  const PG_Evaluator& eC, size_t iC) const
    {
        return eC.itsPG.radials[iC]->Integrate(t, *itsPG.radials[iA], *eB.itsPG.radials[iB],
                                               itsPG.pols[iA], eB.itsPG.pols[iB], eC.itsPG.pols[iC])
             * itsPG.ns[iA] * eB.itsPG.ns[iB] * eC.itsPG.ns[iC];
    }

    // (ab|cd) 4-centre electron-repulsion integral (M&D Integrate4C, called on the D radial).
    double FourC(size_t iA,
                 const PG_Evaluator& eB, size_t iB,
                 const PG_Evaluator& eC, size_t iC,
                 const PG_Evaluator& eD, size_t iD) const
    {
        return eD.itsPG.radials[iD]->Integrate(*itsPG.radials[iA], *eB.itsPG.radials[iB],
                                               *eC.itsPG.radials[iC],
                                               itsPG.pols[iA], eB.itsPG.pols[iB],
                                               eC.itsPG.pols[iC], eD.itsPG.pols[iD])
             * itsPG.ns[iA] * eB.itsPG.ns[iB] * eC.itsPG.ns[iC] * eD.itsPG.ns[iD];
    }

private:
    double Integrate(PG::IType t, size_t i, size_t j, const Cluster* cl=0) const
    {
        return itsPG.radials[i]->Integrate(t, *itsPG.radials[j], itsPG.pols[i], itsPG.pols[j], cl)
             * itsPG.ns[i] * itsPG.ns[j];
    }

    const PG::PGData& itsPG;
    const Cluster*    itsCluster;
};

static_assert(is1E_Evaluator<PG_Evaluator>, "PG_Evaluator must satisfy is1E_Evaluator");

} //namespace
