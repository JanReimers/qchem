// File: BasisSet/Lattice_3D/Imp/GPW_IBS.C  GPW_IBS implementation (ctors + identity).
module;
#include <iostream>
#include <memory>
#include <string>

module qchem.BasisSet.Lattice_3D.GPW_IBS;
import qchem.Symmetry.Factory;              // BlochFactory (the convenience ctor + the k=0 fit-basis irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (pry k out of the abstract Bloch irrep)

namespace qchem::BasisSet::Lattice_3D
{

GPW_IBS::GPW_IBS(const UnitCell& cell, const sym_t& irrep,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, double Rcut)
    : BasisSet::IrrepBasisSetImp<dcmplx>(irrep)
    , GPW_Evaluator(std::move(mol), cell, densityEcut, Symmetry::Lattice_3D::Getk(irrep), Rcut) // Bloch irrep IS k
{}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the primary constructor.
GPW_IBS::GPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, double Rcut)
    : GPW_IBS(cell, Symmetry::BlochFactory(N,kIndex), std::move(mol), densityEcut, Rcut)
{}

// The DFT fit-basis factory: a plane-wave fit basis over GPW's OWN density grid (k=0 Gamma; the density is
// cell-periodic).  Both CD and Vxc share the one density grid this increment (relCutoff refinement deferred).
BasisSet::cFIT_CD_ABS* GPW_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams&) const
{
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}
BasisSet::cFIT_SF_ABS* GPW_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const
{
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

std::string GPW_IBS::BasisSetID() const
{
    return Name()+GPW_Evaluator::IDFragment();   // Name + "|mol=..|k=..|nR=.."
}

std::ostream& GPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " periodic Gaussians, " << GetSymmetry();
}

} //namespace
