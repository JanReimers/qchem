// File: BasisSet/Lattice_3D/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
module qchem.BasisSet.Lattice_3D.BasisSet;
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
import qchem.Symmetry.Factory;                // BlochFactory (the Gamma irrep)
import qchem.Types;

namespace BasisSet::Lattice_3D
{

namespace
{
//! A BasisSet<dcmplx> holding the plane-wave Bloch block(s); owns the IBS list (deleted with the basis).
class PW_BasisSet : public ::BasisSet::BasisSetImp<dcmplx>
{
public:
    PW_BasisSet(const ::Lattice_3D& lat, double Ecut,
                const LocalPotential* loc, const SeparablePotential* nl)
    {
        // Single-k: one Bloch block at Gamma (kIndex 0).  Phase 2 loops the BZ k-mesh here -- the basis
        // ctor is intended to be the one place that enumerates k (the framework's irrep = the Bloch k).
        ReciprocalLattice recip=lat.Reciprocal();
        sym_t gamma=Symmetry::BlochFactory(lat.GetLimits(), ivec3_t(0,0,0));
        auto* pw=new PlaneWave_IBS(recip, gamma, Ecut);
        if (loc) pw->SetPseudopotential(loc,nl);   // pseudo lives on the basis until the Phase-4 wall
        Insert(pw);                                 // BasisSetImp takes ownership
    }
};
} //anon

Complex_BS* Factory(Type type, const ::Lattice_3D& lat, double Ecut,
                    const LocalPotential* loc, const SeparablePotential* nl)
{
    assert(type==Type::PW);
    return new PW_BasisSet(lat, Ecut, loc, nl);
}

} //namespace
