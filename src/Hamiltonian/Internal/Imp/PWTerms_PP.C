// File: Hamiltonian/Internal/Imp/PWTerms_PP.C  the three external-PP terms: V_loc short, V_loc long, and the KB separable projectors.
//
// One implementation unit of module qchem.Hamiltonian.Internal.PWTerms.  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
module;
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cassert>
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.RunPolicy;   // theRunPolicy().XCFromDM() -- the declared XC-feed deviation (N5)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform: the fit basis's FFT pair (RhoOnGrid, the BALL route)
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)


namespace qchem::Hamiltonian
{


// The dropped-G=0 alignment, PER ELECTRON, evaluated ONCE at construction (R2.16).
//
// E_alpha = N * (1/Omega) Sum_a alpha_a, with alpha_a the model's finite G->0 limit
// (FormFactorG0 = integral[V_loc^a + Z/r]).  It is kept in the total energy but NOT in the matrix, so it
// stays out of the band-structure cross-check.  The alignment is a PERIODIC neutralising-background
// artifact: a finite/molecular Structure has no G=0 background, so its coefficient is exactly ZERO (even
// though its atoms DO have form factors -- the physics decision lives here, not in the geometry).
//
// Both the isFinite() question and the SumFormFactors sum are answered HERE rather than per energy call:
// the structure and the model are fixed at construction, so the answer cannot change during a run.  Each
// GetEnergy then just scales by the current electron count.
namespace
{
double G0AlignmentPerElectron(const Structure& st, const std::function<double(int)>& formFactorG0)
{
    if (st.isFinite()) return 0.0;              // no neutralising background => no alignment, exactly 0
    return st.SumFormFactors(formFactorG0);     // periodic: folds in 1/Omega
}
} // namespace

Ven_PP_Short::Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
{
    assert(st->GetNumAtoms()>0);
    assert(loc && "Ven_PP_Short: the term owns the local pseudopotential model (must be non-null)");
    itsAlphaZ=G0AlignmentPerElectron(*st, [loc](int Z){return loc->FormFactorG0Short(Z);});
}

// Assemble the external matrix from the MODEL the term owns: hand the basis the abstract local model and
// let it assemble <i|V_loc,short|j>.  The dynamic_cast is the sanctioned abstract->abstract move
// (cobs_t = Orbital_1E_IBS<dcmplx> ACROSS to the Integrals_Pseudo capability); only a basis that supports
// reciprocal-space PP assembly answers it.
template <class U> hmat_t<U> Ven_PP_Short::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pw && "Ven_PP_Short requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    // SHORT-range local only.  The LONG (softened-Coulomb) half is Ven_PP_Long and the KB projectors are
    // Ven_PP_NonLocal (the CP2K local-PP split, doc/GPWPlan.md 0e-PP).
    return pw->MakeLocalPotentialShort(&*theStructure, *itsLocal);
}
chmat_t Ven_PP_Short::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Short::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_Short::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Een is the band expectation over the (G!=0) matrix (== the prototype's electron-ion energy).
    te.Een     += cd->DM_Contract(this);                 // integral rho V_loc,short (G!=0)
    te.E_alphaZ+= cd->GetTotalCharge()*itsAlphaZ;        // SHORT G=0 alignment (0 for a finite structure)
}

std::ostream& Ven_PP_Short::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: SHORT-range local PP, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//--------------------------------------------------------------- Ven_PP_NonLocal (the KB projectors)
Ven_PP_NonLocal::Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsSep(nl)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a local-only PP omits this TERM rather than adding an all-zeros one.
    if (!nl)
        throw std::runtime_error("Ven_PP_NonLocal: the KB projector term requires a SeparablePotential "
                                 "model.  A local-only pseudopotential must not add this term.");
}

template <class U> hmat_t<U> Ven_PP_NonLocal::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pw && "Ven_PP_NonLocal requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    return pw->MakeSeparablePotential(&*theStructure, *itsSep);
}
chmat_t Ven_PP_NonLocal::MakeMatrix(const cobs_t* bs, const Spin& s) const
{
    if (std::getenv("GPW_NL_PER_L"))
    {   // I0 diagnostic (doc/SphericalLatticePlan.md): bank the per-l blocks once per irrep block.
        // Complex path only -- the itsByL bank is chmat_t; extend if the diagnostic ever needs real blocks.
        auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
        assert(pw);
        const std::string id=bs->BasisSetID();
        if (itsByLSeen.insert(id).second)
            for (auto& lH : pw->MakeSeparablePotentialByL(&*theStructure, *itsSep))
                itsByL[lH.first].emplace(id, std::move(lH.second));
    }
    return MakeMatrixT<dcmplx>(bs,s);
}
rsmat_t Ven_PP_NonLocal::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_NonLocal::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ion, and short-ranged by construction, so no G=0 alignment of its own.
    const double eNL = cd->DM_Contract(this);            // Tr(D V_NL)
    te.Een   += eNL;
    te.EenNL += eNL;                                     // the diagnostic V_loc/V_NL split (Een keeps the total)
    if (!itsByL.empty())
    {   // GPW_NL_PER_L: the per-channel decomposition of eNL (l=-1 = a basis that answered lumped).
        std::cout << "[NL per-l]";
        double sum=0;
        for (const auto& lB : itsByL)
        {
            const double el=cd->DM_ContractBlocks(lB.second);
            sum+=el;
            std::cout << "  l=" << lB.first << ": " << el;
        }
        std::cout << "  (sum=" << sum << "  Tr(D V_NL)=" << eNL << ")" << std::endl;
    }
}

std::ostream& Ven_PP_NonLocal::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: KB separable nonlocal projectors, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

// Kinetic is now the shared Kinetic<dcmplx> term (qchem.Hamiltonian.Internal.Kinetic).
// Ion-ion (Ewald) is now the shared IonIon<dcmplx> term (qchem.Hamiltonian.Internal.IonIon).

//------------------------------------------------------------------- Ven_PP_Long (the LONG range half)
// The long-range (softened-Coulomb / Gaussian core-charge) local-PP matrix.  DENSITY-INDEPENDENT --
// MakeLocalPotentialLong takes only (structure, model) -- so this is an ordinary static term and rides the
// standard per-Irrep static cache, exactly as Ven_PP_Short does.  It was previously a side-block inside the
// Hartree term (summed into that term's matrix, then subtracted back out of its energy) purely because the
// two are solved through the same G-space Poisson machinery; that is a COMPUTATIONAL kinship, not a
// physical one, and it cost a nullable model, a second block cache, and a "Hartree" term that contributed
// to E_een.  Assembled through the SAME Integrals_Pseudo cross-cast Ven_PP_Short uses.
Ven_PP_Long::Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a run with no local PP omits this TERM rather than constructing an
    // all-zeros one.  (Absence of a capability belongs in the term LIST, not in a per-call branch.)
    if (!loc)
        throw std::runtime_error("Ven_PP_Long: the long-range local-PP term requires a LocalPotential "
                                 "model.  A run without a local pseudopotential must not add this term.");
    itsAlphaZ=G0AlignmentPerElectron(*st, [loc](int Z){return loc->FormFactorG0Long(Z);});
}

template <class U> hmat_t<U> Ven_PP_Long::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pp=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pp && "Ven_PP_Long requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    return pp->MakeLocalPotentialLong(&*theStructure, *itsLocal);
}
chmat_t Ven_PP_Long::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Long::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_Long::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ION, so NO 1/2: the double-counting factor belongs to the electron-electron Hartree alone.
    te.Een     += cd->DM_Contract(this);                 // Tr(D V_long) = E_een,long
    te.E_alphaZ+= cd->GetTotalCharge()*itsAlphaZ;        // LONG G=0 alignment (0 for a finite structure)
}

std::ostream& Ven_PP_Long::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: LONG-range local PP (Gaussian core charge, G-space), "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//----------------------------------------------------------------------------------- Hartree
// Holds its CD fit basis (from Ham_PW_DFT::BuildTerms via the orbital basis's factory, never assuming
// orbital==fit) and hands it to the density's GetRepulsion3C each SCF cycle -- mirrors FittedVee holding its

} //namespace
