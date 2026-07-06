// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Band_FT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (PW_Pseudo)
import qchem.Fitting.FunctionFitter;        // Make{Density,Scalar}Fitter + Projected{Density,Scalar}_G (both PW fitters)
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Ewald;                           // NuclearRepulsion (Ewald lattice sum for the crystal)
import qchem.Blaze;                            // blazem::zeroH (PW_IonIon's zero matrix)

namespace qchem::Hamiltonian
{

PW_Pseudo::PW_Pseudo(const st_t& st, const Pseudopotential::LocalPotential* loc,
                         const Pseudopotential::SeparablePotential* nl)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
    , itsSep(nl)
{
    assert(st->GetNumAtoms()>0);
    assert(loc && "PW_Pseudo: the term owns the local pseudopotential model (must be non-null)");
}

// Assemble the external matrix from the MODELS the term owns: hand the basis the abstract local +
// optional KB nonlocal models and let it assemble <i|V_ext|j>.  The dynamic_cast is the sanctioned
// abstract->abstract move (cobs_t = Orbital_1E_IBS<dcmplx> ACROSS to the Integrals_Pseudo capability); only a
// basis that supports reciprocal-space PP assembly answers it.  V = V_loc + V_NL.
chmat_t PW_Pseudo::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
    assert(pw && "PW_Pseudo requires an Integrals_Pseudo<dcmplx> (e.g. plane-wave) basis");
    chmat_t V=pw->MakeLocalPotential(&*theStructure, *itsLocal);
    if (itsSep) V += pw->MakeSeparablePotential(&*theStructure, *itsSep);
    return V;
}

void PW_Pseudo::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Een stays the band expectation over the (G!=0) external matrix (== the prototype's electron-ion
    // energy).  The dropped-G=0 alignment alpha is a separate constant in te.Ealign -- kept in the total
    // energy but NOT the matrix, and out of the band-structure cross-check.  It is E_alpha = (N/Omega)
    // Sum_a alpha_a with alpha_a = the model's finite G->0 limit (FormFactorG0 = integral[V_loc^a+Z/r]).
    // The alignment is a PERIODIC neutralising-background artifact: it exists only when the Structure is
    // periodic (!isFinite(), Omega finite); a finite/molecular Structure has no G=0 background, so no
    // alignment term.  The term owns the model (alpha) and reads Omega straight off the geometry -- no basis
    // (the old basis-side PseudoG0Energy was a scalar PP formula leaking into the neutral basis interface).
    te.Een=cd->DM_Contract(this);                                          // integral rho V_ext (G!=0)
    // A finite/molecular structure has no G=0 neutralising background, so NO alignment (even though its atoms
    // DO have form factors -- the physics decision lives here, via isFinite(), not in the geometry).  For a
    // periodic cell the structure folds in 1/Omega: SumFormFactors returns (1/Omega) Sum_a alpha_a.
    if (!theStructure->isFinite())                                         // periodic only (Omega finite)
        te.Ealign = cd->GetTotalCharge() *
                    theStructure->SumFormFactors([this](int Z){return itsLocal->FormFactorG0(Z);});
}

std::ostream& PW_Pseudo::Write(std::ostream& os) const
{
    return os << "    PW external (pseudo)potential with " << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//----------------------------------------------------------------------------------- Kinetic
chmat_t PW_Kinetic::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    chmat_t T=bs->MakeKinetic();   // <p^2> block (no 1/2)
    T*=0.5;                        // T = 1/2 <p^2>
    return T;
}

void PW_Kinetic::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);   // <T> = integral rho (1/2 p^2)
}

std::ostream& PW_Kinetic::Write(std::ostream& os) const
{
    return os << "    PW kinetic energy 1/2<p^2>." << std::endl;
}

//----------------------------------------------------------------------------------- Ion-Ion (Ewald)
PW_IonIon::PW_IonIon(const st_t& st, std::function<double(int)> zionOf)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsZionOf(std::move(zionOf))
{
    assert(st->GetNumAtoms()>0);
    assert(itsZionOf && "PW_IonIon: a Z->Zion charge map is required");
}

// Constant energy term: no Hamiltonian-matrix contribution (the ion-ion energy is independent of the
// electronic state), so the matrix is zero.
chmat_t PW_IonIon::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    return blazem::zeroH<dcmplx>(bs->GetNumFunctions());
}

void PW_IonIon::GetEnergy(EnergyBreakdown& te, const cDM_CD*) const
{
    te.Enn=NuclearRepulsion(*theStructure, itsZionOf);   // Ewald with the PP's ion charges (Zion), not itsZ
}

std::ostream& PW_IonIon::Write(std::ostream& os) const
{
    return os << "    PW ion-ion (Ewald) over " << theStructure->GetNumAtoms() << " ions per cell." << std::endl;
}

//----------------------------------------------------------------------------------- Hartree
// Built once (in Ham_PW_DFT::BuildTerms) from the basis's density-fit basis -- exactly as FittedVee is
// built with its CD fit basis -- so the fitter is created ONCE, not per SCF cycle.
PW_Hartree::PW_Hartree(fbs_t fb)
    : itsFitter(Fitting::MakeDensityFitter(fb))   // the ortho (G-space) density fitter, through the factory
{}
PW_Hartree::~PW_Hartree() = default;   // itsFitter's abstract type is complete here

chmat_t PW_Hartree::CalcMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    newCD(cd);   // dirty the Irrep cache if cd is new (the cross-iteration freshness mechanism)
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "PW_Hartree requires a FourierDensity (periodic) charge density");
    // Reuse the pre-built ortho fitter (mirrors FittedVee): on the orthonormal {G} basis the projection IS
    // the fit, so DoFit just RECEIVES the density's rho-tilde (= Sum_k w_k rho_k), and Repulsion delegates
    // the FFT-free G-space Poisson solve to the basis (Band_FT_IBS).
    itsFitter->DoFit(Fitting::ProjectedDensity_G(fd->GetFourierDensity()));
    return itsFitter->Repulsion(bs);
}

void PW_Hartree::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    te.Eee=0.5*cd->DM_Contract(this,cd);   // E_H = 1/2 integral rho V_H[rho]
}

std::ostream& PW_Hartree::Write(std::ostream& os) const
{
    return os << "    PW Hartree potential ro(r_2)/r_12 (G-space)." << std::endl;
}

//----------------------------------------------------------------------------------- XC
// Built once (in Ham_PW_DFT::BuildTerms) with its Vxc fit basis -- the overlap-metric sibling of PW_Hartree.
PW_XC::PW_XC(const xc_t& xc, fbs_t fb)
    : itsXc(xc)
    , itsScalarFitter(Fitting::MakeScalarFitter(fb))   // the ortho (G-space) scalar fitter, through the factory
{}
PW_XC::~PW_XC() = default;   // itsScalarFitter's abstract type is complete here

// XC through the pre-built ortho scalar fitter, mirroring FittedVxc: build v_xc on the grid (rho(r) via
// inverse FFT, then the functional pointwise), forward-FFT it to its V-tilde (the projection = the fit on
// the orthonormal {G}), DoFit that, then ask for the <i|v_xc|j> assembly (delegated to the basis).  No
// O(Npts*n^2) pointwise density sampling.
chmat_t PW_XC::CalcMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    newCD(cd);
    auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(pw && "PW_XC requires a Band_FT_IBS (plane-wave) basis");
    assert(fd && "PW_XC requires a FourierDensity (periodic) charge density");
    itsBasis=pw;                               // captured for GetEnergy (which has no basis parameter)
    rvec_t rho=pw->RhoOnGrid(fd->GetFourierDensity());
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=itsXc->GetVxc(rho[q]);   // V(r) = v_xc(rho(r))
    itsScalarFitter->DoFit(Fitting::ProjectedScalar_G(pw->ForwardGrid(vxc)));   // V-tilde IS the fit
    return itsScalarFitter->Overlap(bs);                                        // <i|v_xc|j> (no kernel)
}

void PW_XC::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    assert(itsBasis && "PW_XC::GetEnergy called before GetMatrix");
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd);
    rvec_t rho=itsBasis->RhoOnGrid(fd->GetFourierDensity());
    rvec_t exc(rho.size());
    for (size_t q=0;q<rho.size();q++) {double ro=rho[q]; exc[q]=itsXc->GetEpsXc(ro)*ro;}
    te.Exc += itsBasis->Integral(exc);     // E_xc = integral eps_xc(rho) rho  (fresh with this cd)
}

std::ostream& PW_XC::Write(std::ostream& os) const
{
    return os << "    PW exchange-correlation potential v_xc(rho(r))." << std::endl;
}

} //namespace
