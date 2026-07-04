// File:: Hamiltonian/Internal/Imp/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <memory>
#include <string>
#include <utility>
#include <vector>
module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.PWTerms;        // PW_Kinetic/External/Hartree/XC (the plane-wave KS terms)
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Types;
import qchem.Structure;
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH + GTH_PP + HGH_*/MultiSpecies_* (re-exported)
import qchem.PeriodicTable;                       // PeriodicTable::GetZ(symbol) -> atomic number (the composite key)

namespace qchem::Hamiltonian
{

Ham_1E::Ham_1E(const st_t& st)
{
    InsertStandardTerms(st);
}

Ham_HF_U::Ham_HF_U(const st_t& st) 
{
    InsertStandardTerms(st);
    Add(new Vee);
    Add(new Vxc(-0.5));
}


Ham_DFT_U::Ham_DFT_U(const st_t& st,double alpha_ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFT_U(st,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const st_t& st,ExFunctional* ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    InsertStandardTerms(st);
       
    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
    Add(new FittedVxc(XFitBasis, XcFunct));
}

// Dirac exchange + VWN5 correlation: the parameter-free in-house LSDA (delegates to the generic ctor).
Ham_DFTcorr_U::Ham_DFTcorr_U(const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFTcorr_U(st, new SlaterExchange(2.0/3.0), new VWN_Correlation(), mp, bs)
{}

// Generic separate-terms LSDA: exchange via FittedVxc (3/4 virial energy, exact for exchange) + correlation
// via FittedVcorr (E_c = integral eps_c rho -- needs the functional's GetEpsXc), sharing ONE Vxc fit basis
// so the 3-centre integrals are computed once.  Used by both the in-house (Slater+VWN) and libxc paths, so
// the correct correlation energy is shared -- no path lumps X+C into a single 3/4-virial term.
Ham_DFTcorr_U::Ham_DFTcorr_U(const st_t& st, ExFunctional* exchange, ExFunctional* correlation,
                             const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    InsertStandardTerms(st);

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp)); // ONE Vxc fit basis, shared X and C
    FittedVxc::ex_t exch(exchange);
    Add(new FittedVxc  (XFitBasis, exch));
    FittedVxc::ex_t corr(correlation);
    Add(new FittedVcorr(XFitBasis, corr));
}

// Spin-native polarized LSDA: mirror Ham_DFTcorr_U but with the polarized exchange (FittedVxcPol, Dirac)
// and the polarized correlation (FittedVcorrPol, spin-native VWN5) terms, sharing one Vxc fit basis.  The
// unpolarized Ham_DFTcorr_U is the zeta=0 collapse of this.
Ham_DFTcorr_P::Ham_DFTcorr_P(const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    InsertStandardTerms(st);

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxcPol::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp)); // ONE Vxc fit basis, shared X and C
    FittedVxcPol::ex_t exch(new SlaterExchange(2.0/3.0, Spin(Spin::Up))); // Dirac exchange (alpha = 2/3), polarized
    Add(new FittedVxcPol  (XFitBasis, exch));
    FittedVcorrPol::corr_t corr(new VWN_Correlation());                  // spin-native VWN5 correlation
    Add(new FittedVcorrPol(XFitBasis, corr));
}

// PSEUDOPOTENTIAL LSDA: like Ham_DFTcorr_U but with the bare nuclear attraction (Ven) replaced by the
// mesh-quadratured local pseudopotential V_loc(r) + the KB-separable non-local projectors, and NO ion-ion.
// Kinetic + PP_Local [+ PP_NonLocal] + Hartree + Dirac exchange + VWN5.
Ham_PP_U::Ham_PP_U(const st_t& st, std::shared_ptr<const Pseudopotential::LocalPotential_R> vloc,
                   std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep,
                   const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    Add(new Kinetic);
    Add(new PP_Local(st, std::move(vloc), mp));      // pseudized replacement for Ven; NO Vnn
    if (sep) Add(new PP_NonLocal(st, std::move(sep), mp));   // KB separable projectors (null => local-only)

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis, st->GetNumElectrons()));

    FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp)); // ONE Vxc fit basis, shared X and C
    FittedVxc::ex_t exch(new SlaterExchange(2.0/3.0));            // Dirac exchange (alpha = 2/3)
    Add(new FittedVxc  (XFitBasis, exch));
    FittedVxc::ex_t corr(new VWN_Correlation());                 // VWN5 correlation
    Add(new FittedVcorr(XFitBasis, corr));
}

Ham_PP_U::Ham_PP_U(const st_t& st, const std::string& element, int q, const qcMesh::MeshParams& mp,
                   const rbs_t* bs)
    : Ham_PP_U(st,
               std::make_shared<const Pseudopotential::HGH_LocalPotential>(Pseudopotential::GetGTH(element,"LDA",q).local),
               std::make_shared<const Pseudopotential::HGH_SeparablePotential>(Pseudopotential::GetGTH(element,"LDA",q).nonlocal),
               mp, bs)
{}

// Plane-wave LDA Kohn-Sham: the five G-space framework terms.  Exchange and correlation are SEPARATE
// PW_XC terms (Dirac + VWN5), mirroring Ham_DFTcorr_U, so the correlation energy is the correct
// E_c = integral eps_c rho.  No fit basis / mesh: the plane-wave basis owns the integration, and the
// pseudopotential is carried by the basis (the external term just supplies the structure factor).
void Ham_PW_DFT::BuildTerms(const st_t& st, const Pseudopotential::LocalPotential* loc,
                            const Pseudopotential::SeparablePotential* nl)
{
    Add(new PW_Kinetic);
    Add(new PW_Pseudo(st, loc, nl));                           // electron-ion (incl. G=0 alignment)
    Add(new PW_Hartree);
    Add(new PW_XC(std::make_shared<SlaterExchange>(2.0/3.0)));    // Dirac exchange (alpha = 2/3)
    Add(new PW_XC(std::make_shared<VWN_Correlation>()));          // VWN5 correlation
    Add(new PW_IonIon(st, loc->ZionFn()));                       // ion-ion Ewald: Zion from the PP, not itsZ
}

// Explicit-models ctor: the caller owns the models (itsOwnedLocal/Sep stay null).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const Pseudopotential::LocalPotential* loc,
                       const Pseudopotential::SeparablePotential* nl)
{
    BuildTerms(st, loc, nl);
}

// Single-species convenience ctor: the 1-species case of the multi-species build.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const std::string& element,
                       const std::string& functional, int valence)
{
    BuildFromGTH(st, {{element, valence}}, functional);
}

// Multi-species convenience ctor.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, std::initializer_list<std::pair<std::string,int>> species,
                       const std::string& functional)
{
    BuildFromGTH(st, std::vector<std::pair<std::string,int>>(species), functional);
}

// Look up each (element, valence) from the GTH database and build + OWN a per-Z router model (one
// MultiSpecies_Local + one MultiSpecies_Separable, keyed by atomic number so the assembly's per-atom
// FormFactor(a->itsZ,...) dispatches to the right species).  The owned models outlive the terms (members,
// destroyed after the cHamiltonian base that holds them), so each term's &loc/&nl stays valid for the run.
void Ham_PW_DFT::BuildFromGTH(const st_t& st, const std::vector<std::pair<std::string,int>>& species,
                              const std::string& functional)
{
    auto loc=std::make_shared<Pseudopotential::MultiSpecies_LocalPotential>();
    auto sep=std::make_shared<Pseudopotential::MultiSpecies_SeparablePotential>();
    for (const auto& [element, valence] : species)
    {
        int Z=thePeriodicTable().GetZ(element);          // atomic number = the atoms' itsZ key
        Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH(element, functional, valence);
        loc->Add(Z, std::make_shared<Pseudopotential::HGH_LocalPotential>(pp.local));
        sep->Add(Z, std::make_shared<Pseudopotential::HGH_SeparablePotential>(pp.nonlocal));
    }
    itsOwnedLocal=loc;
    itsOwnedSep  =sep;
    BuildTerms(st, loc.get(), sep.get());
}

Ham_HF_P::Ham_HF_P(const st_t& st)
{
    InsertStandardTerms(st);
    Add(new Vee);
    Add(new VxcPol);
}


Ham_DFT_P::Ham_DFT_P(const st_t& st,double alpha_ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFT_P(st,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const st_t& st,ExFunctional* ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    InsertStandardTerms(st);
    FittedVee::fbs_t CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
    Add(new FittedVxcPol(XFitBasis, XcFunct));
    
}


Ham_DHF_1E::Ham_DHF_1E(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Ven(st));
    //Add(new Vnn(st));
}

Ham_DHF_U::Ham_DHF_U(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(st));
    Add(new Ven(st));
    Add(new Vee());
    Add(new Vxc(-0.5));
}
Ham_DHF_P::Ham_DHF_P(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(st));
    Add(new Ven(st));
    Add(new Vee());
    Add(new VxcPol());
}

} //namespace
