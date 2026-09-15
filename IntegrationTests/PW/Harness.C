// File: IntegrationTests/PW/Harness.C  THE PLANE-WAVE SCF HARNESS: the two drivers the PW grid runs on --
// the standalone PROTOTYPE loop (RunSCF / RunSCF_kpoints: the pre-framework reference, kept because
// PW_Si.Γ_eqPrototype pins the framework against it) and the framework driver (RunFrameworkGamma over
// SolidSCFIterator).  PW has NO facade yet (SolidCalculation is built over a Gaussian orbital basis), which
// is why these exist; when the facade grows its basis-family axis they retire like the GPW drivers did.
// Its own library (qcPW_Harness) so IntegrationTests/PW/*.C compile one copy.  (2026-09-15, TE phase 5.)
module;
#include <map>
#include <set>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cstdio>
export module qchem.Tests.PW_Harness;
export import qchem.Tests.PW_Fields;
export import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
export import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
export import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
export import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
export import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
export import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
export import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
export import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
export import qchem.Blaze;          // mat_t<dcmplx>
export import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // Ven_PP_Short/Long, Vee_Hartree (dcmplx Hamiltonian terms)
export import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
export import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Hamiltonian.Internal.IonIon;           // IonIon<double> (ion-ion term: pair sum / Ewald via isFinite)
import qchem.Hamiltonian.Internal.Kinetic;          // Kinetic<dcmplx> (the shared kinetic term)
export import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
export import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // PeriodicIrrepCD<dcmplx> (the periodic leaf -- 3c-2b)
export import qchem.ChargeDensity.SeedCD;                  // SeedCD + PolarizedSeedCD (the spin-SAD staggering gate)
export import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
export import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
export import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
export import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
export import qchem.Symmetry.Irrep;                        // Irrep
export import qchem.LASolver;                              // complex Hermitian eigensolver
export import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
export import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
export import qchem.SCFIterator;                           // cSCFIterator (the real framework SCF driver)
export import qchem.SCFParams;                             // SCFParams
export import qchem.WaveFunction;                          // cWaveFunction (read view of the converged state)
export import qchem.ElectronConfiguration;                 // ElectronConfiguration base
export import qchem.ElectronConfiguration.Crystal;         // Crystal_EC (single-k Bloch configuration)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;      // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

export namespace qchem::tests::pw
{
using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

// Vee_Hartree/Vxc_Quadrature now take their fit basis (from the basis's own factory) at construction, like
// FittedVee/FittedVxc.  These low-level term tests build one straight from the plane-wave basis at hand.
qchem::Hamiltonian::Vee_Hartree* NewPWHartree(const PlaneWave_IBS& pw)
{
    // Pure V_H[rho]: the Hartree term takes ONLY its fit basis now -- the long-range core-charge fold that
    // used to need a structure + local model here is its own term, Ven_PP_Long (doc/GPWPlan.md 0e-PP).
    return new qchem::Hamiltonian::Vee_Hartree(
        qchem::Hamiltonian::Vee_Hartree::fbs_t(pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{})));
}
// The XC term takes a QUADRATURE, and MakeDensitySampler reads off the strategy the fit basis supports --
// here a plane-wave (raster) basis, so the pair/collocation one.
qchem::Hamiltonian::Vxc_Quadrature* NewPWXC(const PlaneWave_IBS& pw, const qchem::Hamiltonian::Vxc_Quadrature::xc_t& xc)
{
    return new qchem::Hamiltonian::Vxc_Quadrature(xc, qchem::ChargeDensity::MakeDensitySampler(
        qchem::ChargeDensity::fitbasis_t(pw.CreateVxcFitBasisSet(nullptr, qcMesh::MeshParams{}))), SpinGroup::UnPolarized);
}

//! Evaluate v_xc(rho(r)) on the uniform grid (returns the field + its coordinates \a rf), and by
//! reference E_xc = integral eps_xc rho and the XC double-count VxcDc = integral rho v_xc.
//! \a vxcOf / \a epsxcOf are rho->v_xc and rho->eps_xc (sum exchange+correlation by the caller).
std::vector<double> XcGridField(const RhoG& rho, const ivec3_t& Ng, double Omega,
                                const std::function<double(double)>& vxcOf,
                                const std::function<double(double)>& epsxcOf,
                                std::vector<rvec3_t>& rf, double& Exc, double& VxcDc)
{
    rf=UniformGrid(Ng);
    size_t Npts=rf.size();
    std::vector<double> vxc(Npts);
    Exc=0.0; VxcDc=0.0;
    for (size_t q=0;q<Npts;q++)
    {
        double rr=RhoOfR(rho,rf[q]);
        if (rr<0.0) rr=0.0;                 // clamp tiny numerical negatives before the functional
        vxc[q]=vxcOf(rr);
        Exc   += epsxcOf(rr)*rr;
        VxcDc += vxc[q]*rr;
    }
    Exc   *= Omega/double(Npts);
    VxcDc *= Omega/double(Npts);
    return vxc;
}

//! Build Vxc~(dm) over a single basis's difference set (single-k convenience).
RhoG BuildXcVtilde(const PlaneWave_IBS& pw, const RhoG& rho, const ivec3_t& Ng, double Omega,
                   const std::function<double(double)>& vxcOf,
                   const std::function<double(double)>& epsxcOf, double& Exc, double& VxcDc)
{
    std::vector<rvec3_t> rf;
    std::vector<double>  vxc=XcGridField(rho,Ng,Omega,vxcOf,epsxcOf,rf,Exc,VxcDc);
    RhoG vtil;                              // forward-transform once per distinct dm in the diff set
    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=pw.GetGIndex(i)-pw.GetGIndex(j);
            if (vtil.find(dm)==vtil.end()) vtil[dm]=ForwardDFT(vxc,rf,dm);
        }
    return vtil;
}

// --- the SCF loop (prototype) -----------------------------------------------------------------
// Assemble H(k) = 1/2 <p^2> + V_ext + V_H[rho] + V_xc[rho], diagonalise, occupy the lowest bands,
// rebuild rho, linearly mix, iterate to self-consistency.  This is the hand-rolled driver that the
// templated Hamiltonian/SCFIterator<dcmplx> will eventually replace (V_H/V_xc -> tDynamic_HT<dcmplx>).

struct SCFResult
{
    RhoG   rho;                                   //!< converged self-consistent density rho~(dm)
    double Ekin=0, EH=0, Exc=0, Eext=0;           //!< energy components at the fixed point
    double Eband=0, VxcDc=0;                      //!< band sum, and XC double-count integral rho v_xc
    double Etot_band=0;                           //!< band route:  Eband - E_H + E_xc - integral rho v_xc
    double Etot_direct=0;                         //!< direct route: T + integral rho V_ext + E_H + E_xc
    double gap=0;                                 //!< CBM-VBM over the sampled k-points (multi-k only)
    int    iters=0;
    bool   converged=false;
};

// \a Vext is the FULL static external block (local PP + KB nonlocal): density-independent, assembled
// once by the caller.  E_ext is taken as the band expectation Sum_n f_n <psi|Vext|psi> -- valid for the
// nonlocal projector part too (where there is no real-space integral rho V_ext).
SCFResult RunSCF(const PlaneWave_IBS& pw, const UnitCell& B, double Omega,
                 const chmat_t& Vext, int Nelec, const ivec3_t& Ng,
                 const std::function<double(double)>& vxcOf, const std::function<double(double)>& epsxcOf,
                 double alpha=0.5, double tol=1e-9, int maxiter=400)
{
    size_t n=pw.GetNumFunctions();
    int    Nocc=Nelec/2;
    chmat_t K=pw.MakeKinetic();                   // diagonal <p^2>, constant across iterations
    chmat_t S=pw.MakeOverlap();                   // identity

    SCFResult R;
    R.rho[ivec3_t(0,0,0)] = dcmplx(double(Nelec)/Omega, 0.0);   // uniform initial guess

    mat_t<dcmplx> U; rvec_t e; rvec_t f; RhoG rhoNew;
    for (int it=0; it<maxiter; it++)
    {
        R.iters=it+1;
        auto   Vh=HartreeVtilde(R.rho,B);
        double Exc,VxcDc;
        RhoG   vxcmap=BuildXcVtilde(pw,R.rho,Ng,Omega,vxcOf,epsxcOf,Exc,VxcDc);
        auto   Vloc=[&](const ivec3_t& dm)->dcmplx { return Vh(dm)+RhoAt(vxcmap,dm); };  // Hartree+XC

        hmat_t<dcmplx> H(n);                       // 1/2 <p^2> (diag) + Vext(i,j) + V_H+V_xc (Gi-Gj)
        for (size_t i=0;i<n;i++)
            for (size_t j=i;j<n;j++)
                H(i,j) = (i==j ? dcmplx(0.5*std::real(dcmplx(K(i,i)))) : dcmplx(0.0))
                       + dcmplx(Vext(i,j))
                       + Vloc(pw.GetGIndex(i)-pw.GetGIndex(j));

        LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
        las->SetBasisOverlap(S);
        auto sol=las->Solve(H);
        delete las;
        U=std::get<0>(sol); e=std::get<1>(sol);

        std::vector<size_t> idx(e.size());         // occupy the Nocc lowest-energy bands (2 each)
        for (size_t i=0;i<idx.size();i++) idx[i]=i;
        std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return e[a]<e[b];});
        f=rvec_t(e.size(),0.0);
        for (int k=0;k<Nocc;k++) f[idx[k]]=2.0;

        rhoNew=BuildDensity(pw,U,f,Omega);
        double drho=0.0;                           // max |rho_new - rho_in| over all components
        for (const auto& kv : rhoNew) drho=std::max(drho, std::abs(kv.second-RhoAt(R.rho,kv.first)));

        RhoG mixed;                                // linear mix
        for (const auto& kv : rhoNew) mixed[kv.first]=(1-alpha)*RhoAt(R.rho,kv.first)+alpha*kv.second;
        R.rho=mixed;

        if (drho<tol) { R.converged=true; break; }
    }

    // Final energies at the self-consistent density (rhoNew = output of the converged solve).
    R.EH=HartreeEnergy(rhoNew,B,Omega);
    double dummyExc,dummyDc;
    BuildXcVtilde(pw,rhoNew,Ng,Omega,vxcOf,epsxcOf,dummyExc,dummyDc);
    R.Exc=dummyExc; R.VxcDc=dummyDc;
    for (size_t b=0;b<e.size();b++)                          // E_ext = Sum_n f_n <psi|Vext|psi>
        if (f[b]!=0.0)
        {
            dcmplx ev(0.0);
            for (size_t i=0;i<n;i++)
                for (size_t j=0;j<n;j++) ev += std::conj(dcmplx(U(i,b)))*dcmplx(Vext(i,j))*dcmplx(U(j,b));
            R.Eext += f[b]*std::real(ev);
        }
    for (size_t b=0;b<e.size();b++) R.Eband += f[b]*e[b];
    for (size_t b=0;b<e.size();b++)
        if (f[b]!=0.0)
            for (size_t i=0;i<n;i++) R.Ekin += f[b]*0.5*std::real(dcmplx(K(i,i)))*std::norm(dcmplx(U(i,b)));
    R.Etot_band   = R.Eband - R.EH + R.Exc - R.VxcDc;
    R.Etot_direct = R.Ekin + R.Eext + R.EH + R.Exc;
    return R;
}

// --- multi-k (Brillouin-zone) SCF ---------------------------------------------------------------
// The BZ sum Sum_k w_k is the irrep loop the framework will eventually own (k IS a Bloch irrep).  Each
// k has its OWN plane-wave basis (cutoff set {G:1/2|k+G|^2<Ecut} depends on k) and its OWN external
// block (the KB projectors beta~(|k+G|) are k-dependent); both are static, so built once.  The density
// rho~(dm) (k-independent, lattice-periodic) accumulates Sum_k w_k Sum_n f c c* over all k.  For an
// insulator (Si: 8 valence e-) exactly Nelec/2 bands are filled at EVERY k -> trivial occupation.
SCFResult RunSCF_kpoints(const ReciprocalLattice& recip, const UnitCell& B, double Omega,
                         const Structure& st, const HGH_LocalPotential& loc, const HGH_SeparablePotential& nl,
                         const ivec3_t& Nmp, double Ecut, int Nelec, const ivec3_t& Ng,
                         const std::function<double(double)>& vxcOf, const std::function<double(double)>& epsxcOf,
                         double alpha, double tol, int maxiter)
{
    std::vector<std::unique_ptr<PlaneWave_IBS>> pw;     // per-k static data (basis, external block, kinetic)
    std::vector<chmat_t> Vext, K;
    for (int kx=0;kx<Nmp.x;kx++)
        for (int ky=0;ky<Nmp.y;ky++)
            for (int kz=0;kz<Nmp.z;kz++)
            {
                auto p=std::make_unique<PlaneWave_IBS>(recip, Nmp, ivec3_t(kx,ky,kz), Ecut);
                Vext.push_back(p->MakeSpeciesFieldMatrix(&st, loc, qchem::BasisSet::FieldRange::Full)+p->MakeProjectorMatrix(&st,nl));
                K.push_back(p->MakeKinetic());
                pw.push_back(std::move(p));
            }
    int    Nk=int(pw.size());
    double wk=1.0/Nk;
    int    Nocc=Nelec/2;

    std::set<ivec3_t,IVecLess> diff;                    // union of all k difference sets (for the XC transform)
    for (int k=0;k<Nk;k++)
        for (size_t i=0;i<pw[k]->GetNumFunctions();i++)
            for (size_t j=0;j<pw[k]->GetNumFunctions();j++)
                diff.insert(pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j));

    SCFResult R;
    R.rho[ivec3_t(0,0,0)] = dcmplx(double(Nelec)/Omega, 0.0);
    std::vector<mat_t<dcmplx>> Uk(Nk); std::vector<rvec_t> ek(Nk), fk(Nk); RhoG rhoNew;

    for (int it=0; it<maxiter; it++)
    {
        R.iters=it+1;
        auto Vh=HartreeVtilde(R.rho,B);
        std::vector<rvec3_t> rf; double Exc,VxcDc;
        std::vector<double> vxcField=XcGridField(R.rho,Ng,Omega,vxcOf,epsxcOf,rf,Exc,VxcDc);
        RhoG vxcmap; for (const ivec3_t& dm : diff) vxcmap[dm]=ForwardDFT(vxcField,rf,dm);
        auto Vloc=[&](const ivec3_t& dm)->dcmplx { return Vh(dm)+RhoAt(vxcmap,dm); };

        rhoNew.clear();
        for (int k=0;k<Nk;k++)
        {
            size_t n=pw[k]->GetNumFunctions();
            hmat_t<dcmplx> H(n), S=pw[k]->MakeOverlap();
            for (size_t i=0;i<n;i++)
                for (size_t j=i;j<n;j++)
                    H(i,j) = (i==j ? dcmplx(0.5*std::real(dcmplx(K[k](i,i)))) : dcmplx(0.0))
                           + dcmplx(Vext[k](i,j)) + Vloc(pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j));

            LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
            las->SetBasisOverlap(S);
            auto sol=las->Solve(H);
            delete las;
            Uk[k]=std::get<0>(sol); ek[k]=std::get<1>(sol);

            std::vector<size_t> idx(ek[k].size());                 // occupy the Nocc lowest bands
            for (size_t i=0;i<idx.size();i++) idx[i]=i;
            std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return ek[k][a]<ek[k][b];});
            fk[k]=rvec_t(ek[k].size(),0.0);
            for (int o=0;o<Nocc;o++) fk[k][idx[o]]=2.0;

            for (size_t b=0;b<fk[k].size();b++)                    // accumulate weighted density
                if (fk[k][b]!=0.0)
                    for (size_t i=0;i<n;i++)
                        for (size_t j=0;j<n;j++)
                            rhoNew[pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j)]
                                += (wk*fk[k][b]/Omega)*dcmplx(Uk[k](i,b))*std::conj(dcmplx(Uk[k](j,b)));
        }

        double drho=0.0;
        for (const auto& kv : rhoNew) drho=std::max(drho, std::abs(kv.second-RhoAt(R.rho,kv.first)));
        RhoG mixed;
        for (const auto& kv : rhoNew) mixed[kv.first]=(1-alpha)*RhoAt(R.rho,kv.first)+alpha*kv.second;
        R.rho=mixed;
        if (drho<tol) { R.converged=true; break; }
    }

    R.EH=HartreeEnergy(rhoNew,B,Omega);                            // energies at the converged density
    {
        std::vector<rvec3_t> rfe; double exc,dc;
        XcGridField(rhoNew,Ng,Omega,vxcOf,epsxcOf,rfe,exc,dc);
        R.Exc=exc; R.VxcDc=dc;
    }
    double vbm=-1e30, cbm=1e30;
    for (int k=0;k<Nk;k++)
    {
        for (size_t b=0;b<fk[k].size();b++)
        {
            R.Eband += wk*fk[k][b]*ek[k][b];
            if (fk[k][b]!=0.0)
            {
                for (size_t i=0;i<pw[k]->GetNumFunctions();i++)
                    R.Ekin += wk*fk[k][b]*0.5*std::real(dcmplx(K[k](i,i)))*std::norm(dcmplx(Uk[k](i,b)));
                dcmplx ev(0.0);
                size_t n=pw[k]->GetNumFunctions();
                for (size_t i=0;i<n;i++)
                    for (size_t j=0;j<n;j++)
                        ev += std::conj(dcmplx(Uk[k](i,b)))*dcmplx(Vext[k](i,j))*dcmplx(Uk[k](j,b));
                R.Eext += wk*fk[k][b]*std::real(ev);
            }
        }
        std::vector<double> es;                                    // sorted bands at this k -> gap
        for (size_t b=0;b<ek[k].size();b++) es.push_back(ek[k][b]);
        std::sort(es.begin(),es.end());
        vbm=std::max(vbm, es[Nocc-1]);
        cbm=std::min(cbm, es[Nocc]);
    }
    R.gap=cbm-vbm;
    R.Etot_band   = R.Eband - R.EH + R.Exc - R.VxcDc;
    R.Etot_direct = R.Ekin + R.Eext + R.EH + R.Exc;
    return R;
}

// the lattice, seed a UNIFORM density, run the full SCFIterator (complex DIIS), and report.  Takes
// ownership of \a ham (the SCFIterator deletes it).  Mirrors the Si framework test body.
struct FwResult { bool converged; double charge; qchem::EnergyBreakdown E; size_t iters; };
// Takes a Ham FACTORY (not a pre-built Ham) so the Hamiltonian is built WITH the basis -- the fit-basis
// seam: Ham_PW_DFT now needs the composite basis to create the Hartree density-fit basis (like the
// molecular Ham DFT ctors take bs for FittedVee).
FwResult RunFrameworkGamma(const Lattice_3D& lat, double Ecut, int Nelec,
                           std::function<qchem::Hamiltonian::cHamiltonian*(const BasisSet::Complex_BS*)> mkHam,
                           const char* label,
                           qchem::ChargeDensity::SeedStrategy seed=qchem::ChargeDensity::SeedStrategy::Uniform)
{
    namespace L3=BasisSet::Lattice;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, Ecut));
    qchem::Hamiltonian::cHamiltonian* ham=mkHam(bs.get());   // build the Ham WITH the basis (fit-basis seam)
    std::unique_ptr<qchem::Hamiltonian::cHamiltonian> hamOwner(ham);   // R2.22: the iterator borrows; this scope owns
    Irrep  irr=bs->GetIrreps(Spin::None)[0];
    size_t n  =bs->GetNumFunctions();
    Crystal_EC ec(irr, Nelec);
    using qchem::SCFAccelerators::DIISParams;
    // EMax must exceed the ionic [F,D] error (~1.4-3) or DIIS bails ("En>EMax") and linear mixing
    // oscillates on the strong Madelung field.  Engage DIIS from the start (EMax large).
    auto* acc=new qchem::SCFAccelerators::SCFAcceleratorDIIS(DIISParams{10, 8.0, 1e-10, 1e-9});
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> accOwner(acc);   // R2.22: the iterator borrows; this scope owns
    // \a seed defaults to Uniform rho(r)=N/V; IonicSAD pre-bakes the formal-charge transfer (Na+ + F-).
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc, seed, lat.GetStructure().get());
    SCFParams par;
    par.NMaxIter=120; par.MinΔρ=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=false;
    scf.Iterate(par);
    const qchem::WaveFunction::cWaveFunction* wf=scf.GetWaveFunction();
    auto cd=wf->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    qchem::EnergyBreakdown E=scf.GetEnergy();
    std::cout << "["<<label<<"] nPW="<<n<<" iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy() << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]
              << " Eee="<<E["Eee"]<<" Exc="<<E["Exc"]<<" Enn="<<E["Enn"]<<" E_alphaZ="<<E["E_alphaZ"]<<")" << std::endl;
    return {scf.Converged(), charge, E, scf.GetIterationCount()};
}

} // namespace qchem::tests::pw
