// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <cassert>
#include <iostream>
#include <memory>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.FourierDFT_IBS;         // cast bs UP to the G-space (FFT Hartree/XC) capability
import qchem.BasisSet.DFTPotential_IBS;       // PW_External still uses the real-space external matrix

namespace qchem::Hamiltonian
{

PW_External::PW_External(const st_t& st)
    : cStatic_HT_Imp()
    , theStructure(st)
{
    assert(st->GetNumAtoms()>0);
}

// Ask the basis for the configured external (pseudo)potential matrix.  The dynamic_cast is the
// sanctioned abstract->abstract move (cobs_t = Orbital_1E_IBS<dcmplx> down to the richer abstract
// DFTPotential_IBS<dcmplx>); only a basis that supports G-space DFT assembly answers it.
chmat_t PW_External::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const BasisSet::DFTPotential_IBS<dcmplx>*>(bs);
    assert(pw && "PW_External requires a DFTPotential_IBS (e.g. plane-wave) basis");
    return pw->MakeExternalPotential(&*theStructure);
}

void PW_External::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    te.Een=cd->DM_Contract(this);   // integral rho V_ext (the density contracts our matrix)
}

std::ostream& PW_External::Write(std::ostream& os) const
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

//----------------------------------------------------------------------------------- Hartree
chmat_t PW_Hartree::CalcMatrix(const cobs_t* bs, const Spin&, const cDM_CD* cd) const
{
    newCD(cd);   // dirty the Irrep cache if cd is new (the cross-iteration freshness mechanism)
    // G-space Poisson solve: take the density's reciprocal-space coefficients rho-tilde (composite =
    // Sum_k w_k rho_k) and assemble V_H(dm)=4pi rho-tilde/G^2 directly -- no O(Npts*n^2) pointwise rho(r).
    auto pw=dynamic_cast<const BasisSet::FourierDFT_IBS*>(bs);
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(pw && "PW_Hartree requires a FourierDFT_IBS (plane-wave) basis");
    assert(fd && "PW_Hartree requires a FourierDensity (periodic) charge density");
    double Eh;
    return pw->IntegralHartree(fd->GetFourierDensity(), Eh);
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
PW_XC::PW_XC(const xc_t& xc)
    : itsXc(xc)
{}

// XC via FFT: get rho(r) on the grid (inverse FFT of the density's rho-tilde), apply v_xc pointwise,
// forward-FFT to the matrix.  No O(Npts*n^2) pointwise density sampling.
chmat_t PW_XC::CalcMatrix(const cobs_t* bs, const Spin&, const cDM_CD* cd) const
{
    newCD(cd);
    auto pw=dynamic_cast<const BasisSet::FourierDFT_IBS*>(bs);
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(pw && "PW_XC requires a FourierDFT_IBS (plane-wave) basis");
    assert(fd && "PW_XC requires a FourierDensity (periodic) charge density");
    itsBasis=pw;                               // captured for GetEnergy (which has no basis parameter)
    rvec_t rho=pw->RhoOnGrid(fd->GetFourierDensity());
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=itsXc->GetVxc(rho[q]);   // V(r) = v_xc(rho(r))
    return pw->IntegralPotentialGrid(vxc);
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
    te.Exc += itsBasis->IntegralGrid(exc);     // E_xc = integral eps_xc(rho) rho  (fresh with this cd)
}

std::ostream& PW_XC::Write(std::ostream& os) const
{
    return os << "    PW exchange-correlation potential v_xc(rho(r))." << std::endl;
}

} //namespace
