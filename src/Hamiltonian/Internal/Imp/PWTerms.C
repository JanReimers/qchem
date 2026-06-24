// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <cassert>
#include <iostream>
#include <memory>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ScalarFunction;

namespace qchem::Hamiltonian
{

namespace
{
// Compose an LDA functional with the density to make a real-space scalar field to hand the basis.
// The density (cDM_CD) IS-A ScalarFunction<double>, so rho(r) = (*itsRho)(r).
class VxcField : public virtual ScalarFunction<double>
{
public:
    VxcField(const ExFunctional* xc, const ScalarFunction<double>* rho) : itsXc(xc), itsRho(rho) {}
    virtual double  operator()(const rvec3_t& r) const {return itsXc->GetVxc((*itsRho)(r));}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
private:
    const ExFunctional* itsXc;  const ScalarFunction<double>* itsRho;
};
class EpsXcRhoField : public virtual ScalarFunction<double>   // the XC ENERGY density eps_xc(rho) rho
{
public:
    EpsXcRhoField(const ExFunctional* xc, const ScalarFunction<double>* rho) : itsXc(xc), itsRho(rho) {}
    virtual double  operator()(const rvec3_t& r) const {double ro=(*itsRho)(r); return itsXc->GetEpsXc(ro)*ro;}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
private:
    const ExFunctional* itsXc;  const ScalarFunction<double>* itsRho;
};
} //anon

PW_External::PW_External(const cl_t& cl)
    : cStatic_HT_Imp()
    , theStructure(cl)
{
    assert(cl->GetNumAtoms()>0);
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
    auto pw=dynamic_cast<const BasisSet::DFTPotential_IBS<dcmplx>*>(bs);
    assert(pw && "PW_Hartree requires a DFTPotential_IBS (e.g. plane-wave) basis");
    double Eh;
    return pw->IntegralHartree(*cd, Eh);   // cd IS-A ScalarFunction<double> (the density rho(r))
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

chmat_t PW_XC::CalcMatrix(const cobs_t* bs, const Spin&, const cDM_CD* cd) const
{
    newCD(cd);
    auto pw=dynamic_cast<const BasisSet::DFTPotential_IBS<dcmplx>*>(bs);
    assert(pw && "PW_XC requires a DFTPotential_IBS (e.g. plane-wave) basis");
    itsBasis=pw;                               // captured for GetEnergy (which has no basis parameter)
    VxcField vxc(itsXc.get(), cd);             // V(r) = v_xc(rho(r))
    return pw->IntegralPotential(vxc);
}

void PW_XC::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    assert(itsBasis && "PW_XC::GetEnergy called before GetMatrix");
    EpsXcRhoField exc(itsXc.get(), cd);
    te.Exc += itsBasis->Integral(exc);         // E_xc = integral eps_xc(rho) rho  (fresh with this cd)
}

std::ostream& PW_XC::Write(std::ostream& os) const
{
    return os << "    PW exchange-correlation potential v_xc(rho(r))." << std::endl;
}

} //namespace
