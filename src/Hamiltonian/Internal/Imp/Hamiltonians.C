// File:: Hamiltonian/Internal/Imp/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <memory>
module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Types;
import qchem.Structure;

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
    Add(new Vxc);
}


Ham_DFT_U::Ham_DFT_U(const st_t& st,double alpha_ex, const MeshParams& mp, const bs_t* bs)
    : Ham_DFT_U(st,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const st_t& st,ExFunctional* ex, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(st);
       
    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get()));
    FittedVee::mesh_t m(st->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,st->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get()));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

// Dirac exchange + VWN5 correlation as separate terms (correct E_c = integral eps_c rho), sharing one
// Vxc fit basis so the 3-centre integrals are computed once.
Ham_DFTcorr_U::Ham_DFTcorr_U(const st_t& st, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(st);

    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get()));
    FittedVee::mesh_t m(st->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,st->GetNumElectrons()));

    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get())); // ONE Vxc fit basis, shared X and C
    FittedVxc::ex_t exch(new SlaterExchange(2.0/3.0));            // Dirac exchange (alpha = 2/3)
    Add(new FittedVxc  (XFitBasis, exch, m));
    FittedVxc::ex_t corr(new VWN_Correlation());                 // VWN5 correlation
    Add(new FittedVcorr(XFitBasis, corr, m));
}

Ham_HF_P::Ham_HF_P(const st_t& st)
{
    InsertStandardTerms(st);
    Add(new Vee);
    Add(new VxcPol);
}


Ham_DFT_P::Ham_DFT_P(const st_t& st,double alpha_ex, const MeshParams& mp, const bs_t* bs)
    : Ham_DFT_P(st,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const st_t& st,ExFunctional* ex, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(st);
    FittedVee::bs_t CFitBasis(bs->CreateCDFitBasisSet(st.get()));
    FittedVee::mesh_t  m(st->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,st->GetNumElectrons()));
    
    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::bs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get()));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
    
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
    Add(new Vxc());
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
