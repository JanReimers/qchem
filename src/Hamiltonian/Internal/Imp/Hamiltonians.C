// File:: Hamiltonians.C  Create fully implemented Hamiltonians
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

Ham_1E::Ham_1E(const cl_t& cl) 
{
    InsertStandardTerms(cl);
}

Ham_HF_U::Ham_HF_U(const cl_t& cl) 
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new Vxc);
}


Ham_DFT_U::Ham_DFT_U(const cl_t& cl,double alpha_ex, const MeshParams& mp, const bs_t* bs)
    : Ham_DFT_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(cl);
       
    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

// Dirac exchange + VWN5 correlation as separate terms (correct E_c = integral eps_c rho), sharing one
// Vxc fit basis so the 3-centre integrals are computed once.
Ham_DFTcorr_U::Ham_DFTcorr_U(const cl_t& cl, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(cl);

    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));

    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get())); // ONE Vxc fit basis, shared X and C
    FittedVxc::ex_t exch(new SlaterExchange(2.0/3.0));            // Dirac exchange (alpha = 2/3)
    Add(new FittedVxc  (XFitBasis, exch, m));
    FittedVxc::ex_t corr(new VWN_Correlation());                 // VWN5 correlation
    Add(new FittedVcorr(XFitBasis, corr, m));
}

Ham_HF_P::Ham_HF_P(const cl_t& cl)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new VxcPol);
}


Ham_DFT_P::Ham_DFT_P(const cl_t& cl,double alpha_ex, const MeshParams& mp, const bs_t* bs)
    : Ham_DFT_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const bs_t* bs)
{
    InsertStandardTerms(cl);
    FittedVee::bs_t CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t  m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));
    
    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
    
}


Ham_DHF_1E::Ham_DHF_1E(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Ven(cl));
    //Add(new Vnn(cl));
}

Ham_DHF_U::Ham_DHF_U(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(cl));
    Add(new Ven(cl));
    Add(new Vee());
    Add(new Vxc());
}
Ham_DHF_P::Ham_DHF_P(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(cl));
    Add(new Ven(cl));
    Add(new Vee());
    Add(new VxcPol());
}

} //namespace
