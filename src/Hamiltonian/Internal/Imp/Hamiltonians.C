// File:: Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <memory>
module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.BasisSet;
import qchem.Fit_IBS;
import qchem.Cluster;

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


Ham_DFT_U::Ham_DFT_U(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
       
    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

Ham_HF_P::Ham_HF_P(const cl_t& cl)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new VxcPol);
}


Ham_DFT_P::Ham_DFT_P(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
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
    //Add(new Vnn(cl));
}

Ham_DHF::Ham_DHF(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(cl));
    Add(new Ven(cl));
    // Add(new Vee());
    // Add(new VxcPol());
}