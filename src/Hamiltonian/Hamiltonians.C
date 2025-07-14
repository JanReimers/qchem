// File:: Hamiltonians.C  Create fully implemented Hamiltonians

#include <memory>

#include "Hamiltonians.H"
#include "Vee.H"
#include "Vxc.H"
#include <BasisSet/BasisSet.H>
#include <BasisSet/Fit_IBS.H>

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

#include "SlaterExchange.H"
#include "FittedVxc.H"
#include "ExchangeFunctional.H" 
#include "FittedVee.H"
#include <Cluster/Cluster.H>

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

#include "VxcPol.H"
Ham_HF_P::Ham_HF_P(const cl_t& cl)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new VxcPol);
}

#include "FittedVxcPol.H"

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

#include "DiracKinetic.H"
#include "RestMass.H"
#include "Vnn.H"
#include "Ven.H"

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