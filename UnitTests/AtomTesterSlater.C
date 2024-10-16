// File AtomTester.C Member functions for the atom tester class.

#include "AtomTesterSlater.H"
#include "Cluster.H"
#include "BasisSet.H"
#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "Mesh/RadialMesh/LogRadialMesh.H"
//#include "Mesh/AngularMesh/GaussLegendreAngularMesh.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"


AtomTesterSlater::AtomTesterSlater(size_t n, double emin, double emax)
    : BaseTester()
    , itsNbasis(n)
    , itsEmin(emin)
    , itsEmax(emax)
    , itsLmax(0)
    , itsSCFIParams({40,1e-3,1.0,0.0}) //MaxITer, MinDeltaRo, StartingRelaxRo, kT
{
}

void AtomTesterSlater::Init(Atom* atom)
{
    assert(itsCluster->GetNumAtoms()==0);
    itsCluster->Insert(atom);
}

void AtomTesterSlater::Init(Atom* atom, int Lmax,double spin)
{
    Init(atom);
    itsLmax=Lmax;
    auto bg=new Slater::BasisSet(itsLAParams,itsNbasis,itsEmin,itsEmax,itsLmax);
    BaseTester::Init(bg,spin);
}

void AtomTesterSlater::Init(int NBasis, int Lmax, double spin, const LAParams& lap)
{
    itsLmax=Lmax;
    itsNbasis=NBasis;
    auto bg=new Slater::BasisSet(itsLAParams,itsNbasis,itsEmin,itsEmax,itsLmax);
    BaseTester::Init(bg,spin,lap);
}

IrrepBasisSet* AtomTesterSlater::GetCbasisSet() const
{
    return new Slater::IrrepBasisSet(itsLAParams,GetDatabase(),itsNbasis,itsEmin*2.0,itsEmax*2.0,0);
}

IrrepBasisSet* AtomTesterSlater::GetXbasisSet() const
{
    return new Slater::IrrepBasisSet(itsLAParams,GetDatabase(),itsNbasis,itsEmin*2.0/3.0,itsEmax*2.0/3.0,0);
}

rc_ptr<Mesh> AtomTesterSlater::GetIntegrationMesh() const
{
    RadialMesh*            rm=new MHLRadialMesh(50,2U,2.0); //mem leak
    AngularMesh*           am=new GaussAngularMesh(1);      //mem leak
    return new AtomMesh(*rm,*am); //why not own?
}

