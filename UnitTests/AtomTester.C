// File AtomTester.C Member functions for the atom tester class.

#include "AtomTester.H"
#include "Cluster.H"
#include "BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "Mesh/RadialMesh/LogRadialMesh.H"
//#include "Mesh/AngularMesh/GaussLegendreAngularMesh.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"


AtomTester::AtomTester(size_t n, double emin, double emax)
    : BaseTester()
    , itsNbasis(n)
    , itsEmin(emin)
    , itsEmax(emax)
    , itsLmax(0)
    , itsSCFIParams({40,1e-3,1.0,0.0}) //MaxITer, MinDeltaRo, StartingRelaxRo, kT
{
}

void AtomTester::Init(Atom* atom)
{
    assert(itsCluster->GetNumAtoms()==0);
    itsCluster->Insert(atom);
}

void AtomTester::Init(Atom* atom, int Lmax,double spin)
{
    Init(atom);
    itsLmax=Lmax;
    auto bg=new SphericalGaussian::BasisSet(itsLAParams,itsNbasis,itsEmin,itsEmax,itsLmax);
    BaseTester::Init(bg,spin);
}

void AtomTester::Init(int NBasis, int Lmax, double spin, const LAParams& lap)
{
    itsLmax=Lmax;
    itsNbasis=NBasis;
    auto bg=new SphericalGaussian::BasisSet(itsLAParams,itsNbasis,itsEmin,itsEmax,itsLmax);
    BaseTester::Init(bg,spin,lap);
}

IrrepBasisSet* AtomTester::GetCbasisSet() const
{
    return new SphericalGaussian::IrrepBasisSet(itsLAParams,new HeapDB<double>,itsNbasis,itsEmin*2.0,itsEmax*2.0,0);
}

IrrepBasisSet* AtomTester::GetXbasisSet() const
{
    return new SphericalGaussian::IrrepBasisSet(itsLAParams,new HeapDB<double>,itsNbasis,itsEmin*2.0/3.0,itsEmax*2.0/3.0,0);
}

Mesh* AtomTester::GetIntegrationMesh() const
{
    RadialMesh*            rm=new MHLRadialMesh(50,2U,2.0);
    AngularMesh*           am=new GaussAngularMesh(1);
    return new AtomMesh(*rm,*am);
}

