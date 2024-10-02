// File AtomTester.C Member functions for the atom tester class.

#include "AtomTester.H"
#include "Cluster.H"
#include "BasisGroup.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
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
    BaseTester::Init(spin);
}

void AtomTester::Init(int NBasis, int Lmax, double spin, const LinearAlgebraParams& lap)
{
    itsLmax=Lmax;
    itsNbasis=NBasis;
    BaseTester::Init(spin,lap);
}

void AtomTester::LoadOrbitalBasisSet()
{
    assert(itsBasisGroup);
    for (int l=0; l<=itsLmax; l++)
    {
        BasisSet* bs=new SphericalGaussianBS(itsLAParams,new HeapDB<double>,itsNbasis,itsEmin,itsEmax,l);
        itsBasisGroup->Insert(bs);
    }
}

BasisSet* AtomTester::GetCbasisSet() const
{
    return new SphericalGaussianBS(itsLAParams,new HeapDB<double>,itsNbasis,itsEmin*2.0,itsEmax*2.0,0);
}

BasisSet* AtomTester::GetXbasisSet() const
{
    Mesh* mesh=GetIntegrationMesh();
    assert(mesh);
    return new SphericalGaussianBS(itsLAParams,new HeapDB<double>,itsNbasis,itsEmin*2.0/3.0,itsEmax*2.0/3.0,0,mesh);
}

Mesh* AtomTester::GetIntegrationMesh() const
{
    RadialMesh*            rm=new MHLRadialMesh(50,2U,2.0);
    AngularMesh*           am=new GaussAngularMesh(1);
    return new AtomMesh(*rm,*am);
}

