
#include "QchemTester.H"
#include <MeshParams.H>
#include <SCFIterator.H>
#include <WaveFunction.H>
#include <Hamiltonian.H>
#include <Cluster.H>
#include <BasisSet.H>
#include <Irrep_BS.H>
#include <TotalEnergy.H>
#include <Orbital.H>
#include <ChargeDensity.H>
#include <memory>

PeriodicTable QchemTester::itsPT;

QchemTester::QchemTester()
: itsCluster(0)
, itsBasisSet(0)
, itsHamiltonian(0)
, itsWaveFunction(0)
, itsSCFIterator(0)
, MaxRelErrE(0)
{
    //Cannot call virtual functions from here.
}

QchemTester::~QchemTester()
{
    delete itsBasisSet;
    delete itsWaveFunction;
    delete itsSCFIterator;
}

void QchemTester::Init(double eps)
{
    assert(eps>0.0);
    MaxRelErrE=eps;
    
    assert(itsCluster);
    assert(&*itsCluster);
    itsBasisSet=GetBasisSet(); //SG, PG, Slater
    assert(itsBasisSet);
    itsHamiltonian=GetHamiltonian(itsCluster); //HF,semi HF, DFT all Pol or un-polarized.
    assert(itsHamiltonian);
    itsWaveFunction=GetWaveFunction(itsBasisSet); //Polarized or un-polarized
    assert(itsWaveFunction);
    itsSCFIterator=new SCFIterator(itsWaveFunction,itsHamiltonian);
}

void QchemTester::Iterate(const SCFIterationParams& ipar)
{
    assert(itsSCFIterator);
    itsSCFIterator->Iterate(ipar);
}


double QchemTester::TotalEnergy() const
{
    return GetEnergyBreakdown().GetTotalEnergy();
}

EnergyBreakdown QchemTester::GetEnergyBreakdown() const
{
    assert(itsHamiltonian);
    DM_CD* cd=itsWaveFunction->GetChargeDensity();
    EnergyBreakdown te=itsHamiltonian->GetTotalEnergy(cd);
    delete cd;
    return te;
}

double QchemTester::TotalCharge() const
{
    return itsWaveFunction->GetChargeDensity()->GetTotalCharge();
}

const Orbitals* QchemTester::GetOrbitals(const Irrep_QNs& qns) const
{
    return itsWaveFunction->GetOrbitals(qns);
}

const Orbital* QchemTester::GetOrbital(size_t index, const Irrep_QNs& qns) const
{
    const Orbitals* orbs=GetOrbitals(qns);
    assert(index<(size_t)orbs->GetNumOrbitals());
    const Orbital* o=0;
    for (auto oi:orbs->Iterate<Orbital>())
    {
        if (index==0) 
        {
            o=oi;
            break;
        }
        index--;
    }
    return o;
}

#include <cmath> //fabs
double QchemTester::RelativeError(double E,bool quiet) const
{
    double error=(E-TotalEnergy())/E;
    if (!quiet)
    {
        std::cout.precision(9);
        std::cout << "E relative error=" << error*100.0 << "%, ";
        std::cout.precision(2);
        if (fabs(error)>1e-7)
            std::cout << error*1e6 << "(ppm)" << std::endl;
        else
            std::cout << error*1e9 << "(ppb)" << std::endl;
    }
    return fabs(error);
}



double QchemTester::RelativeHFError(bool quiet) const
{
    double E_HF=itsPT.GetEnergyHF(itsCluster->GetNuclearCharge());
    return RelativeError(E_HF,quiet);
}

double QchemTester::RelativeDFTError(bool quiet) const
{
    double E_DFT=itsPT.GetEnergyDFT(itsCluster->GetNuclearCharge());
    return RelativeError(E_DFT,quiet);
}

int QchemTester::GetZ() const
{
    return GetCluster()->GetNuclearCharge();
}

QchemTester::symv_t QchemTester::GetSymmetries() const
{
    return itsBasisSet->GetSymmetries();
}

// Fit_IBS* QchemTester::GetCBasisSet() const
// {
//     return itsBasisSet->CreateCDFitBasisSet(itsCluster.get());
// }

// Fit_IBS* QchemTester::GetXBasisSet() const
// {
//     return itsBasisSet->CreateVxcFitBasisSet(itsCluster.get());    
// }


#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
BasisSet* SG_OBasis::GetBasisSet () const
{
    BasisSet* bs=new Atoml::Gaussian::BasisSet(N,emin,emax,Lmax);
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    //std::cout << *bs << std::endl;
    return  bs;
}

#include "Imp/BasisSet/Atom/l/Slater_BS.H"
BasisSet* SL_OBasis::GetBasisSet () const
{
    BasisSet* bs=new Atoml::Slater::BasisSet(N,emin,emax,Lmax);
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    return bs;
}

#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
BasisSet* SLm_OBasis::GetBasisSet () const
{
    // BasisSet* bs=new Atom_ml::Slater::BasisSet(N,emin,emax,Lmax);
    BasisSet* bs=new Atom_ml::Slater::BasisSet(N,emin,emax,*GetElectronConfiguration());
    
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    return bs;
}

#include "Imp/BasisSet/Atom/kappa/Slater_BS.H"
BasisSet* SLmj_OBasis::GetBasisSet () const
{
    assert(N>0);
    BasisSet* bs=new Atom_kappa::Slater::BasisSet(N,emin,emax,Lmax);
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    assert(bs->GetNumFunctions()>0);
    return bs;
}

#include "Imp/BasisSet/Atom/kappa/Gaussian_BS.H"
BasisSet* SG_RKB_OBasis::GetBasisSet () const
{
    assert(N>0);
    BasisSet* bs=new Atom_kappa::Gaussian::BasisSet(N,emin,emax,Lmax);
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    assert(bs->GetNumFunctions()>0);
    return bs;
}

#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/BasisSet/Molecule/PolarizedGaussian/Readers/Gaussian94.H"

BasisSet* SGm_OBasis::GetBasisSet () const
{
//    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
//    const Cluster* cl=GetCluster();
//    Atom* a=*cl->begin();
//    SphericalGaussian_m::BasisSet* bs=new SphericalGaussian_m::BasisSet(lap,&reader,a);
    BasisSet* bs=new Atom_ml::Gaussian::BasisSet(N,emin,emax,*GetElectronConfiguration());
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    return  bs;
}


#include "Imp/BasisSet/Molecule/PolarizedGaussian/BasisSet.H"

BasisSet* PG_OBasis::GetBasisSet () const
{
    if (N==0)
    {
        PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
        BasisSet* bs=new PolarizedGaussian::BasisSet(&reader,GetCluster());  
        bs->Set(lap);
        // StreamableObject::SetToPretty();
        //std::cout << *bs << std::endl;
        return bs;
//        return new PolarizedGaussian::BasisSet(lap, &reader,GetCluster());        
    }
    else
    {
        BasisSet* bs=new PolarizedGaussian::BasisSet(N,emin,emax,LMax,GetCluster());
        bs->Set(lap);
        // StreamableObject::SetToPretty();
        //std::cout << *bs << std::endl;
        return bs;
//        return new PolarizedGaussian::BasisSet(lap, N,emin,emax,LMax,GetCluster());   
    }
}

#include "Imp/BasisSet/Atom/l/BSpline_BS.H"
BasisSet* BS_OBasis::GetBasisSet () const
{
    BasisSet* bs=new Atoml::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    //std::cout << *bs << std::endl;
    return  bs;
}

#include "Imp/BasisSet/Atom/ml/BSpline_BS.H"
BasisSet* BSm_OBasis::GetBasisSet () const
{
    BasisSet* bs=new Atom_ml::BSpline::BasisSet<6>(N,rmin,rmax,*GetElectronConfiguration());
    bs->Set(lap);
    // StreamableObject::SetToPretty();
    // std::cout << *bs << std::endl;
    return  bs;
}



#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
TestAtom::TestAtom(int Z, int q) : ec(Z-q) //Pass in # of electrons.
{
    Cluster* cl=new Molecule;
    cl->Insert(new Atom(Z,q,Vector3D<double>(0,0,0)));
    itsCluster=cl_t(cl);
};

MeshParams TestAtom::GetMeshParams() const
{
    return MeshParams({qchem::MHL,50,3,2.0,qchem::Gauss,1,0,0,2});
}

void TestMolecule::Init(Molecule* m)
{
    assert(m);
    itsCluster=cl_t(m);
    ec=Molecule_EC(m->GetNumElectrons());
}

MeshParams TestMolecule::GetMeshParams() const
{
    return MeshParams({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,2});
}

    

#include "Imp/WaveFunction/UnPolarized_WF.H"
WaveFunction* TestUnPolarized::GetWaveFunction(const BasisSet* bs) const
{
    return new UnPolarized_WF(bs,GetElectronConfiguration());
}

#include "Imp/WaveFunction/Polarized_WF.H"

WaveFunction* TestPolarized::GetWaveFunction(const BasisSet* bs) const
{   
    assert(bs->GetNumFunctions()>0);
    return new Polarized_WF(bs,GetElectronConfiguration() );
}

