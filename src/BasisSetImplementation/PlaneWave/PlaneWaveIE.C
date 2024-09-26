// File: PlaneWaveIE.C  Here is where all the integral get calculated.

#pragma implementation

#include "BasisSet/TBasisSet.H"
#include "BasisSet/TBasisSetBrowser.H"
#include "BasisSet/IntegralDataBase.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveIE.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveBF.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveBS.H"
#include "BasisSetImplementation/PlaneWave/BlochQN.H"
#include "Cluster/ClusterBrowser.H"
#include "Cluster/Atom.H"
#include "oml/list.h"
#include "oml/io3d.h"
#include "oml/list_io.h"
#include "Misc/ptr_vector.h"
#include "Misc/MatrixList.H"
#include "Misc/DFTDefines.H"

#include <cmath>
#include <cassert>
#include <cstdlib>

const double epsilon=0.001;//Vectors with |G|<epsilon are considered zero vectors.

//-----------------------------------------------------------------
//
//  Construction zone.
//

PlaneWaveIE::PlaneWaveIE()
    : IntegralEngineImplementation<std::complex<double> >()
{};

PlaneWaveIE::PlaneWaveIE(const PlaneWaveIE& pgie)
    : IntegralEngineImplementation<std::complex<double> >(pgie)
    , itsGs(pgie.itsGs)
{};

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
IntegralEngine<std::complex<double> >* PlaneWaveIE::Clone() const
{
    return new PlaneWaveIE(*this);
}

std::ostream& PlaneWaveIE::Write(std::ostream& os) const
{
    IntegralEngineImplementation<std::complex<double> >::Write(os);
    return os << itsGs << itsK << itsRLCell;
}

std::istream& PlaneWaveIE::Read (std::istream& is)
{
    IntegralEngineImplementation<std::complex<double> >::Read(is);
    return is >> itsGs >> itsK >> itsRLCell;
}

void PlaneWaveIE::Insert(TBasisSet<std::complex<double> >* bs)
{
    assert(bs);
    IntegralEngineImplementation<std::complex<double> >::Insert(bs);
    for (BasisSetBrowser b(*bs); b; b++)
    {
        const PlaneWaveBF* pwbf=dynamic_cast<const PlaneWaveBF*>(&*b);
        assert(pwbf);
        itsGs.push_back(pwbf->itsG);
    }
    const BlochQN* bqn=dynamic_cast<const BlochQN*>(&bs->GetQuantumNumber());
    assert(bqn);
    itsK=bqn->GetK();

    const PlaneWaveBS* pwbs=dynamic_cast<const PlaneWaveBS*>(bs);
    assert(pwbs);
    itsRLCell=pwbs->itsRLCell;
}
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
PlaneWaveIE::SMat PlaneWaveIE::MakeOverlap() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Unit(ret);
    return ret;
}

PlaneWaveIE::Mat PlaneWaveIE::MakeOverlap(const TBasisSet<std::complex<double> >& theBasisSet) const
{
    CheckBasisSet();
    const PlaneWaveIE* OtherIE;
    const IntegralEngine<std::complex<double> >*      theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const PlaneWaveIE*>(theOtherIE);
    assert(OtherIE);

    Mat ret(itsN,OtherIE->itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
    {
        int i2=1;
        for (List<RVec3>::const_iterator b2(OtherIE->itsGs.begin()); b2!=OtherIE->itsGs.end(); b2++,i2++)
            if (*b1==*b2) ret(i1,i2)=1.0;
    }

    return ret;
}

PlaneWaveIE::Vec PlaneWaveIE::MakeOverlap(const ScalarFunction<double>& f) const
{
    CheckBasisSet();
    const PlaneWaveBF* pwbf= dynamic_cast<const PlaneWaveBF*>(&f);
    assert(pwbf);

    Vec ret(itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
        if (*b1==pwbf->itsG) ret(i1)=1.0;

    return ret;
}

//------------------------------------------------------------------------------------
//
//  Calculates overlap matricies <ab|c> for a whole set basis functions c.
//  Everything gets normalized and the elements below the diagonal are defined.
//
void PlaneWaveIE::MakeOverlap3C(MList& ret, const TBasisSet<std::complex<double> >& otherBS) const
{
    ret.Empty();
    for (TBasisSetBrowser<std::complex<double> > b(otherBS); b; b++)
    {
        SMat m(itsN,itsN);
        Fill(m,std::complex<double>(0.0));
        const PlaneWaveBF* pwbf=dynamic_cast<const PlaneWaveBF*>(&*b);
        assert(pwbf);
        RVec3 G3=pwbf->itsG;

        int i1=1;
        for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
        {
            int i2=i1;
            for (List<RVec3>::const_iterator b2(b1); b2!=itsGs.end(); b2++,i2++)
                if (itsRLCell.GetDistance(*b1-*b2+G3)<epsilon) m(i1,i2)=1.0;
        }
        ret.Add(m);
    }

    ret.Clear();
}


//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//
PlaneWaveIE::SMat PlaneWaveIE::MakeRepulsion() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
    {
        double G=itsRLCell.GetDistance(*b1);
        if (G>epsilon) ret(i1,i1)=4*Pi/(G*G);
    }

    return ret;
}

PlaneWaveIE::Mat PlaneWaveIE::MakeRepulsion(const TBasisSet<std::complex<double> >& theBasisSet) const
{
    CheckBasisSet();
    const PlaneWaveIE* OtherIE;
    const IntegralEngine<std::complex<double> >*      theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const PlaneWaveIE*>(theOtherIE);
    assert(OtherIE);

    Mat ret(itsN,OtherIE->itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
    {
        int i2=1;
        double G=itsRLCell.GetDistance(*b1);
        for (List<RVec3>::const_iterator b2(OtherIE->itsGs.begin()); b2!=OtherIE->itsGs.end(); b2++,i2++)
            if (*b1==*b2 && G>epsilon) ret(i1,i2)=4*Pi/(G*G);
    }

    return ret;
}

PlaneWaveIE::Vec PlaneWaveIE::MakeRepulsion(const ScalarFunction<double>& f) const
{
    const PlaneWaveBF* pwbf= dynamic_cast<const PlaneWaveBF*>(&f);
    assert(pwbf);

    Vec ret(itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    double G=itsRLCell.GetDistance(pwbf->itsG);

    if (G>epsilon)
        for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
            if (*b1==pwbf->itsG) ret(i1)=4*Pi/(G*G);

    return ret;
}


//------------------------------------------------------------------------------------
//
//  Calculates repulsion matricies <ab|1/r12|c> for a whole set basis functions c.
//  Everything gets normalized and the elements below the diagonal are defined.
//
void PlaneWaveIE::MakeRepulsion3C(MList& ret, const TBasisSet<std::complex<double> >& otherBS) const
{
    ret.Empty();
    for (TBasisSetBrowser<std::complex<double> > b(otherBS); b; b++)
    {
        SMat m(itsN,itsN);
        Fill(m,std::complex<double>(0.0));
        const PlaneWaveBF* pwbf=dynamic_cast<const PlaneWaveBF*>(&*b);
        assert(pwbf);
        RVec3 G3=pwbf->itsG;
        double G=itsRLCell.GetDistance(G3);

        int i1=1;
        if (G>epsilon)
            for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
            {
                int i2=i1;
                for (List<RVec3>::const_iterator b2(b1); b2!=itsGs.end(); b2++,i2++)
                    if (!(*b1-*b2+G3)<epsilon) m(i1,i2)=4*Pi/(G*G);
            }

        ret.Add(m);
    }
    ret.Clear();
}

void PlaneWaveIE::MakeRepulsion4C(ERIList& eris) const
{
    std::cerr << "PlaneWaveIE::MakeRepulsion4C 4 center plane wave integrals are not supported" << std::endl;
    exit(-1);
}


//----------------------------------------------------------------------------------------
//
//  Special integrals
//
PlaneWaveIE::SMat PlaneWaveIE::MakeKinetic() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
    {
        double G=itsRLCell.GetDistance(*b1+itsK);
        ret(i1,i1)=G*G;
    }
    return ret;
}

double FormFactor(const Cluster&, const RVec3& G);

PlaneWaveIE::SMat PlaneWaveIE::MakeNuclear(const Cluster& theCluster) const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,std::complex<double>(0.0));
    int i1=1;
    for (List<RVec3>::const_iterator b1(itsGs.begin()); b1!=itsGs.end(); b1++,i1++)
    {
        int i2=i1;
        for (List<RVec3>::const_iterator b2(b1); b2!=itsGs.end(); b2++,i2++)
        {
            RVec3 G=*b1-*b2;
            double mG=itsRLCell.GetDistance(G);
            if (mG > epsilon) ret(i1,i2)=4*Pi*FormFactor(theCluster,G)/(mG*mG);
        }
    }
    return ret;
}

double FormFactor(const Cluster& cl, const RVec3& G)
{
    double ret=0;
    for (ClusterBrowser b(cl); b; b++)
        ret+=b->itsZ*cos(2*Pi*G*b->itsR);

    return ret;
}



