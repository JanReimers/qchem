// File: ContractedGaussianRF.C  Contracted Gaussian in 3D space.


#include "Misc/Polarization.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianH3.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite1.H"
#include "Mesh/MeshBrowser.H"
#include "Misc/ptr_vector.h"
#include "Misc/ptrvector_io.h"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
//#include "oml/vector_io.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

//#######################################################################
//
//  Contracted gaussian implementation
//
typedef optr_vector<RadialFunction*>::      iterator ITER;
typedef optr_vector<RadialFunction*>::const_iterator CITER;
typedef Vector<double>              ::const_iterator CVITER;
typedef Vector<double>              ::      iterator VITER;

ContractedGaussianRF::ContractedGaussianRF()
    : RadialFunctionImplementation()
    , itsCoeffs             ( )
    , itsPrimatives         ( )
{};

ContractedGaussianRF::ContractedGaussianRF(const Vector<double>& theCoeffs, ptr_vector<RadialFunction*>& its_rfs)
    : RadialFunctionImplementation( its_rfs[0]->GetCenter (), its_rfs[0]->GetL())
    , itsCoeffs    (theCoeffs)
    , itsPrimatives(         )
{
//
//  Plant all the radials in the list.
//
    for (ptr_vector<RadialFunction*>::iterator b(its_rfs.begin()); b!=its_rfs.end(); b++) itsPrimatives.push_back(&b);
//
//  Absorb normalizations into the contraction coefficients.
//
    CITER b(itsPrimatives.begin());
    VITER c(itsCoeffs.begin());
    for (; b!=itsPrimatives.end()&&c!=itsCoeffs.end(); b++,c++) *c *= b->GetNormalization(Polarization(GetL(),0,0));
    Check(); //Make sure all L's and centers are the same.
};


Hermite1* ContractedGaussianRF::MakeH1() const
{
    Hermite1* ret=new Hermite1;
    CITER b(itsPrimatives.begin());
    CVITER c(itsCoeffs.begin());
    for (; b!=itsPrimatives.end()&&c!=itsCoeffs.end(); b++,c++) ret->Add(b->GetH1(),*c);
    return ret;
}

bool ContractedGaussianRF::operator==(const RadialFunction& rf) const
{
    bool ret=false;
    const ContractedGaussianRF* cg = dynamic_cast<const ContractedGaussianRF*>(&rf);
    if (cg)
    {
        bool base     = RadialFunctionImplementation::operator==(rf);
        bool coeff    = itsCoeffs.size()==cg->itsCoeffs.size() && itsCoeffs==cg->itsCoeffs;
        bool index    = true;

        CITER b1(    itsPrimatives.begin());
        CITER b2(cg->itsPrimatives.begin());
        for (; b1!=itsPrimatives.end()&&b2!=cg->itsPrimatives.end(); b1++,b2++)  index = index && b1->GetID() == b2->GetID();

        ret = base && coeff && index;
    }
    return ret;
}

double ContractedGaussianRF::GetNormalization(const Polarization& p) const
{
    BasisFunctionBlock block(Clone(),1);
    block.Add(p);
    SMatrix<double> ret(1,1);
    Fill(ret,0.0);
    BasisFunctionBlockPair bfbp(&block,&block);
    this->Get2CenterIntegrals(RadialFunction::Overlap2C,bfbp,ret,NULL,1.0);
    return 1.0/sqrt(ret(1,1));
}

double ContractedGaussianRF::GetCharge(const Polarization& p) const
{
    double ret=0.0;
    CITER b(itsPrimatives.begin());
    CVITER c(itsCoeffs.begin());
    for (; b!=itsPrimatives.end()&&c!=itsCoeffs.end(); b++,c++)
    {
        ret += *c * b->GetCharge(p);
    }
    return ret;
}

void ContractedGaussianRF::Get2CenterIntegrals(Types2C type, BFBP& p, SMat& ret, const Cluster* cl, double scale) const
{
    CITER ri(itsPrimatives.begin());
    CVITER ci(itsCoeffs.begin());
    for (; ri!=itsPrimatives.end()&&ci!=itsCoeffs.end(); ri++,ci++)
        ri->Get2CenterIntegrals(type,p,ret,cl,*ci*scale);
}
void ContractedGaussianRF::Get2CenterIntegrals(Types2C type, BFBP& p, Mat& ret, double scale) const
{
    CITER ri(itsPrimatives.begin());
    CVITER ci(itsCoeffs.begin());
    for (; ri!=itsPrimatives.end()&&ci!=itsCoeffs.end(); ri++,ci++)
        ri->Get2CenterIntegrals(type,p,ret,*ci*scale);
}

void ContractedGaussianRF::Get3CenterIntegrals(Types3C type, BFBT& t, std::vector<SMat>& ret, double scale)
{
    CITER ri(itsPrimatives.begin());
    CVITER ci(itsCoeffs.begin());
    for (; ri!=itsPrimatives.end()&&ci!=itsCoeffs.end(); ri++,ci++)
        ri->Get3CenterIntegrals(type,t,ret,*ci*scale);
}

void ContractedGaussianRF::GetRepulsion4C(BFBQ& q, ERIList& eris, double scale)
{
    CITER ri(itsPrimatives.begin());
    CVITER ci(itsCoeffs.begin());
    for (; ri!=itsPrimatives.end()&&ci!=itsCoeffs.end(); ri++,ci++)
        ri->GetRepulsion4C(q,eris,*ci*scale);
}


Hermite3* ContractedGaussianRF::GetH3(const RadialFunction& r1, const RadialFunction& r2) const
{
    ContractedGaussianH3* ret = new ContractedGaussianH3(itsCoeffs);
    for (CITER b(itsPrimatives.begin()); b!=itsPrimatives.end(); b++) ret->Insert( b->GetH3(r1,r2) );
    return ret;
}

Matrix<double> ContractedGaussianRF::GetAux(const std::vector<Polarization>& N, const std::vector<Polarization>& Pc,
        int LP, double AlphaP,const RVec3& P) const
{
    Matrix<double> ret(N.size(),Pc.size());
    Fill(ret,0.0);

    CITER b(itsPrimatives.begin());
    CVITER c(itsCoeffs.begin());
    for (; b!=itsPrimatives.end()&&c!=itsCoeffs.end(); b++,c++)
    {
        ret += *c * b->GetAux(N,Pc,LP,AlphaP,P);
    }
    return ret;
}

std::ostream& ContractedGaussianRF::Write(std::ostream& os) const
{
    if (Binary()) os << itsCoeffs << itsPrimatives;
    if (Ascii ()) os << itsCoeffs << std::endl << itsPrimatives << std::endl;
    if (Pretty())
    {
        os << "Contracted {";
        for(CITER p(itsPrimatives.begin()); p!=itsPrimatives.end(); p++)
        {
            std::ostringstream s;
            s << *p << std::ends;
            std::string ss(s.str());
            std::string num=ss.substr(10);
            os << num << " ";
        }
        os << "}";
    }
    if (!Pretty()) RadialFunctionImplementation::Write(os);
    return os;
}

std::istream& ContractedGaussianRF::Read(std::istream& is)
{
    is >> itsCoeffs >> itsPrimatives;
    RadialFunctionImplementation::Read(is);
    Check(); //Make all L's and centers are the same.
    return is;
}

double ContractedGaussianRF::operator()(const RVec3& r) const
{
    double ret=0;
    CVITER c(itsCoeffs.begin());
    CITER p(itsPrimatives.begin());
    for (; p!=itsPrimatives.end()&&c!=itsCoeffs.end(); c++,p++) ret+=*c * p->operator()(r);
    return ret;
}

RVec3 ContractedGaussianRF::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    CVITER c(itsCoeffs.begin());
    CITER p(itsPrimatives.begin());
    for (; p!=itsPrimatives.end()&&c!=itsCoeffs.end(); c++,p++) ret+=*c * p->Gradient(r);
    return ret;
}

void ContractedGaussianRF::Eval(const Mesh& mesh, Vector<double>& vec) const
{
    CVITER c(itsCoeffs.begin());
    CITER p(itsPrimatives.begin());
    for (; p!=itsPrimatives.end()&&c!=itsCoeffs.end(); c++,p++)
    {
        Vector<double>::iterator i(vec.begin());
        Vector<double>::const_iterator  v((*p)(mesh).begin());
        for (; i!=vec.end()&&v; i++,v++) *i += (*c) * (*v);
    }
}

void ContractedGaussianRF::EvalGrad(const Mesh& mesh, Vector<RVec3>& vec) const
{
    CVITER c(itsCoeffs.begin());
    CITER p(itsPrimatives.begin());
    for (; p!=itsPrimatives.end()&&c!=itsCoeffs.end(); c++,p++)
    {
        Vector<RVec3>::iterator i(vec.begin());
        Vector<RVec3>::const_iterator  v(p->Gradient(mesh).begin());
        for (; i&&v; i++,v++) *i += (*c) * (*v);
    }
}


RadialFunction* ContractedGaussianRF::Clone() const
{
    return new  ContractedGaussianRF(*this);
}

RadialFunction* ContractedGaussianRF::Clone(const RVec3& newCenter) const
{
    ptr_vector<RadialFunction*> newList;
    for (CITER p(itsPrimatives.begin()); p!=itsPrimatives.end(); p++) newList.push_back(p->Clone(newCenter));
    return new  ContractedGaussianRF(itsCoeffs,newList);
}

//
//  Make sure all primatives are self consistent.
//
void ContractedGaussianRF::Check() const
{
    const RVec3& center=itsPrimatives[0]->GetCenter();
    int          L     =itsPrimatives[0]->GetL     ();
    for(unsigned i=1; i<itsPrimatives.size(); i++)
    {
        if (center!=itsPrimatives[i]->GetCenter() || L!=itsPrimatives[i]->GetL())
        {
            std::cerr
                << "ContractedGaussianRF ID=" << GetID()
                << "Primatives have different L or Center" << std::endl
                << "center(0)=" << center << "  L(0)=" << L << std::endl
                << "center(" << i << ")=" << itsPrimatives[i]->GetCenter()
                << "  L(" << i << ")=" << itsPrimatives[i]->GetL() << std::endl;
        }
    }
}




