// File: PolarizedGaussianBF.C  Polarized Gaussian in 3D space.



#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBF.H"
#include "Mesh/MeshBrowser.H"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

//#######################################################################
//
//  Polarized gaussian implementation
//

PolarizedGaussianBF::PolarizedGaussianBF()
    : itsRadial       (0)
    , itsNormalization(0)
    , itsCharge       (0)
{};

PolarizedGaussianBF::PolarizedGaussianBF(const RadialFunction* theRF,const Polarization& thePol)
    : itsRadial(theRF )
    , itsPol   (thePol)
    , itsNormalization(0)
    , itsCharge       (0)
{
    assert(itsRadial);
};

void PolarizedGaussianBF::Init(double norm, double charge)
{
    assert(!std::isnan(norm));
    assert(norm>0);

    itsNormalization=norm;
    itsCharge=charge;
}

bool PolarizedGaussianBF::operator==(const BasisFunction& bf) const
{
    assert(itsRadial);
    const PolarizedGaussianBF& pgbf = dynamic_cast<const PolarizedGaussianBF&>(bf);
    assert(&pgbf);
    return (*itsRadial)==(*pgbf.itsRadial) && (itsPol==pgbf.itsPol);
}

double PolarizedGaussianBF::GetNormalization() const
{
    return itsNormalization;
}

double PolarizedGaussianBF::GetCharge() const
{
    return itsCharge;
}

//------------------------------------------------------------------------
//
//  These basis functions are really just proxys, they dont own the radial.
//  So no information needs be pickled.
//
std::ostream& PolarizedGaussianBF::Write(std::ostream& os) const
{
    assert(itsRadial);
    UniqueID::Write(os);

    if (StreamableObject::Pretty())
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os <<  std::setw(5) << std::setprecision(2) << itsRadial->GetCenter() << " "
        << itsPol << " " << *itsRadial;
    }
    return os;
}

std::istream& PolarizedGaussianBF::Read(std::istream& is)
{
    UniqueID::Read(is);
    return is;
}

void PolarizedGaussianBF::Insert(const RadialFunction* theRF,const Polarization& thePol)
{
    itsRadial       =theRF;
    itsPol          =thePol;
    assert(itsRadial);
}

double PolarizedGaussianBF::operator()(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*itsPol(dr) * rf(r);
}

PolarizedGaussianBF::RVec3 PolarizedGaussianBF::Gradient(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*(itsPol.Gradient(dr) * rf(r) + itsPol(dr) * rf.Gradient(r));
}

void PolarizedGaussianBF::Eval(const Mesh& mesh, Vec& vec) const
{
    assert(itsRadial);
    Vector<double> radial((*itsRadial)(mesh));
    Vector<double>::const_iterator  rf(radial.begin());
    Vec            ::iterator i(vec.begin());
    MeshBrowser               m(mesh);
    for(; i!=vec.end()&&m&&rf!=radial.end(); i++,m++,rf++)
    {
//        std::cout << itsNormalization << " " << *rf << " " << itsPol( m.R()-itsRadial->GetCenter() ) << std::endl;
        *i += itsNormalization * (*rf) * itsPol( m.R()-itsRadial->GetCenter() );
    }
}

void PolarizedGaussianBF::EvalGrad(const Mesh& mesh, Vec3Vec& vec) const
{
    assert(itsRadial);
    Vector<double> radial((*itsRadial)(mesh));
    Vector<double>::const_iterator  rf(radial.begin());
    Vector<RVec3 > gradrf(itsRadial->Gradient(mesh));
    Vector<RVec3 >::const_iterator grf(gradrf.begin());
    Vec3Vec        ::iterator  i(vec.begin());
    MeshBrowser                m(mesh);

    for (; rf!=radial.end()&&grf!=gradrf.end()&&i!=vec.end()&&m; rf++,grf++,i++,m++)
    {
        RVec3 dr = m.R()-itsRadial->GetCenter();
        *i += itsNormalization*(itsPol(dr) * (*grf) + itsPol.Gradient(dr) * (*rf));
    }
}

BasisFunction* PolarizedGaussianBF::Clone() const
{
    assert(itsRadial);
    return new  PolarizedGaussianBF(*this);
}

