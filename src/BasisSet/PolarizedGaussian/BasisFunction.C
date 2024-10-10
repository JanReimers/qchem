// File: BasisFunction.C  Polarized Gaussian in 3D space.



#include "Imp/BasisSet/PolarizedGaussian/BasisFunction.H"
#include "Mesh/MeshBrowser.H"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

namespace PolarizedGaussian
{

//#######################################################################
//
//  Polarized gaussian implementation
//

BasisFunction::BasisFunction()
    : itsRadial       (0)
    , itsNormalization(0)
{};

BasisFunction::BasisFunction(const RadialFunction* theRF,const Polarization& thePol, double norm)
    : itsRadial(theRF )
    , itsPol   (thePol)
    , itsNormalization(norm)
{
    assert(itsRadial);
};

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    assert(itsRadial);
    const BasisFunction& pgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&pgbf);
    return (*itsRadial)==(*pgbf.itsRadial) && (itsPol==pgbf.itsPol);
}

//------------------------------------------------------------------------
//
//  These basis functions are really just proxys, they dont own the radial.
//  So no information needs be pickled.
//
std::ostream& BasisFunction::Write(std::ostream& os) const
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

std::istream& BasisFunction::Read(std::istream& is)
{
    UniqueID::Read(is);
    return is;
}

void BasisFunction::Insert(const RadialFunction* theRF,const Polarization& thePol)
{
    itsRadial       =theRF;
    itsPol          =thePol;
    assert(itsRadial);
}

double BasisFunction::operator()(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*itsPol(dr) * rf(r);
}

BasisFunction::RVec3 BasisFunction::Gradient(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*(itsPol.Gradient(dr) * rf(r) + itsPol(dr) * rf.Gradient(r));
}

void BasisFunction::Eval(const Mesh& mesh, Vec& vec) const
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

void BasisFunction::EvalGrad(const Mesh& mesh, Vec3Vec& vec) const
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

BasisFunction* BasisFunction::Clone() const
{
    assert(itsRadial);
    return new  BasisFunction(*this);
}

} //namespace PolarizedGaussian
