// File: PlaneWaveBF.C  Polarized Gaussian in 3D space.



#include "BasisSetImplementation/PlaneWave/PlaneWaveBF.H"
#include "Mesh/MeshBrowser.H"
#include "Misc/DFTDefines.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include "oml/vector.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

//#######################################################################
//
//  Polarized gaussian implementation
//

PlaneWaveBF::PlaneWaveBF()
    : itsG            (0,0,0)
    , itsNormalization(0)
{};

PlaneWaveBF::PlaneWaveBF(const RVec3& G, double Volume)
    : itsG            (G )
    , itsNormalization(1/sqrt(Volume))
{};

void PlaneWaveBF::Init(double norm, double charge)
{
    assert(!std::isnan(norm));
    assert(norm>0);
    assert(charge>0);

    itsNormalization=norm;
    itsCharge=charge;
}

bool PlaneWaveBF::operator==(const BasisFunction& bf) const
{
    const PlaneWaveBF* pgbf = dynamic_cast<const PlaneWaveBF*>(&bf);
    bool ret=pgbf!=0;
    if (ret) ret=ret&&(itsG==pgbf->itsG);
    return ret;
}

double PlaneWaveBF::GetNormalization() const
{
    return itsNormalization;
}

double PlaneWaveBF::GetCharge() const
{
    return itsG==RVec3(0,0,0) ? 1/itsNormalization : 0;
}

//------------------------------------------------------------------------
//
//  These basis functions are really just proxys, they dont own the radial.
//  So no information needs be pickled.
//
std::ostream& PlaneWaveBF::Write(std::ostream& os) const
{
    UniqueID::Write(os);

    if (Pretty())
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os <<  std::setw(5) << std::setprecision(2) << itsG;
    }

    if (Binary())
    {
        os << itsG;
        BinaryWrite(itsNormalization,os);
    }

    if (Ascii()) os << itsG << " " << itsNormalization << " ";

    return os;
}

std::istream& PlaneWaveBF::Read(std::istream& is)
{
    UniqueID::Read(is);

    if (Binary())
    {
        is >> itsG;
        BinaryRead(itsNormalization,is);
    }

    if (Ascii()) is >> itsG >> itsNormalization;

    return is;
}

std::complex<double> PlaneWaveBF::operator()(const RVec3& r) const
{
    return itsNormalization*exp(std::complex<double>(0,2*Pi*itsG*r));
}

PlaneWaveBF::Vec3 PlaneWaveBF::Gradient(const RVec3& r) const
{
    return itsNormalization*exp(std::complex<double>(0,2*Pi*itsG*r))*2.0*Pi*
           Vec3(std::complex<double>(0,itsG.x),std::complex<double>(0,itsG.y),std::complex<double>(0,itsG.z));
}

void PlaneWaveBF::Eval(const Mesh& mesh, Vec& vec) const
{
    Vec::iterator i(vec.begin());
    MeshBrowser            m(mesh);
    for(; i!=vec.end()&&m; i++,m++) *i += itsNormalization*exp(std::complex<double>(0,2*Pi*itsG*m.R()));
}

void PlaneWaveBF::EvalGrad(const Mesh& mesh, Vec3Vec& vec) const
{
    Vec3Vec::iterator i(vec.begin());
    MeshBrowser                m(mesh);
    for (; i!=vec.end()&&m; i++,m++)
        *i += itsNormalization*exp(std::complex<double>(0,2*Pi*itsG*m.R()))*2.0*Pi*
           Vec3(std::complex<double>(0,itsG.x),std::complex<double>(0,itsG.y),std::complex<double>(0,itsG.z));
}

BasisFunction* PlaneWaveBF::Clone() const
{
    return new  PlaneWaveBF(*this);
}

