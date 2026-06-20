// File: GaussianRF.C  Implementation of the single radial type (contracted Gaussian) + the internal
// primitive helper PrimGaussian that carries the M&D integral kernels.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <memory>

module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Radial.GaussianH3;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite3;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.RNLM;

import qchem.Blaze;     // rvec_t (itsCoeff)
import qchem.Cluster;
import qchem.Math;

namespace BasisSet::Molecule::PolarizedGaussian1
{
//#######################################################################
//
//   PrimGaussian: primitive Gaussian + M&D kernels
//

PrimGaussian::PrimGaussian(double theExponent, const rvec3_t& theCenter, int theL)
    : itsExponent(theExponent)
    , itsCenter  (theCenter)
    , itsL       (theL)
    , itsH1      (0)
{
    if (itsExponent < 0)
    {
        std::cerr << "PrimGaussian exponent < 0" << std::endl;
        exit(-1);
    }
}

PrimGaussian::~PrimGaussian()
{
    delete itsH1;
}

const Hermite1& PrimGaussian::GetH1() const
{
    if (!itsH1) itsH1 = new Hermite1(itsExponent, itsL);
    return *itsH1;
}

//
//  Messy helper function for Laplacian operator.
//
static double GetGrad2(const Polarization& p1, const Polarization& p2, const GaussianCD& ab)
{
    static Polarization p0(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    const Polarization& P1 = p1;
    const Polarization& P2 = p2;

    double txx =   P1.n * P2.n * ab.H2(p0,P1-x,P2-x)
                 - 2*P1.n * ab.b * ab.H2(p0,P1-x,P2+x)
                 - 2*P2.n * ab.a * ab.H2(p0,P1+x,P2-x)
                 + 4*ab.ab        * ab.H2(p0,P1+x,P2+x);

    double tyy =   P1.l * P2.l * ab.H2(p0,P1-y,P2-y)
                 - 2*P1.l * ab.b * ab.H2(p0,P1-y,P2+y)
                 - 2*P2.l * ab.a * ab.H2(p0,P1+y,P2-y)
                 + 4*ab.ab        * ab.H2(p0,P1+y,P2+y);

    double tzz =   P1.m * P2.m * ab.H2(p0,P1-z,P2-z)
                 - 2*P1.m * ab.b * ab.H2(p0,P1-z,P2+z)
                 - 2*P2.m * ab.a * ab.H2(p0,P1+z,P2-z)
                 + 4*ab.ab        * ab.H2(p0,P1+z,P2+z);

    return txx+tyy+tzz;
}

double PrimGaussian::Integrate2C(IType type, const PrimGaussian* a, const PrimGaussian* b,
                                 const Polarization& pa, const Polarization& pb, CDCache& cache, const Cluster* cl)
{
    double s = 0.0;
    Polarization zero(0,0,0);
    const GaussianCD& ab = cache.findCD(a->GetGData(), b->GetGData());
    switch (type)
    {
        case Overlap2C :
            s = pow(Pi/ab.AlphaP,1.5)*ab.Eij*ab.H2(zero,pa,pb);
            break;
        case Repulsion2C :
            {
                auto NLMs = GaussianCD::GetNMLs(a->GetL());
                const Hermite1& H1a = a->GetH1();
                const Hermite1& H1b = b->GetH1();
                const RNLM& R = cache.find(ab);

                double factor = 1.0/(ab.ab*sqrt(ab.AlphaP));
                factor = (pb.GetTotalL()%2) ? -factor : factor;

                for (auto bNLM:NLMs)
                {
                    if (bNLM > pa) continue;
                    double ha = H1a(bNLM,pa);
                    if (ha==0.0) continue;
                    double RR = 0.0;
                    for (int n=0; n<=pb.n; n++)
                        for (int l=0; l<=pb.l; l++)
                            for (int m=0; m<=pb.m; m++)
                            {
                                Polarization NLMp(n,l,m);
                                double hb = H1b(NLMp,pb);
                                if (hb!=0.0)
                                    RR += hb*R(bNLM+NLMp);
                            }
                    s += ha*RR;
                }
                s *= 2*Pi52*factor;
            }
            break;
        case Grad2 :
            {
                double factor = pow(Pi/ab.AlphaP,1.5)*ab.Eij;
                double h = GetGrad2(pa,pb,ab);
                if (h!=0) s = factor*h;
            }
            break;
        case Nuclear :
            {
                assert(cl);
                RNLM R;
                for (auto& atom:*cl)
                    R.Add(RNLM(ab.Ltotal,ab.AlphaP,ab.P-atom->itsR), -1.0*(atom->itsZ));

                auto NLMs = GaussianCD::GetNMLs(ab.Ltotal);
                const Polarization Pab = pa + pb;
                for (auto bNLM:NLMs)
                {
                    if (bNLM > Pab) continue;
                    if (double h = ab.H2(bNLM,pa,pb); h!=0)
                        s += h*R(bNLM);
                }
                s *= 2*Pi/ab.AlphaP*ab.Eij;
            }
            break;
    }
    return s;
}

double PrimGaussian::Integrate3C(qchem::IType3C type, const PrimGaussian* ga, const PrimGaussian* gb,
                                 const Polarization& pa, const Polarization& pb, const Polarization& pc,
                                 CDCache& cache, const PrimGaussian* gc)
{
    double s = 0.0;
    switch (type)
    {
        case qchem::Overlap3C :
            {
                Hermite3* H3 = gc->GetH3(*ga,*gb);
                s = (*H3)(pa,pb,pc);
                delete H3;
            }
            break;
        case qchem::Repulsion3C :
            {
                const GaussianCD& ab(cache.findCD(ga->GetGData(), gb->GetGData()));
                const RNLM&        R(cache.find(ab.GetGData(), gc->GetGData()));

                auto NLMs = GaussianCD::GetNMLs(ab.Ltotal);
                const Hermite1& Hc = gc->GetH1();
                const Polarization Pab = pa+pb;
                for (auto nlm:NLMs)
                {
                    if (nlm > Pab) continue;
                    double hab = ab.H2(nlm,pa,pb);
                    if (hab==0.0) continue;
                    double Rs = 0.0;
                    for (int n=0; n<=pc.n; n++)
                        for (int l=0; l<=pc.l; l++)
                            for (int m=0; m<=pc.m; m++)
                            {
                                Polarization NLMp(n,l,m);
                                if (double h = Hc(NLMp,pc); h!=0.0)
                                    Rs += h*R(nlm+NLMp);
                            }
                    if (Rs!=0) s += hab*Rs;
                }
                double factor = 1.0/(ab.AlphaP*gc->GetExponent()*sqrt(ab.AlphaP+gc->GetExponent()));
                factor = (pc.GetTotalL()%2) ? -factor : factor;
                s *= 2*Pi52 * ab.Eij*factor;
            }
            break;
    }
    return s;
}

double PrimGaussian::Integrate4C(const PrimGaussian* ga, const PrimGaussian* gb,
                                 const Polarization& pa, const Polarization& pb,
                                 const Polarization& pc, const Polarization& pd,
                                 CDCache& cache, const PrimGaussian* gc, const PrimGaussian* gd)
{
    const GaussianCD& ab(cache.findCD(ga->GetGData(), gb->GetGData()));
    const GaussianCD& cd(cache.findCD(gc->GetGData(), gd->GetGData()));

    const std::vector<Polarization>& abNLMs = GaussianCD::GetNMLs(ab.Ltotal);
    const std::vector<Polarization>& cdNLMs = GaussianCD::GetNMLs(cd.Ltotal);

    double lambda = 2*Pi52/(ab.AlphaP*cd.AlphaP*sqrt(ab.AlphaP+cd.AlphaP)); //M&D 3.31
    lambda *= ab.Eij*cd.Eij; //M&D 2.25
    const RNLM& rnlm(cache.find(ab.GetGData(), cd.GetGData())); //M&D section 4A

    double s = 0.0;
    const Polarization Pab = pa + pb;
    const Polarization Pcd = pc + pd;
    for (auto abNLM:abNLMs)
    {
        if (abNLM > Pab) continue;
        double hab = ab.H2(abNLM,pa,pb);
        if (hab==0.0) continue;
        for (auto cdNLM:cdNLMs)
        {
            if (cdNLM > Pcd) continue;
            double hcd = cd.H2(cdNLM,pc,pd);
            if (hcd==0) continue;
            double r = rnlm(abNLM + cdNLM);
            if (r!=0)
                s += hab*hcd*r*cdNLM.GetSign();
        }
    }
    return s*lambda;
}

//
//  Calculate 3 center hermite functions.  Here *this is treated as the third argument, center C.
//
Hermite3* PrimGaussian::GetH3(const PrimGaussian& g1, const PrimGaussian& g2) const
{
    const double  a = g1.itsExponent, b = g2.itsExponent, c = itsExponent;
    const rvec3_t A = g1.itsCenter,   B = g2.itsCenter,   C = itsCenter;
    const int    La = g1.itsL,       Lb = g2.itsL,       Lc = itsL;

    double alphaQ = a+b+c;

    rvec3_t AB = A-B;
    rvec3_t AC = A-C;
    rvec3_t BC = B-C;
    rvec3_t Q  = (a*A+b*B+c*C)/alphaQ;

    double Eabc = pow(Pi/alphaQ,1.5)*exp( -(a*b*AB*AB + a*c*AC*AC + b*c*BC*BC) / alphaQ );

    return new GaussianH3(alphaQ,Q-A,Q-B,Q-C,La,Lb,Lc,Eabc);
}

//
// double factorial table, starts at -1, so you have to add 1 to the index.
//
static double DoubleFactData[14] = {1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080};
static inline double DoubleFact(int i) {return DoubleFactData[i+1];}

double PrimGaussian::GetNormalization(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    double s = pow(Pi/(2*itsExponent),1.5);
    double t = intpow(4*itsExponent,p.GetTotalL());
    double f = DoubleFact(2*p.n-1) * DoubleFact(2*p.l-1) * DoubleFact(2*p.m-1);
    return sqrt(t/(s*f));
}

double PrimGaussian::GetCharge(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    if ((p.n%2) || (p.l%2) || (p.m%2)) return 0;
    double s = pow(Pi/itsExponent,1.5);
    double t = intpow(2*itsExponent,p.GetTotalL()/2);
    double f = DoubleFact(p.n-1) * DoubleFact(p.l-1) * DoubleFact(p.m-1);
    return s*f/t;
}

double PrimGaussian::operator()(const rvec3_t& r) const
{
    rvec3_t dr = itsCenter-r;
    return exp(-itsExponent*dr*dr);
}

rvec3_t PrimGaussian::Gradient(const rvec3_t& r) const
{
    rvec3_t dr = itsCenter-r;
    return -2*itsExponent* (*this)(r) * dr;
}

//#######################################################################
//
//   GaussianRF: the single radial function (a contraction of >=1 primitives)
//

GaussianRF::GaussianRF()
    : itsCenter(0,0,0)
    , itsL     (0)
{}

GaussianRF::GaussianRF(double theExponent, const rvec3_t& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL)
{
    itsPrims.push_back(std::make_unique<PrimGaussian>(theExponent,theCenter,theL));
    itsCoeff.resize(1);
    itsCoeff[0]=1.0;
}

GaussianRF::GaussianRF(const rvec_t& coeffs, const rvec_t& exponents, const rvec3_t& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL)
    , itsCoeff (coeffs)
{
    assert(coeffs.size()==exponents.size());
    for (size_t i=0;i<exponents.size();++i)
        itsPrims.push_back(std::make_unique<PrimGaussian>(exponents[i],theCenter,theL));
    // Absorb each primitive's normalization into its contraction coefficient.
    for (size_t i=0;i<itsPrims.size();++i)
        itsCoeff[i] *= itsPrims[i]->GetNormalization(Polarization(theL,0,0));
}

GaussianRF::~GaussianRF() {}

// Deep-copy the primitives (a copy is a distinct object and gets its own UniqueID).
GaussianRF::GaussianRF(const GaussianRF& o)
    : itsCenter(o.itsCenter)
    , itsL     (o.itsL)
    , itsCoeff (o.itsCoeff)
{
    for (auto& g:o.itsPrims)
        itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),o.itsCenter,o.itsL));
}

GaussianRF& GaussianRF::operator=(const GaussianRF& o)
{
    if (this!=&o)
    {
        itsCenter = o.itsCenter;
        itsL      = o.itsL;
        itsCoeff  = o.itsCoeff;
        itsPrims.clear();
        for (auto& g:o.itsPrims)
            itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),o.itsCenter,o.itsL));
    }
    return *this;
}

// The same radial (exponents + contraction) placed at a new centre.
GaussianRF GaussianRF::AtCenter(const rvec3_t& newCenter) const
{
    GaussianRF ret;
    ret.itsCenter = newCenter;
    ret.itsL      = itsL;
    ret.itsCoeff  = itsCoeff;
    for (auto& g:itsPrims)
        ret.itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),newCenter,itsL));
    return ret;
}

bool GaussianRF::operator==(const GaussianRF& g) const
{
    if (norm(itsCenter-g.itsCenter) >= 0.01) return false;          // same centre (0.01 a.u.)
    if (itsPrims.size()!=g.itsPrims.size()) return false;
    for (size_t i=0;i<itsPrims.size();++i)
    {
        double ea=itsPrims[i]->GetExponent(), eb=g.itsPrims[i]->GetExponent();
        if (fabs((ea-eb)/(ea+eb)*2.0) >= 0.001) return false;       // exponents within 0.1%
        if (itsCoeff[i]!=g.itsCoeff[i]) return false;
    }
    return true;
}

double GaussianRF::GetNormalization(const Polarization& p) const
{
    // Only meaningful for an uncontracted function; contracted normalization comes from the
    // self-overlap (see PGData).
    assert(itsPrims.size()==1);
    return itsPrims[0]->GetNormalization(p);
}

double GaussianRF::GetCharge(const Polarization& p) const
{
    double ret = 0.0;
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*itsPrims[i]->GetCharge(p);
    return ret;
}

// Centre-independent identity (L + per-primitive exponent,coeff): symmetry-equivalent shells on
// different atoms share it, so the symmetry code can match them without exposing raw exponents.
std::string GaussianRF::TypeID() const
{
    std::string key = std::to_string(itsL);
    for (size_t i=0;i<itsPrims.size();++i)
        key += ":" + std::to_string(itsPrims[i]->GetExponent()) + "," + std::to_string(itsCoeff[i]);
    return key;
}

double GaussianRF::Integrate(IType type, rf_t& rb, po_t& pa, po_t& pb, CDCache& cache, const Cluster* cl) const
{
    double s = 0.0;
    for (size_t i=0;i<itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            s += itsCoeff[i]*rb.itsCoeff[j]
                 * PrimGaussian::Integrate2C(type, itsPrims[i].get(), rb.itsPrims[j].get(), pa, pb, cache, cl);
    return s;
}

// this is centre C: <ab|c>
double GaussianRF::Integrate(qchem::IType3C type, rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc, CDCache& cache) const
{
    double s = 0.0;
    for (size_t i=0;i<ra.itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            for (size_t k=0;k<itsPrims.size();++k)
                s += ra.itsCoeff[i]*rb.itsCoeff[j]*itsCoeff[k]
                     * PrimGaussian::Integrate3C(type, ra.itsPrims[i].get(), rb.itsPrims[j].get(),
                                                 pa, pb, pc, cache, itsPrims[k].get());
    return s;
}

// this is centre D: (ab|cd)
double GaussianRF::Integrate(rf_t& ra, rf_t& rb, rf_t& rc, po_t& pa, po_t& pb, po_t& pc, po_t& pd, CDCache& cache) const
{
    double s = 0.0;
    for (size_t i=0;i<ra.itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            for (size_t k=0;k<rc.itsPrims.size();++k)
                for (size_t l=0;l<itsPrims.size();++l)
                    s += ra.itsCoeff[i]*rb.itsCoeff[j]*rc.itsCoeff[k]*itsCoeff[l]
                         * PrimGaussian::Integrate4C(ra.itsPrims[i].get(), rb.itsPrims[j].get(),
                                                     pa, pb, pc, pd, cache,
                                                     rc.itsPrims[k].get(), itsPrims[l].get());
    return s;
}

std::ostream& GaussianRF::Write(std::ostream& os) const
{
    if (itsPrims.size()==1)
        os << "Primative  " << std::setw(8) << itsPrims[0]->GetExponent();
    else
    {
        os << "Contracted {";
        for (auto& g:itsPrims) os << g->GetExponent() << " ";
        os << "}";
    }
    return os;
}

double GaussianRF::operator()(const rvec3_t& r) const
{
    double ret = 0.0;
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*(*itsPrims[i])(r);
    return ret;
}

rvec3_t GaussianRF::Gradient(const rvec3_t& r) const
{
    rvec3_t ret(0,0,0);
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*itsPrims[i]->Gradient(r);
    return ret;
}

} //namespace BasisSet::Molecule::PolarizedGaussian1
