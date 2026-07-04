// File: ElectronConfigurations/Imp/AtomDirac_EC.C  Electron configuration for atoms.
module;
#include <cassert>
#include <cstdlib>
#include <algorithm>

module qchem.ElectronConfiguration.AtomDirac;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Atom.Spherical;
import qchem.Symmetry.Factory;

namespace qchem {


AtomDirac_EC::AtomDirac_EC(int Z) : Atom_EC(Z, NsOnly_t{})
{
    double s=0.5; //electron spin.
    for (size_t l:iv_t(0,itsLMax+1))
    {
        int Nf=itsNs.Nf[l],Nv=itsNs.Nv[l],Nu=abs(itsNs.Nu[l]);
        size_t gl=2*l+1; //l degeneracy for one spin.
        int NCoreLevels=Nf/(2*gl);
        { // fill the j=l-1/2 levels first
            double j=Symmetry::Atom::SphericalSpinor::j(l,-s);
            int    κ=Symmetry::Atom::SphericalSpinor::κ(l,-s);
            if (κ!=0)
            {
                int g=(2*j+1)/2; //degeneracy for one spin state
                assert(Nf%2==0);
                assert((Nv-Nu)%2==0);
                int Npair=std::min(g,(Nv-Nu)/2);
                Nu=std::min(Nu,g-Npair);
                Nf=std::min(Nf,2*g);
                if (Nu==0 || Nu==g) //No mj splitting
                {
                    sym_t s=Symmetry::ΩFactory(κ);
                    itsOccupations[Irrep(Spin::Up  ,s)]=NCoreLevels*g+Npair+Nu;
                    itsOccupations[Irrep(Spin::Down,s)]=NCoreLevels*g+Npair;
                    itsNs.Nf[l]-=2*NCoreLevels*g;
                    itsNs.Nv[l]-=2*Npair+Nu;
                    itsNs.Nu[l]-=Nu;
                }
                else // mj splitting
                {
                    assert(Nf/2%g==0);
                    int gu=Nu,gp=g-gu;
                    assert(gu+gp<=g);
                    rvec_t mj_p(gp),mj_u(gu);
                    double mj=-j;
                    for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                    for (size_t i=0;i<gp;i++) mj_p[i]=mj++;
                    sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                    sym_t su=Symmetry::ΩFactory(κ,mj_u);
                    SetSplitOccupations(sp,su,NCoreLevels,gp,gu,Npair,Nu);
                    itsNs.Nf[l]-=2*NCoreLevels*g;
                    itsNs.Nv[l]-=2*Npair;
                    itsNs.Nu[l]-=Nu;
                }
            }
        }
        Nf=itsNs.Nf[l];Nv=itsNs.Nv[l];Nu=abs(itsNs.Nu[l]);
        { //Then fill the j=l+1/2 levels with any left over electrons.
            double j=Symmetry::Atom::SphericalSpinor::j(l,+s);
            int    κ=Symmetry::Atom::SphericalSpinor::κ(l,+s);
            int g=(2*j+1)/2; //degeneracy for one spin state
            assert(Nf%2==0);
            assert(Nu<=g);
            if (Nu==0 || Nu==g) //No mj splitting
            {
                sym_t s=Symmetry::ΩFactory(κ);
                itsOccupations[Irrep(Spin::Up  ,s)]=NCoreLevels*g+Nu;
                itsOccupations[Irrep(Spin::Down,s)]=NCoreLevels*g;
                itsNs.Nf[l]-=2*NCoreLevels*g;
                itsNs.Nv[l]-=Nu;
                itsNs.Nu[l]-=Nu;
            }
            else // mj splitting — j=l+1/2 Kramers structure differs from NR
            {
                assert(Nf/2%g==0);
                assert(Nv<2*g);
                int Npair=(Nv-Nu)/2;
                int gu=Nu;
                int gp=NCoreLevels>0 ? 2*g-gu : 2*Npair;
                rvec_t mj_p(gp),mj_u(gu);
                double mj=-j;
                for (size_t i=0;i<gp;i++) mj_p[i]=mj++;
                for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                sym_t su=Symmetry::ΩFactory(κ,mj_u);
                gp/=2;
                itsOccupations[Irrep(Spin::Up  ,sp)]=NCoreLevels*gp + Npair*(1 + 2*NCoreLevels);
                itsOccupations[Irrep(Spin::Down,sp)]=NCoreLevels*(gp+Nu) + Npair*(1 - NCoreLevels);
                itsOccupations[Irrep(Spin::Up  ,su)]=NCoreLevels*gu + Nu - 2*NCoreLevels*Npair;
                itsOccupations[Irrep(Spin::Down,su)]=NCoreLevels*Npair;
                itsNs.Nf[l]-=2*(NCoreLevels*g);
                itsNs.Nv[l]-=2*Npair+Nu;
                itsNs.Nu[l]-=Nu;
            }
        }
    }
}


   

} // namespace qchem