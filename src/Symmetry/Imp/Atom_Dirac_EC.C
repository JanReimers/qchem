// File: Symmetry/Imp/Atom_Dirac_EC.C  Electron configuration for atoms.
module;
#include <cassert>
#include <iostream>

module qchem.Symmetry.Atom_Dirac_EC;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
using std::cout;
using std::endl;


Atom_Dirac_EC::Atom_Dirac_EC(int Z) : Atom_EC(Z)
{
    itsOccupations.clear();
    for (size_t l:iv_t(0,itsLMax+1))
    {
        { // fill the j=l-1/2 levels first
            double ms=-0.5;
            double j=Symmetry::SphericalSpinor::j(l,ms);
            int    κ=Symmetry::SphericalSpinor::κ(l,ms);
            if (κ!=0)
            {
                Display();
                int Nf=itsNs.Nf[l],Nv=itsNs.Nv[l],Nu=abs(itsNs.Nu[l]);
                int g=(2*j+1)/2; //degeneracy for one spin state
                assert(Nf%2==0); //Full shells better be even!
                //assert(Nu<=g); //No the unpaired electrons could, for example be spread over p1/2 and p3/2
                assert((Nv-Nu)%2==0);
                int Npair=(Nv-Nu)/2;
                Nu=std::min(Nu,g); //All unpaired may not fit in this irrep
                Npair=std::min(g-Nu,Npair); //All paired may also not fit in this irrep
                if (Nu==0 || Nu==g) //No m splitting, but g unapired electrons, or all paired
                {
                    sym_t s=Symmetry::ΩFactory(κ);
                    itsOccupations[Irrep(Spin::Up  ,s)]=Nf/2+Nu;
                    itsOccupations[Irrep(Spin::Down,s)]=Nf/2;
                    itsNs.Nf[l]-=Nf;
                    itsNs.Nv[l]-=Nu;
                    itsNs.Nu[l]-=Nu;
                }
                else // M splitting.  All shells get the same M splitting.
                {
                    assert(Nf/2%g==0);
                    // assert(Nv<2*g); //Again no, the unpaired electrons could, for example be spread over p1/2 and p3/2
                    int Nlevel=Nf/2/g;
                    int gu=Nu,gp=g-gu; //unpaired and paired degeneracies
                    assert(gu+gp<=g);
                    rvec_t mj_p(gp),mj_u(gu);
                    double mj=-j;
                    for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                    for (size_t i=0;i<gp;i++) mj_p[i]=mj++;
                    sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                    sym_t su=Symmetry::ΩFactory(κ,mj_u);
                    itsOccupations[Irrep(Spin::Up  ,sp)]=Nlevel*gp+Npair;
                    itsOccupations[Irrep(Spin::Down,sp)]=Nlevel*gp+Npair;
                    itsOccupations[Irrep(Spin::Up  ,su)]=Nlevel*gu+Nu;
                    itsOccupations[Irrep(Spin::Down,su)]=Nlevel*gu;
                    itsNs.Nf[l]-=2*(Nlevel*g);
                    itsNs.Nv[l]-=2*Npair;
                    itsNs.Nu[l]-=Nu;
                }
                Display();
            }
        }
        { //Then fill the j=l+1/2 levels with any left over electrons.
            double ms=0.5;
            double j=Symmetry::SphericalSpinor::j(l,ms);
            int    κ=Symmetry::SphericalSpinor::κ(l,ms);
            
            int Nf=itsNs.Nf[l],Nv=itsNs.Nv[l],Nu=abs(itsNs.Nu[l]);
            int g=(2*j+1)/2; //degeneracy for one spin state
            assert(Nf%2==0); //Full shells better be even!
            assert(Nu<=g);
            assert((Nv-Nu)%2==0);
            if (Nu==0 || Nu==g) //No m splitting, but g unapired electrons, or all paired
            {
                sym_t s=Symmetry::ΩFactory(κ);
                itsOccupations[Irrep(Spin::Up  ,s)]=Nf/2+Nu;
                itsOccupations[Irrep(Spin::Down,s)]=Nf/2;
                itsNs.Nf[l]-=Nf;
                itsNs.Nv[l]-=Nu;
                itsNs.Nu[l]-=Nu;
            }
            else // M splitting.  All shells get the same M splitting.
            {
                assert(Nf/2%g==0);
                assert(Nv<2*g);
                int Nlevel=Nf/2/g;
                int Npair=(Nv-Nu)/2;
                int gu=Nu,gp=g-gu; //unpaired and paired degeneracies
                rvec_t mj_p(gp),mj_u(gu);
                double mj=-j;
                for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                for (size_t i=0;i<gp;i++) mj_p[i]=mj++;
                // assert(mj==j+1);
                sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                sym_t su=Symmetry::ΩFactory(κ,mj_u);
                // Irrep isp=Irrep(Spin::Up  ,sp);
                // Irrep isu=Irrep(Spin::Up  ,su);
                itsOccupations[Irrep(Spin::Up  ,sp)]=Nlevel*gp+Npair;
                itsOccupations[Irrep(Spin::Down,sp)]=Nlevel*gp+Npair;
                itsOccupations[Irrep(Spin::Up  ,su)]=Nlevel*gu+Nu;
                itsOccupations[Irrep(Spin::Down,su)]=Nlevel*gu;
                itsNs.Nf[l]-=2*(Nlevel*g);
                itsNs.Nv[l]-=2*Npair;
                itsNs.Nu[l]-=Nu;
                // cout << "isp=" << isp << " seqn=" 
                // << isp.SequenceIndex() 
                // << " occ=" << Nlevel*gp+Npair << endl;
                // cout << "isu=" << isu << " seqn=" 
                // << isu.SequenceIndex() 
                // << " occ=" << Nlevel*gp+Nu << endl;
            }
        }
    }

    
}


   
