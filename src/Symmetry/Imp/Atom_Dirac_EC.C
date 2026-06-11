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
    double s=0.5; //electron spin.
    for (size_t l:iv_t(0,itsLMax+1))
    {
        int Nf=itsNs.Nf[l],Nv=itsNs.Nv[l],Nu=abs(itsNs.Nu[l]);
        size_t gl=2*l+1; //l degeneracy for one spin.
        int NCoreLevels=Nf/(2*gl);
        { // fill the j=l-1/2 levels first
            double j=Symmetry::SphericalSpinor::j(l,-s);
            int    κ=Symmetry::SphericalSpinor::κ(l,-s);
            if (κ!=0)
            {
                int g=(2*j+1)/2; //degeneracy for one spin state
                assert(Nf%2==0); //Full shells better be even!
                //assert(Nu<=g); //No the unpaired electrons could, for example be spread over p1/2 and p3/2
                assert((Nv-Nu)%2==0);
                int Npair=std::min(g,(Nv-Nu)/2);
                Nu=std::min(Nu,g-Npair); //All unpaired may not fit in this irrep
                Nf=std::min(Nf,2*g);
                // Npair=std::min(g-Nu,Npair); //All paired may also not fit in this irrep
                if (Nu==0 || Nu==g) //No m splitting, but g unapired electrons, or all paired
                {
                    sym_t s=Symmetry::ΩFactory(κ);
                    itsOccupations[Irrep(Spin::Up  ,s)]=NCoreLevels*g+Npair+Nu;
                    itsOccupations[Irrep(Spin::Down,s)]=NCoreLevels*g+Npair;
                    itsNs.Nf[l]-=2*NCoreLevels*g;
                    itsNs.Nv[l]-=2*Npair+Nu;
                    itsNs.Nu[l]-=Nu;
                }
                else // M splitting.  All shells get the same M splitting.
                {
                    assert(Nf/2%g==0);
                    // assert(Nv<2*g); //Again no, the unpaired electrons could, for example be spread over p1/2 and p3/2
                    int gu=Nu,gp=g-gu; //unpaired and paired degeneracies
                    assert(gu+gp<=g);
                    rvec_t mj_p(gp),mj_u(gu);
                    double mj=-j;
                    for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                    for (size_t i=0;i<gp;i++) mj_p[i]=mj++;
                    sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                    sym_t su=Symmetry::ΩFactory(κ,mj_u);
                    itsOccupations[Irrep(Spin::Up  ,sp)]=NCoreLevels*gp+Npair;
                    itsOccupations[Irrep(Spin::Down,sp)]=NCoreLevels*gp+Npair;
                    itsOccupations[Irrep(Spin::Up  ,su)]=NCoreLevels*gu+Nu;
                    itsOccupations[Irrep(Spin::Down,su)]=NCoreLevels*gu;
                    itsNs.Nf[l]-=2*NCoreLevels*g;
                    itsNs.Nv[l]-=2*Npair;
                    itsNs.Nu[l]-=Nu;
                }
                Display();
            }
        }
        Nf=itsNs.Nf[l];Nv=itsNs.Nv[l];Nu=abs(itsNs.Nu[l]);
        { //Then fill the j=l+1/2 levels with any left over electrons.
            
            double j=Symmetry::SphericalSpinor::j(l,+s);
            int    κ=Symmetry::SphericalSpinor::κ(l,+s);
            
            int g=(2*j+1)/2; //degeneracy for one spin state
            assert(Nf%2==0); //Full shells better be even!
            assert(Nu<=g);
            if (Nu==0 || Nu==g) //No m splitting, but g unapired electrons, or all paired
            {
                sym_t s=Symmetry::ΩFactory(κ);
                itsOccupations[Irrep(Spin::Up  ,s)]=NCoreLevels*g+Nu;
                itsOccupations[Irrep(Spin::Down,s)]=NCoreLevels*g;
                itsNs.Nf[l]-=2*NCoreLevels*g;
                itsNs.Nv[l]-=Nu;
                itsNs.Nu[l]-=Nu;
            }
            else // M splitting.  All shells get the same M splitting.
            {
                assert(Nf/2%g==0);
                assert(Nv<2*g);
                int Npair=(Nv-Nu)/2; //Valance pairs.
                int gu=Nu; //# mj states for the un paired sector.
                int gp=NCoreLevels>0 ? 2*g-gu : 2*Npair; //number of mj states for the paired sector.
                rvec_t mj_p(gp),mj_u(gu);
                double mj=-j;
                for (size_t i=0;i<gp;i++) mj_p[i]=mj++; //convention, fill the paired electrons first.
                for (size_t i=0;i<gu;i++) mj_u[i]=mj++;
                // assert(mj==j+1);
                sym_t sp=Symmetry::ΩFactory(κ,mj_p);
                sym_t su=Symmetry::ΩFactory(κ,mj_u);
                cout << "  paired=" << *sp << endl;
                cout << "unpaired=" << *su << endl;
                // Irrep isp=Irrep(Spin::Up  ,sp);
                // Irrep isu=Irrep(Spin::Up  ,su);
                gp/=2; //Single spin version of degeneracy.
                itsOccupations[Irrep(Spin::Up  ,sp)]=NCoreLevels*gp + Npair*(1 + 2*NCoreLevels);
                itsOccupations[Irrep(Spin::Down,sp)]=NCoreLevels*(gp+Nu) + Npair*(1 - NCoreLevels);
                itsOccupations[Irrep(Spin::Up  ,su)]=NCoreLevels*gu + Nu - 2*NCoreLevels*Npair;
                itsOccupations[Irrep(Spin::Down,su)]=NCoreLevels*Npair;
                itsNs.Nf[l]-=2*(NCoreLevels*g);
                itsNs.Nv[l]-=2*Npair+Nu;
                itsNs.Nu[l]-=Nu;
                // cout << "isp=" << isp << " seqn=" 
                // << isp.SequenceIndex() 
                // << " occ=" << NCoreLevels*gp+Npair << endl;
                // cout << "isu=" << isu << " seqn=" 
                // << isu.SequenceIndex() 
                // << " occ=" << NCoreLevels*gp+Nu << endl;
            }
            Display();
        }
    }

    
}


   
