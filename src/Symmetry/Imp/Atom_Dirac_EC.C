// File: Symmetry/Imp/Atom_Dirac_EC.C  Electron configuration for atoms.
module;
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <vector>

module qchem.Symmetry.Atom_Dirac_EC;
import Common.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Ylm;
import qchem.stl_io;

using std::cout;
using std::endl;


Atom_Dirac_EC::Atom_Dirac_EC(int Z)
: itsLMax(0), itsLValance(0)
{
    assert(Z>0);
    assert(Z<=N_Elements);

    size_t l=0;
    { //l=0
        int N=itsNs.N[l],Nf=itsNs.Nf[l],Nv=itsNs.Nv[l],Nu=abs(itsNs.Nu[l]);
        int κ=-1;
        double j=Omega_k_Sym::j(κ);
        size_t g=2*j+1;
        int g2=2*g;
        assert(Nv<g2);
        assert(Nu<=g);
        assert(Nf%g2==0);
        int Nlevel=Nf/g2;
        assert((Nv-Nu)%2==0);
        int Npair=(Nv-Nu)/2;
        if (Nu==0) //No unpaired electrons
        {
            assert(N%2==0);
            sym_t Ω(new Omega_k_Sym(κ));
            itsOccupations[Irrep_QNs(Spin::Up  ,Ω)]=N/2;
            itsOccupations[Irrep_QNs(Spin::Down,Ω)]=N/2;
        }
        else if (Nu==g) //exactly half filled all up.
        {
            sym_t Ω(new Omega_k_Sym(κ));
            itsOccupations[Irrep_QNs(Spin::Up  ,Ω)]=N;
            itsOccupations[Irrep_QNs(Spin::Down,Ω)]=N-Nu;
        }
        else //Partially filled, mj splitting
        {
            std::vector<double> mj_p,mj_u;
            double mj=-j;
            assert((Nv+Nu)%2==0);
            size_t nup=(Nv+Nu)/2;
            size_t ndn=(Nv-Nu)/2;
            int gu=Nu,gp=g-gu; //unpaired and paired degeneracies
            for (size_t i=0;i<gu;i++) mj_u.push_back(mj++);
            for (size_t i=0;i<gp;i++) mj_p.push_back(mj++);
            assert(mj==j+1);
            sym_t Ω_p(new Omega_kmj_Sym(κ,mj_p));
            sym_t Ω_u(new Omega_kmj_Sym(κ,mj_u));
            itsOccupations[Irrep_QNs(Spin::Up  ,Ω_p)]=Nlevel*gp+Npair;
            itsOccupations[Irrep_QNs(Spin::Down,Ω_p)]=Nlevel*gp+Npair;
            itsOccupations[Irrep_QNs(Spin::Up  ,Ω_u)]=Nlevel*gu+Nu;
            itsOccupations[Irrep_QNs(Spin::Down,Ω_u)]=Nlevel*gu;
        }
       
        
    }
    
}


int Atom_Dirac_EC::GetN() const
{
    int ne=0;
    for (auto n:itsNs.N) ne+=n; //Sum over l
    return ne;
}

   
