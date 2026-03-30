module;
#include <cassert>
#include <iostream>
#include <blaze/Math.h>
module qchem.Cluster.MoleculeMesh;
import qchem.Atom;

double          Poly             (double,int m_mu);

//
//  Use Becke's fuzzy polyedra algorithm for integrating over a molecule.
//  See A. D. Becke, J. Chem. Phys, 88(4), page 2547 (1988).
//
MoleculeMesh::MoleculeMesh(const Cluster& cl, const MeshParams& mp)
{
    assert(mp.m_mu>=0);

    int natom=cl.GetNumAtoms();
    size_t ia=0;
    for (auto& a:cl) 
    {
        Mesh* mesh_a= a->CreateMesh(mp);
        for (auto rw: *mesh_a)
        {
            rvec3_t r=::r(rw);
            // Load up up matrix of cutoff profiles.
            rmat_t s(natom,natom);
            size_t ib=0;
            for (auto& b:cl)
            {
                double Rb=norm(r-b->itsR);
                size_t ic=0;
                for (auto& c:cl)
                {
                    
                    if (ib!=ic) 
                    {
                        double Rc=norm(r-c->itsR);
                        double ubc=(Rb-Rc)/norm(b->itsR-c->itsR);
                        s(ib,ic) = Poly(ubc,mp.m_mu);
                    }
                    ic++;
                }
                ib++;
            }
            //cout << "s=" << s << endl;
            //  Load up and array of cell functions
            rvec_t P(natom,1.0);
            for (int i=0; i<natom; i++)
                for (int j=0; j<natom; j++)
                    if (i!=j) P[i]*=s(i,j);
            
            //std::cout << "r,w = "<< r << "," << w(rw) <<    " P=" << P(1) << " " << P(2) << std::endl;
            //std::cout << "r,w = "<< r << "," << w(rw) << std::endl;

            if(natom>1 && P[ia]>0)
            {
                double relativeWeight=P[ia]/blaze::sum(P);
//                cout << "ia,r,w=" << ia << " " << r << " " << relativeWeight << endl;
                push_back(::r(rw),::w(rw)*relativeWeight);
            }
            else if(natom==1)
                push_back(::r(rw),::w(rw));
        }
        delete mesh_a;
        ia++;
    }
}


Mesh* MoleculeMesh::Clone() const
{
    return new MoleculeMesh(*this);
}

using std::cout;
using std::endl;

double Poly(double u,int k)
{
//    assert(u<= 1.000000000001);
//    assert(u>=-1.000000000001);
    assert(k< 10);
    assert(k>=0);
    for (int i=k; i>=0; i--) 
        u=0.5*(3*u-u*u*u);
    return 0.5*(1-u);
}

