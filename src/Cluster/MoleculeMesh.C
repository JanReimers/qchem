// File: MoleculeMesh.C  mesh implementation



#include "Imp/Cluster/MoleculeMesh.H"
#include "Imp/Cluster/Atom.H"
#include <Cluster.H>
#include <MeshParams.H>
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/io3d.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

double          Poly             (double,int m_mu);

//
//  Use Becke's fuzzy polyedra algorithm for integrating over a molecule.
//  See A. D. Becke, J. Chem. Phys, 88(4), page 2547 (1988).
//
MoleculeMesh::MoleculeMesh(const Cluster& cl, const MeshParams& mp)
{
    assert(mp.m_mu>=0);

    int natom=cl.GetNumAtoms();
    size_t ia=1;
    for (auto a:cl) 
    {
        Mesh* mesh_a= a->CreateMesh(mp);
        for (auto rw: *mesh_a)
        {
            RVec3 r=::r(rw);
            // Load up up matrix of cutoff profiles.
            Matrix<double> s(natom,natom);
            size_t ib=1;
            for (auto b:cl)
            {
                double Rb=norm(r-b->itsR);
                size_t ic=0;
                for (auto c:cl)
                {
                    ic++;
                    if (ib==ic) continue;
                    double Rc=norm(r-c->itsR);
                    double ubc=(Rb-Rc)/norm(b->itsR-c->itsR);
                    s(ib,ic) = Poly(ubc,mp.m_mu);
                }
                ib++;
            }
            //cout << "s=" << s << endl;
            //  Load up and array of cell functions
            Vector<double> P(natom);
            Fill(P,1.0);
            for (int i=1; i<=natom; i++)
                for (int j=1; j<=natom; j++)
                    if (i!=j) P(i)*=s(i,j);
            
            //std::cout << "r,w = "<< r << "," << w(rw) <<    " P=" << P(1) << " " << P(2) << std::endl;
            //std::cout << "r,w = "<< r << "," << w(rw) << std::endl;

            if(natom>1 && P(ia)>0)
            {
                double relativeWeight=P(ia)/Sum(P);
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
#include "oml/io3d.h"

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

