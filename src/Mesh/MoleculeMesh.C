// File: MoleculeMesh.C  mesh implementation



#include "Mesh/MoleculeMesh.H"
#include "Mesh/MeshBrowser.H"
#include "Cluster/Cluster.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

Vector<RVec3>  GetPositions     (const Cluster&);
void            LoadFuzzyPoints  (const Atom&, const Cluster&, int, Vector<RVec3>::iterator&, Vector<double>::iterator&,int&);
void            GetCutoffProfiles(Matrix<double>&,const Vector<RVec3>&, const RVec3&, int);
double          CutoffProfile    (const RVec3&,const RVec3&,const RVec3&,int);
double          Poly             (double,int);
void            CalcCellFunctions(Vector<double>&, const Matrix<double>&);

//
//  Use Becke's fuzzy polyedra algorithm for integrating over a molecule.
//  See A. D. Becke, J. Chem. Phys, 88(4), page 2547 (1988).
//
MoleculeMesh::MoleculeMesh(const Cluster& cl, int m)
    : MeshImplementation()
{
//
//  First count total number of mesh points for all atoms.
//
    int nmax=0;
    for (auto atom:cl) nmax+=atom->GetIntegrationMesh()->GetNumPoints();
    std::cout << "Molecular mesh: total points=" << nmax;

    Vector<RVec3 > Points (nmax);
    Vector<double> Weights(nmax);
//
//  Now generate fuzzy points for each atom.
//
    Vector<RVec3 >::iterator ip(Points.begin());
    Vector<double>::iterator iw(Weights.begin());
    int numpoints=0;
    for (auto atom:cl) LoadFuzzyPoints(*atom,cl,m,ip,iw,numpoints);

    std::cout << ", Fuzzy points=" << numpoints << std::endl;
    Points .SetLimits(VecLimits(1,numpoints),true);
    Weights.SetLimits(VecLimits(1,numpoints),true);

    Initialize(Points,Weights);
}


Mesh* MoleculeMesh::Clone() const
{
    return new MoleculeMesh(*this);
}

void LoadFuzzyPoints(const Atom& n, const Cluster& cl, int m, Vector<RVec3>::iterator& p, Vector<double>::iterator& w,int& numpoints)
{
    assert(m>=0);
//
//  First find the index for this Atom.
//
    int index=0,i=1;
    for (auto atom:cl) 
    {
        if (atom->itsR==n.itsR) 
            index=i;
        i++;        
    }
    if (index==0)
    {
        std::cerr << "MoleculeMesh::LoadFuzzy index for Atom not found in cluster" << std::endl;
        exit(-1);
    }

    int na=cl.GetNumAtoms();
    Vector<RVec3>  nuclearPositions=GetPositions(cl);
    Matrix<double> s(na,na);
    Vector<double> P(na);

    MeshBrowser mb(*n.GetIntegrationMesh());
    for (; mb; mb++)
    {
        GetCutoffProfiles(s,nuclearPositions,mb.R(),m);
        CalcCellFunctions(P,s);

        if(P(index)>0)
        {
            double relativeWeight=P(index)/Sum(P);
            *p=mb.R();
            *w=mb.W()*relativeWeight;
            p++;
            w++;
            numpoints++;
        }
    }
}

//
//  Load up a vector with all atom coordinates
//
Vector<RVec3> GetPositions(const Cluster& cl)
{
    Vector<RVec3> ret(cl.GetNumAtoms());
    int i=1;
    for (auto atom:cl) ret(i++)=atom->itsR;
    return ret;
}

void GetCutoffProfiles(Matrix<double>& S, const Vector<RVec3>& R, const RVec3& r, int m)
{
    assert(m>=0);
    int n=R.size();
    Matrix<double>::Subscriptor      Ss(S);
    for (int i=1; i<=n; i++)
    {
        Ss(i,i)=0.0;
        for (int j=1; j<=n; j++)
            if (i!=j) Ss(i,j)=CutoffProfile(R(i),R(j),r,m);
    }
}

double CutoffProfile(const RVec3& Ra,const RVec3& Rb,const RVec3& r,int m)
{
    assert(!(Ra-Rb)!=0);
    assert(m>=0);
    double u=(norm(r-Ra)-norm(r-Rb))/norm(Ra-Rb);
    return Poly(u,m);
}

double Poly(double u,int m)
{
    assert(u<= 1.000000000001);
    assert(u>=-1.000000000001);
    assert(m< 10);
    assert(m>=0);
    for (; m>0; m--) u=0.5*u*(3-u*u);
    return 0.5*(1-u);
}

void CalcCellFunctions(Vector<double>& P, const Matrix<double>& S)
{
    int n=S.GetNumRows();
    Fill(P,1.0);
    Vector<double>::Subscriptor      Ps(P);
    for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
            if (i!=j) Ps(i)*=S(i,j);
}

