// File: MoleculeMesh.C  mesh implementation



#include "Imp/Cluster/MoleculeMesh.H"
#include "Imp/Cluster/Atom.H"
#include <Cluster.H>
#include <MeshParams.H>
#include "oml/matrix.h"
#include "oml/vector.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

Vector<RVec3>  GetPositions     (const Cluster&);
void            GetCutoffProfiles(Matrix<double>&,const Vector<RVec3>&, const RVec3&, int m_mu);
double          CutoffProfile    (const RVec3&,const RVec3&,const RVec3&,int m_mu);
double          Poly             (double,int m_mu);
void            CalcCellFunctions(Vector<double>&, const Matrix<double>&);

//
//  Use Becke's fuzzy polyedra algorithm for integrating over a molecule.
//  See A. D. Becke, J. Chem. Phys, 88(4), page 2547 (1988).
//
MoleculeMesh::MoleculeMesh(const Cluster& cl, const MeshParams& mp)
{
    for (auto atom:cl) LoadFuzzyPoints(*atom,cl,mp);
}


Mesh* MoleculeMesh::Clone() const
{
    return new MoleculeMesh(*this);
}

void MoleculeMesh::LoadFuzzyPoints(const Atom& n, const Cluster& cl,const MeshParams& mp)
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

    Mesh* am= n.CreateMesh(mp);
    for (auto rw: *am)
    {
        GetCutoffProfiles(s,nuclearPositions,::r(rw),mp.m_mu);
        CalcCellFunctions(P,s);

        if(P(index)>0)
        {
            double relativeWeight=P(index)/Sum(P);
            push_back(::r(rw),::w(rw)*relativeWeight);
        }
    }
    delete am;
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

void GetCutoffProfiles(Matrix<double>& S, const Vector<RVec3>& R, const RVec3& r, int m_mu)
{
    assert(m>=0);
    int n=R.size();
    Matrix<double>::Subscriptor      Ss(S);
    for (int i=1; i<=n; i++)
    {
        Ss(i,i)=0.0;
        for (int j=1; j<=n; j++)
            if (i!=j) Ss(i,j)=CutoffProfile(R(i),R(j),r,m_mu);
    }
}

double CutoffProfile(const RVec3& Ra,const RVec3& Rb,const RVec3& r,int m_mu)
{
    assert(norm(Ra-Rb)!=0);
    assert(m_mu>=0);
    double u=(norm(r-Ra)-norm(r-Rb))/norm(Ra-Rb);
    return Poly(u,m_mu);
}

double Poly(double u,int m_mu)
{
    assert(u<= 1.000000000001);
    assert(u>=-1.000000000001);
    assert(m_mu< 10);
    assert(m_mu>=0);
    for (; m_mu>0; m_mu--) u=0.5*u*(3-u*u);
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

