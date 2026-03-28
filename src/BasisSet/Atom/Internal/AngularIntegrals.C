// File: AngularIntegrals.C  Handle the angular part of 2-electron ERIs
export module qchem.BasisSet.Atom.Internal.AngularIntegrals;
export import qchem.Types;
import oml.Vector;

export template <class T> const Vector<T>& operator+=(Vector<T>& a, const Vector<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    return ArrayAdd(a,b);
}

export namespace AngularIntegrals
{
    typedef Vector<double> RVec;

    double Exchange(int k,int la,int lb);
    double Coulomb (int k,int la,int lc,int ma,int mc);            
    double Exchange(int k,int la,int lb,int ma,int mb);            

    RVec Coulomb (int la,int lc);            
    RVec Exchange(int la,int lb);            
    RVec Coulomb (int la,int lc,int ma,int mc);            
    RVec Exchange(int la,int lb,int ma,int mb);            

};
