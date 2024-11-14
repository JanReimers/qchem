
#include "Imp/BasisSet/SphericalGaussian_m/IEClient.H"
#include "Imp/Integrals/GaussianIntegrals.H"

using std::cout;
using std::endl;

template <class T> void FillPower(Vector<T>& arr,T start, T stop);
//template void FillPower(Vector<double>& arr,double start, double stop);

namespace SphericalGaussian_m
{
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=GaussianNorm(es(i),Ls(i));
}

void IrrepIEClient::Init(std::set<double>& exponents,size_t L, int m)
{
    int i=1;
    for (auto& e:exponents) es(i++)=e;
    Fill(Ns,L+1);
    Fill(Ls,L);
    Fill(Ms,m);
    for (auto i:es.indices())  ns(i)=GaussianNorm(es(i),Ls(i));
}



void IEClient::Append(const IrrepIEClient* ic)
{
    itsIrreps.push_back(ic);
    
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    Ms.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        Ms(j)=ic->Ms(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
        
        BFGrouper::Append(es(j),Ls(j),j);
    }
//    for (auto e:unique_esv) cout << e << " ";
//    cout << endl;
//    for (auto e:unique_es) cout << e.first << " ";
//    cout << endl;
//    for (auto i:es_indexes) cout << i << " ";
//    cout << endl;

}


using std::cout;
using std::endl;


const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}




} //namespace
