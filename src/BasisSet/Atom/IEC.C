// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.
module; 
#include <vector>
#include <set>
export module qchem.BasisSet.Atom.IEClient;
export import qchem.BasisSet.Imp.IEClient;
export import qchem.Irrep_BS;
export import oml.Vector;

export struct AtomIrrepIEClient : public virtual ::IrrepIEClient
{
    AtomIrrepIEClient(size_t N) :  es(N), ns(N) {};
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
    std::vector<size_t> es_indices; //Unique exponent index
    
    size_t n,l;
    std::vector<int> ml;
    
    void Init(const Vector  <double>& exponents,const Vector<double>& norms, size_t l);
    void Init(const Vector  <double>& exponents,const Vector<double>& norms, size_t l,  const std::vector<int>& ml);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,double,double> bf_tuple;
    bf_tuple operator()(size_t i) const {return std::make_tuple(n,l,es(i),ns(i));}
    auto indices() const {return es.indices();}
    auto indices(size_t i) const {return es.indices(i);}

    static       AtomIrrepIEClient* dcast(      ::IrrepIEClient*);
    static const AtomIrrepIEClient* dcast(const ::IrrepIEClient*);
    static const AtomIrrepIEClient& dcast(const ::IrrepIEClient& iec) {return *dcast(&iec);}
    static       AtomIrrepIEClient* dcast(      ::IrrepBasisSet*);
    static const AtomIrrepIEClient* dcast(const ::IrrepBasisSet*);
    static const AtomIrrepIEClient& dcast(const ::IrrepBasisSet& ibs) {return *dcast(&ibs);}
};

