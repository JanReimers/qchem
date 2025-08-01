// File: EnergyLevel.C  Energy level with degeneracy and orbital list.
module;
#include <iosfwd>
#include <map>
export module qchem.EnergyLevel;
export import qchem.Symmetry.Orbital;
export import qchem.Orbitals;

export struct EnergyLevel
{
    EnergyLevel(const Orbital* o);
    EnergyLevel(const EnergyLevel&);

    void merge(const EnergyLevel&);
    void Report(std::ostream&) const;
    
    double         e,occ;
    int            degen;
    Orbital_QNs    qns;
  
};

export class EnergyLevels
{
    typedef std::multimap<double,EnergyLevel> el_t;
    typedef std::map<Orbital_QNs,EnergyLevel> oel_t;
    typedef el_t::const_iterator const_iterator;
public:
    void insert(const EnergyLevel& el)
    {
        itsELevels.insert(std::make_pair(el.e,el));
        itsQNLevels.insert(std::make_pair(el.qns,el));
    }
    void clear() 
    {
        itsELevels.clear();
        itsQNLevels.clear();
    }
    void merge(const EnergyLevels& els);
    void merge(const EnergyLevels& els, double tol);
    const EnergyLevel& find(const Orbital_QNs&) const;

    const_iterator begin() const {return itsELevels.begin();}
    const_iterator end  () const {return itsELevels.end  ();}
    size_t         size () const {return itsELevels.size ();}
    void Report(std::ostream&) const;

private:

     el_t  itsELevels;
     oel_t itsQNLevels; //All levels indexed by QNs.
};

