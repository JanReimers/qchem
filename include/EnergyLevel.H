// File: EnergyLevel.H  Energy level with degeneracy and orbital list.
#ifndef _EnergyLevel_H_
#define _EnergyLevel_H_

#include <types.H>
#include <Spin.H>
#include <Orbital_QNs.H>

struct EnergyLevel
{
    EnergyLevel(double _e, double _occ, double _g,const Orbital_QNs& _qns, Orbital* _o)
        : e(_e), occ(_occ), degen(_g), qns(_qns), orbital(_o) 
        {};
    EnergyLevel(const EnergyLevel&);
    ~EnergyLevel();

    void merge(const EnergyLevel&);
    void Report(std::ostream&) const;
    
    double         e,occ;
    int            degen;
    Orbital_QNs    qns;
    Orbital*       orbital;
  
};

#include <map>

class EnergyLevels
{
    typedef std::multimap<double,EnergyLevel> el_t;
    typedef el_t::iterator iterator;
    typedef el_t::const_iterator const_iterator;
public:
    void insert(const EnergyLevel& el)
    {
        itsELevels.insert(std::make_pair(el.e,el));
    }
    void clear() {itsELevels.clear();}
    void merge(const EnergyLevels& els);
    void merge(const EnergyLevels& els, double tol);
    iterator       begin()       {return itsELevels.begin();}
    iterator       end  ()       {return itsELevels.end  ();}
    const_iterator begin() const {return itsELevels.begin();}
    const_iterator end  () const {return itsELevels.end  ();}
    
    void Report(std::ostream&) const;

private:

     el_t itsELevels;
};


#endif //_EnergyLevel_H_
