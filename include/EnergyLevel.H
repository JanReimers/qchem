// File: EnergyLevel.H  Energy level with degeneracy and orbital list.
#ifndef _EnergyLevel_H_
#define _EnergyLevel_H_

#include <types.H>
#include <Spin.H>

struct EnergyLevel
{
    EnergyLevel(double _e, double _occ, double _g,const QuantumNumber& _qn,const Spin& _s, Orbital* _o)
        : e(_e), occ(_occ), degen(_g), qn(_qn), s(_s), orbital(_o) 
        {};

    double               e,occ;
    int                  degen;
    const QuantumNumber& qn;
    Spin                 s;
    Orbital*             orbital;
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
    void merge(const EnergyLevels& els)
    {
        for (auto& el:els) itsELevels.insert(el);
    }
    iterator       begin()       {return itsELevels.begin();}
    iterator       end  ()       {return itsELevels.end  ();}
    const_iterator begin() const {return itsELevels.begin();}
    const_iterator end  () const {return itsELevels.end  ();}
private:
     el_t itsELevels;
};


#endif //_EnergyLevel_H_
