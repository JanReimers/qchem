// File: ClusterBrowser.H  Atom cluster browser.
#ifndef _ClusterBrowser_H_
#define _ClusterBrowser_H_



#include "Cluster/Cluster.H"
#include "Misc/ptr_vector.h"

class Atom;

class ClusterBrowser
{
public:
    ClusterBrowser(const Cluster& cl);

    const Atom* operator++(int)
    {
        current++;
        return &current;
    }
    const Atom& operator *() const
    {
        return *current;
    }
    const Atom* operator->() const
    {
        return &current;
    }

    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector<Atom*>::const_iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};

#endif
