// File: ClusterBrowser.C  Orbital group browser.


#include "Cluster/Cluster.H"
#include "Cluster/ClusterBrowser.H"
#include "Molecule.H"

//#######################################################################
//
//  Browser constructor..
//

ClusterBrowser::ClusterBrowser(const Cluster& cl)
{
    const Molecule* m=dynamic_cast<const Molecule*>(&cl);
    assert(m);
    begin=m->itsAtoms.begin();
    current=begin;
    end=m->itsAtoms.end();
}

