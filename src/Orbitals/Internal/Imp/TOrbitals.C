// File: Imp/TOrbitals.C  
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <algorithm>
module qchem.Orbitals.Internal.OrbitalsImp;
import qchem.Orbitals.Internal.OrbitalImp;
import qchem.ChargeDensity.Factory;
import qchem.Symmetry;
import qchem.Math;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Blaze;

namespace qchem::Orbitals
{

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalsImp<T>::
TOrbitalsImp(const tobs_t<T>* bs, Spin ms)
    : itsBasisSet(bs)
    , itsQNs(bs->GetIrrep(ms))
    , itsD(blazem::zeroH<T>( bs->GetNumFunctions()))
    , itsM(bs->GetNumFunctions())   // default m=n (no truncation); UpdateOrbitals corrects it from U'
{
    assert(itsBasisSet->GetNumFunctions()>0);
};

template <class T> TOrbitalsImp<T>::~TOrbitalsImp()
{
 
}

//-----------------------------------------------------------------
//
//  Orbitals stuff.
//
template <class T> size_t  TOrbitalsImp<T>::GetNumOrbitals() const
{
    return itsOrbitals.size();
}
template <class T> size_t  TOrbitalsImp<T>::GetNumOccOrbitals() const
{
    size_t  n=0;
    for (auto o:this->Iterate())
        if (o->IsOccupied())
            n++;

    return n;
}
template <class T> double TOrbitalsImp<T>::GetEigenValueChange(const Orbitals& og) const
{
    // No UT coverage
    // TODO: OrbitalGroup should return a vector of energies.
    double del=0;
    auto b2=og.Iterate().begin();
    for (auto b1:this->Iterate())
    {
        const Orbital* o2=*b2;
        del+=Square(b1->GetEigenEnergy()-o2->GetEigenEnergy());
        ++b2;
    }
    return sqrt(del);
}

//
//  This is where the real SCF work gets done.
//
template <class T> void TOrbitalsImp<T>::UpdateOrbitals(const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e)
{
    itsOrbitals.clear();
    itsM = UPrime.rows();   // orthonormal-basis dim (n-k after truncation); sizes D' in BuildDensity
    size_t  n=e.size();
    //
    //  Strip out all the positron orbitals.
    //
    static const double e_positron=-c_light*c_light; //Max positron energy = -mc^2.
    size_t index=1;
    for (size_t  i=0; i<n; i++)
    {
 //               std::cout << "o=" << o->GetEigenEnergy() << std::endl;
        if (e[i]<=e_positron) continue; //Strip out all the positron orbitals.
        size_t principle_QN=itsQNs.sym->GetPrincipleOffset() + index++;
        Orbital* o=new TOrbitalImp<T>(itsBasisSet,blazem::column(U,i), blazem::column(UPrime,i), e[i],Orbital_QNs(principle_QN,itsQNs));
        itsOrbitals.push_back(std::unique_ptr<Orbital>(o));

    }
}
template <class T> typename TOrbitalsImp<T>::ds_t TOrbitalsImp<T>::TakeElectrons(double ne)
{
    // Dump electrons into orbitals, starting from the lowest energy.
    for (auto o:this->Iterate())
    {
        ne=o->TakeElectrons(ne);
        if (ne<=0.0) break;
    }
    return BuildDensity(ne);
}

// MOM: occupy the highest-\a priority (max overlap onto the reference) orbitals first, rather than the
// lowest energy -- the within-irrep occupied-subspace continuity that stops a diving diffuse virtual from
// being aufbau-captured (doc/GPWPlan §0b″).  Stable order: ties keep the stored (energy) order.
template <class T> typename TOrbitalsImp<T>::ds_t TOrbitalsImp<T>::TakeElectrons(double ne, const rvec_t& priority)
{
    assert(priority.size()==itsOrbitals.size());
    std::vector<size_t> order(itsOrbitals.size());
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return priority[a]>priority[b];});
    for (size_t idx:order)
    {
        ne=itsOrbitals[idx]->TakeElectrons(ne);
        if (ne<=0.0) break;
    }
    return BuildDensity(ne);
}

// Fermi-Dirac fill (doc/GPWPlan1.md 4b): solve the chemical potential μ by bisection so Σ_i g_i f_i = ne,
// set occ_i = g_i f_i, build D/D', and return the Mermin free-energy term −TS = kT Σ_i g_i[f ln f + (1−f)
// ln(1−f)] ≤ 0 in the leftover slot of the ds_t tuple.  The entropy NEVER touches H: f enters only D (the
// density build below is the SAME AddDensityMatrix path aufbau uses, which already scales |C⟩⟨C| by occ).
// N(μ)=Σ g_i f_i is monotone increasing, so μ is found by plain bisection on an energy window padded by
// many kT (f saturates to 0/1 there).  The cure for near-gapless occupation flapping (NaF Ecut=160).
template <class T> typename TOrbitalsImp<T>::ds_t TOrbitalsImp<T>::TakeElectronsFermi(double ne, double kT)
{
    assert(kT>0.0);
    assert(!itsOrbitals.empty());
    // Fermi occupancy fraction f∈[0,1], overflow-guarded (x=(ε−μ)/kT beyond ±40 saturates).
    auto ffrac=[kT](double e, double mu)->double
    {
        double x=(e-mu)/kT;
        if (x> 40.0) return 0.0;
        if (x<-40.0) return 1.0;
        return 1.0/(1.0+exp(x));
    };
    // Total capacity Σ g_i must cover ne (else the basis is too small -- same contract as aufbau's assert),
    // and the orbital energy range brackets μ.
    double cap=0.0;
    double emin=itsOrbitals.front()->GetEigenEnergy(), emax=emin;
    for (auto o:this->Iterate())
    {
        cap += o->GetDegeneracy();
        emin = std::min(emin,o->GetEigenEnergy());
        emax = std::max(emax,o->GetEigenEnergy());
    }
    assert(ne<=cap+1e-9 && "Fermi fill: too few orbitals for the electron count (basis too small)");
    auto count=[&](double mu)->double
    {
        double N=0.0;
        for (auto o:this->Iterate()) N += o->GetDegeneracy()*ffrac(o->GetEigenEnergy(),mu);
        return N;
    };
    // Bisect μ.  N(lo)≈0≤ne and N(hi)≈cap≥ne by the padding + the capacity assert above.
    double lo=emin-50.0*kT, hi=emax+50.0*kT;
    for (int it=0; it<100 && (hi-lo)>1e-14*(1.0+fabs(hi)); ++it)
    {
        double mu=0.5*(lo+hi);
        if (count(mu)<ne) lo=mu; else hi=mu;
    }
    const double mu=0.5*(lo+hi);
    // Set occupations occ_i=g_i f_i and accumulate −TS = kT Σ g[f ln f + (1−f) ln(1−f)] (x ln x → 0 at 0,1).
    double minusTS=0.0;
    for (auto o:this->Iterate())
    {
        double g=o->GetDegeneracy();
        double f=ffrac(o->GetEigenEnergy(),mu);
        o->SetOccupation(g*f);
        if (f>1e-15 && f<1.0-1e-15) minusTS += kT*g*(f*log(f)+(1.0-f)*log(1.0-f));
    }
    hmat_t<T> DPrime=std::get<1>(BuildDensity(0.0));   // build D/D' from the fractional occupations
    return std::make_tuple(minusTS,std::move(DPrime)); // {−TS, D'} (μ solved => no leftover electrons)
}

// Build D (and the orthonormal-basis D') from the currently-occupied orbitals; shared by both
// TakeElectrons paths (energy-order and MOM-order).  \a ne = leftover electrons, forwarded to the caller.
template <class T> typename TOrbitalsImp<T>::ds_t TOrbitalsImp<T>::BuildDensity(double ne)
{
    itsD=blazem::zeroH<T>(itsD.rows());                  // AO density: full n x n
    hmat_t<T> DPrime(blazem::zeroH<T>(itsM));            // orthonormal density: m x m (m=n-k after truncation)
    for (auto o:Iterate<TOrbital<T>>()) o->AddDensityMatrix(itsD,DPrime);   // was hardcoded TOrbital<double>
    return std::make_tuple(ne,DPrime);
}

template <class T> tDM_CD<T>* TOrbitalsImp<T>::GetChargeDensity() const
{
    // Scale by the irrep's BZ weight (1 for atoms/molecules; w_k for a Bloch k-point) so the composite
    // density sums to the BZ average Sum_k w_k rho_k.  The orbitals' occupations stay physical (the
    // accelerator works on the unscaled D'); only the charge-density object carries the weight.
    double w=GetQNs().sym->GetWeight();
    hmat_t<T> Dw = (w==1.0) ? itsD : hmat_t<T>(w*itsD);
    return qchem::ChargeDensity::IrrepCD_Factory(Dw,itsBasisSet,GetQNs());
}


template <class T>  Irrep TOrbitalsImp<T>::GetQNs() const
{
    return itsQNs;
}


//-----------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> vec_t<T> TOrbitalsImp<T>::operator()(const rvec3_t& r) const
{
    vec_t<T> ret(GetNumOrbitals());
    auto i(ret.begin());
    // No UT coverage
    for (auto b:Iterate<TOrbital<T>>()) 
    {
        *i=(*b)(r);
        i++;
    }
    return ret;
}

template <class T> vec3vec_t<T> TOrbitalsImp<T>::Gradient(const rvec3_t& r) const
{
    // No UT coverage
    vec3vec_t<T> ret(GetNumOrbitals());
    auto i(ret.begin());
    for (auto b:Iterate<TOrbital<T>>()) 
    {
        *i=b->Gradient(r);
        i++;
    } 
    return ret;
}

//-----------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& TOrbitalsImp<T>::Write(std::ostream& os) const
{
    os << "        Orbital group with " << GetNumOrbitals() << " " << itsBasisSet->GetSymmetry() << "orbitals:" << std::endl;
    os << "            Occupation      Energy      Eigenvector" << std::endl;
    os << itsOrbitals;
    return os;
}

template class TOrbitalsImp<double>;
template class TOrbitalsImp<dcmplx>;

} //namespace