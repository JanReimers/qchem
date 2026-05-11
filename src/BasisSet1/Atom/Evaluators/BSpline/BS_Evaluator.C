// File: BasisSet/Atom/BSpline/NR/BS_Evaluator.C
module;
export module BasisSet.Atom.BSpline.NR.BS_Evaluator;
export import qchem.BasisSet1.Atom.BS_Evaluator;
export import BasisSet.Atom.BSpline.NR.IBS_Evaluator_r;
import qchem.BasisSet1.Atom.BSpline.SplineGrouper;
import qchem.BasisSet1.Atom.BSpline.Rk;
import qchem.Basisset.Atom.BSpline.GLQuadrature;
import BasisSet.Atom.BSpline.NR.IBS_Evaluator;

export template <size_t K> class BSpline_BS 
    : public virtual BS_Evaluator
{
public:
    using IBS_Evaluator_t = BSpline_IBS<K>;
    BSpline_BS() : itsRkCache(0) {};
    virtual ~BSpline_BS();
    virtual void Register(IBS_Evaluator *);
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual double loop_4_direct  (size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 
    virtual double loop_4_exchange(size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 
protected:
    void BuildCache(size_t lmax);
private:
    const GLCache* GetGL(size_t l) const;
    SplineGrouper<K> grouper;
    BSpline::RkCache<K>* itsRkCache;
};

export template <size_t K> class BSpline_r_BS 
    : public virtual BS_Evaluator
{
public:
    using IBS_Evaluator_t = BSpline_r_IBS<K>;
    BSpline_r_BS() : itsRkCache(0) {};
    virtual ~BSpline_r_BS();
    virtual void Register(IBS_Evaluator *);
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual double loop_4_direct  (size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 
    virtual double loop_4_exchange(size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 
protected:
    void BuildCache(size_t lmax);
private:
    const GLCache* GetGL(size_t l) const;
    SplineGrouper<K> grouper;
    BSpline::RkCache_r<K>* itsRkCache;
};

