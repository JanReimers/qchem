// File: BasisSet1/Atom/Evaluators/BSpline/BS_Evaluator.C
module;
export module qchem.BasisSet1.Atom.Evaluators.BSpline.BS;
export import qchem.BasisSet1.Atom.Evaluators.BS;
export import qchem.BasisSet1.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet1.Atom.Evaluators.BSpline.IBS;

export template <size_t K> class BSpline_BS_Evaluator 
    : public virtual BS_Evaluator
{
public:
    using IBS_Evaluator_t = BSpline_IBS_Evaluator<K>;
    BSpline_BS_Evaluator() : itsRkCache(0) {};
    virtual ~BSpline_BS_Evaluator();
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

export template <size_t K> class BSpline_r_BS_Evaluator 
    : public virtual BS_Evaluator
{
public:
    using IBS_Evaluator_t = BSpline_r_IBS_Evaluator<K>;
    BSpline_r_BS_Evaluator() : itsRkCache(0) {};
    virtual ~BSpline_r_BS_Evaluator();
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

