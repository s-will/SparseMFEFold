/*
  Zuker with about half of the sugar.  Sparse variant of an
  RNA-folding algorithm "in-between Nussinov and Zuker" (interior loop
  energies, but no proper multi-loop energies)
  
  Demonstration of sparsification and sparse trace-back. For the
  purpose of trace-back, this algorithm is fundamentally more complex
  than the Nussinov algorithm.

  Recursions:
  W(i,j) = min { W(i,k-1),
                 min_i<k<j  W(i,k-1) + V(k,j),
                 V(i,j)
	       }
  V(i,j) = min { HairpinE(i,j),
		 min_kl V(i,j)+ILoopE(i,j,k,l),
		 W(i+1,j-1) + NonILoopPenalty if pair_type(S[i],S[j])>0
               }


  TODO: trace arrows to candidates could be omitted and reconstructed in trace back


*/

/*
  TODO: extend to Zuker-Recursions:
  W(i,j) = min { W(i,j-1),
                 min_i<k<j  W(i,k-1) + V(k,j), <-- CLW
                 V(i,j)
		 0 if i>=j-m
	       }
  V(i,j) = min { HairpinE(i,j),
		 min_kl V(i,j)+ILoopE(i,j,k,l),
		 min_k WM(i+1,k-1) + WM(k,j) + a  <-- CLV
		}
  WM(i,j) = min { WM(i,j-1) + c,
                  WM(i+1,j) + c,
                  V(i,j)+b,
		  min_i<k<j  WM(i,k-1) + V(k,j) + b,  <-- CLWM
		  -inf if i>=j-m
		  }
*/

#include <LocARNA/matrices.hh>

#include <limits>

#include <vector>
//#include <map>
//#include <unordered_map>

#include <cstring>

extern "C" {
#   include "ViennaRNA/pair_mat.h"
#   include "ViennaRNA/loop_energies.h"
}

class HalfZuker {

public:
    //! type of energy
    typedef int energy_t;

private:
    const energy_t NonILoopPenalty_;
    const energy_t ILoopBonus_;

    std::string seq_;
    size_t n_;

    short *S_;
    short *S1_;
    paramT *params_;
    
    LocARNA::Matrix<energy_t> V_;
    std::vector<energy_t> W_;

public:
    typedef std::pair<size_t,energy_t> cand_entry_t;
    typedef std::vector< cand_entry_t > cand_list_t;
private:
    std::vector< cand_list_t > CLW_; //!< candidate list for decomposition in W
    
    std::string structure_;
    
    size_t ta_count_;
    size_t ta_erase_;
    size_t ta_max_;

public:
    /**
     * @brief Trace arrow
     * 
     * Describes a trace arrow from i,j to k_,l_ of give type
     * type_. The sourcce (i,j) is not represented in the data
     * structure. However, each trace arrow is associated to exactly
     * one source.
     */
    class TraceArrow {
	char type_; //!< type of the arrow (pointing to 'C' or 'G')
	unsigned char k_; //!< target row of arrow
	unsigned char l_; //!< target column of arrow
	uint count_; //!< counts how many trace arrows point to the source
    public:
	/**
	 * @brief construct by target coordinates
	 * @param type
	 * @param k
	 * @param l
	 */
	TraceArrow(char type,
		   size_t i,
		   size_t j,
		   size_t k,
		   size_t l
		   )
	    : type_(type),
	      k_(k-i),
	      l_(j-l),
	      count_(0)
	{}

	/**
	 * @brief empty c'tor
	 */ 
	TraceArrow() {}
	
	bool is_W() const {return type_=='G';}
	bool is_V() const {return type_=='C';}
	size_t k(size_t i,size_t j) const {return k_+i;}
	size_t l(size_t i,size_t j) const {return j-l_;}
	size_t source_ref_count() const {return count_;}
	
	void inc_src() {count_++;}
	void dec_src() {count_--;}

    };
private:
    
    /* space saving replacement for map of trace arrows in rows i;
       works for our special case */
    template<class key_t, class val_t>
    class SimpleMap: std::vector<std::pair<key_t, val_t> > {
	typedef std::pair<key_t, val_t> key_val_t;
	typedef std::vector<key_val_t> key_val_vec_t;
	typedef typename key_val_vec_t::iterator iterator;
	typedef typename key_val_vec_t::const_iterator const_iterator;
	
	class  {
	public:
	    bool
	    operator () (const key_val_t &x,
			 const key_t &y) const {
		return x.first < y;
	    }
	} comp;

	iterator
	binsearch (iterator first, iterator last, const key_t& key)
	{
	    first = std::lower_bound(first,last,key,comp);
	    if (first==last || key < first->first) {
		return last;
	    }
	    return first;
	}

	const_iterator
	binsearch (const_iterator first, const_iterator last, const key_t& key) const
	{
	    first = std::lower_bound(first,last,key,comp);
	    if (first==last || key < first->first) {
		return last;
	    }
	    return first;
	}

    public:
	SimpleMap() {}
	
	const_iterator 
	find(const key_t &key) const { 
	    auto it= binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	    assert(it == key_val_vec_t::end() || it->first == key);
	    return it;
	};
	
	iterator
	find(const key_t &key) { 
	    auto it=binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	    assert(it == key_val_vec_t::end() || it->first == key);
	    return it;
	};
    	
	bool
	exists(const key_t &key) const {
	    return find(key) != key_val_vec_t::end();
	}
	
	/** 
	 * @brief push in ascending order of keys 
	 * @param key
	 * @param val
	 *
	 * successive push_ascending must be in ascending order of the key type
	 */
	void
	push_ascending( const key_t &key, const val_t &val ) {
	    assert(size()==0||key > key_val_vec_t::operator[](size()-1).first);
	    key_val_vec_t::push_back(key_val_t(key,val));
	}
	
	void
	erase(iterator it) {
	    //std::cout << key_val_vec_t::size()<<" "<<key_val_vec_t::capacity()<<std::endl;
	    
	    // doing only this, wastes space (due to stl-vector allocation strategy)
	    // key_val_vec_t::erase(it);
	    
	    // instead: copy to new space with exactly the right size
	    // (possible optimization: only mark erased entries at
	    // first and delay the copying until it pays off)
	    key_val_vec_t vec(key_val_vec_t::size()-1);
	    iterator target=copy(key_val_vec_t::begin(),it,vec.begin());
	    copy(it+1,key_val_vec_t::end(),target);
	    
	    vec.swap(*this);
	}

	size_t
	size() const {
	    return key_val_vec_t::size();
	}
	
	size_t
	capacity() const {
	    return key_val_vec_t::capacity();
	}

	void
	reallocate() {
	    key_val_vec_t vec(size());
	    copy(key_val_vec_t::begin(),key_val_vec_t::end(),vec.begin());
	    vec.swap(*this);
	}
    };
    
    
    //typedef std::map< size_t, TraceArrow >        trace_arrow_row_map_t;
    typedef SimpleMap< size_t, TraceArrow >       trace_arrow_row_map_t;
    typedef std::vector< trace_arrow_row_map_t >  trace_arrow_map_t;

    trace_arrow_map_t trace_arrow_;
    
    /**
     * Get target of trace arrow by source (const)
     *
     * @param i source row index
     * @param j source column index
     */
    const TraceArrow &
    trace_arrow_from(size_t i, size_t j) const {
	return trace_arrow_[i].find(j)->second;
    }

    /**
     * Get target of trace arrow by source (non-const)
     *
     * @param i source row index
     * @param j source column index
     */
    TraceArrow &
    trace_arrow_from(size_t i, size_t j) {
	return trace_arrow_[i].find(j)->second;
    }

    /**
     * Check existence of trace arrow by source
     *
     * @param i source row index
     * @param j source column index
     * @returns whether trace arrow exists
     */
    bool
    exists_trace_arrow_from(size_t i, size_t j) const {
	return trace_arrow_[i].exists(j);
    }

    /** 
     * Register trace arrow
     * 
     * @param i source row
     * @param j source column
     * @param type target type
     * @param k target row
     * @param l target column
     */    
    void
    register_trace_arrow(size_t i, size_t j,char type, size_t k, size_t l, energy_t e) {
	// std::cout << "register_trace_arrow "<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
	trace_arrow_[i].push_ascending( j, TraceArrow(type,i,j,k,l) );
	
	inc_source_ref_count(k,l);
	
	ta_count_++;
	ta_max_ = std::max(ta_max_,ta_count_);
    }

    /**
     * Increment the reference count of the source
     *
     * @param i source row index
     * @param j source column index
     *
     * If no trace arrow from source exists, do nothing
     */
    void
    inc_source_ref_count(size_t i, size_t j) {
	// get trace arrow from (i,j) if it exists
	if (! trace_arrow_[i].exists(j)) return;
	
	auto it=trace_arrow_[i].find(j);
	
	TraceArrow &ta=it->second;
	
	ta.inc_src();
    }

    /**
     * Garbage collect trace arrow
     *
     * if count = 0 then remove trace arrow; recursively decrement
     * targets and remove if count drops to 0
     */
    void
    gc_trace_arrow(size_t i, size_t j) {
	
	assert( trace_arrow_[i].exists(j) );

	auto col = trace_arrow_[i].find(j);
	
	const auto &ta = col->second;
	
	if (ta.source_ref_count() == 0) {
	    // get trace arrow from the target if the arrow exists
	    if (exists_trace_arrow_from(ta.k(i,j),ta.l(i,j))) {
		auto &target_ta = trace_arrow_from(ta.k(i,j),ta.l(i,j));
		
		target_ta.dec_src();
		
		gc_trace_arrow(ta.k(i,j),ta.l(i,j));
	    }
	    
	    trace_arrow_[i].erase(col);
	    ta_count_--;
	    ta_erase_++;
	}
    }
    
    /** 
     * Garbage collect all trace arrows in given row
     * 
     * @param i index of the row 
     */
    void
    gc_row_of_trace_arrows( size_t i ) {
	
	assert(i<=n_);
	
	for (size_t j=1; j<=n_ ; j++) {
	    if (! trace_arrow_[i].exists(j)) continue;
	    gc_trace_arrow(i,j);
	}
	
    }

    struct {
	bool operator ()(const cand_entry_t &x, size_t y) const {
	    return x.first > y;
	}
    }
    cand_comp;
    
    /** 
     * Test existence of W-candidate
     * 
     * @param i row
     * @param j col
     * 
     * @return whether (i,j) is candidate for W
     */
    bool
    is_candidateW(size_t i, size_t j) const {
	const cand_list_t &list = CLW_[j];
	auto it = std::lower_bound(list.begin(),list.end(),j,cand_comp); 
	bool res = it != list.end() && it->first==i;
	
	return res;
    }
    
public:
    HalfZuker(const std::string &seq)
	: NonILoopPenalty_(50),
	  ILoopBonus_(-300),
	  seq_(seq),
	  n_(seq.length()),
	  S_(encode_sequence(seq.c_str(),0)),
	  S1_(encode_sequence(seq.c_str(),1)),
	  params_(scale_parameters()),
	  ta_count_(0),
	  ta_erase_(0),
	  ta_max_(0)
    {
	make_pair_matrix();
	
	V_.resize(MAXLOOP,n_+1);
	W_.resize(n_+1,0);
	
	/* init candidate lists */
	CLW_.resize(n_+1);
	for (size_t j=n_; j>0; --j) {
	    CLW_[j].push_back( cand_entry_t(j,0) );
	}

	trace_arrow_.resize(n_+1);
    }

    void
    trace_back() {
	structure_.resize(n_+1,'.');

	/* Traceback */
	trace_W(1,n_);
	structure_ = structure_.substr(1,n_);
    
	std::cout << seq_ << std::endl;
	std::cout << structure_ << std::endl;
    }

    ~HalfZuker() {
	free(params_);
	free(S_);
	free(S1_);
    }

private:
    int
    pair_type(size_t i, size_t j) const {
	return pair[S_[i]][S_[j]];
    }
   
    /* pre: ptype_closing>0 */
    energy_t
    ILoopE(int ptype_closing,size_t i, size_t j, size_t k,  size_t l) const {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here
    
	// note: enclosed bp type 'turned around' for lib call
	int ptype_enclosed = rtype[pair_type(k,l)]; 
    
	if (ptype_enclosed==0) return INF;
    
	return
	    E_IntLoop(k-i-1,j-l-1,
		      ptype_closing,
		      ptype_enclosed,
		      S_[i+1],
		      S_[j-1],
		      S_[k-1],
		      S_[l+1],
		      const_cast<paramT *>(params_));
    }

    energy_t
    HairpinE(size_t i, size_t j) const {
    
	assert(1<=i);
	assert(i<j);
	//assert(j<=len); // don't know len here
    
	int ptype_closing = pair_type(i,j);
    
	if (ptype_closing==0) return INF;
    
	return
	    E_Hairpin(j-i-1, 
		      ptype_closing,
		      S_[i+1],
		      S_[j-1],
		      &seq_.c_str()[i-1],
		      const_cast<paramT *>(params_));

    }


    /**
     * Compute row of G
     *
     * @param i row index
     * @param j column index
     */
    void
    compute_W(size_t i, size_t max_j) {
	//std::cout << "Compute G " <<i<<" "<<max_j<<std::endl; 
	
	for ( size_t j=i-1; j<i+TURN+1; j++ ) { W_[j]=0; }
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
	
	    energy_t d = INF;
	    
	    for ( auto it = CLW_[j].begin(); 
		  CLW_[j].end()!=it && it->first>=i ; ++it ) {
		d = std::min( d, W_[it->first-1] + it->second );
	    }
	
	    W_[j] = d;
	}
	
	// std::cout << "G["<<i<<"]["<<i<<".."<<max_j<<"]: " ;
	// for ( size_t j=i; j<=max_j; j++ ) { 
	//     std::cout << W_[j] << " ";
	// }
	// std::cout << std::endl;
    }

    /** 
     * Trace from W entry
     * 
     * @param i row index
     * @param j column index
     * pre: structure is string of size (n+1)
     * pre: W contains values of row i in interval i..j
     */
    void
    trace_W(size_t i, size_t j) {
	// std::cout << "Trace W "<<i<<" "<<j<<std::endl;
	if (i+TURN+1>=j) return;
	
	size_t k=j+1;

	for ( auto it = CLW_[j].begin(); 
	      CLW_[j].end()!=it && it->first>=i;
	      ++it ) {
	    energy_t d_it = W_[it->first-1] + it->second;
	    
	    if (W_[j] == d_it) {
		k = it->first;
		break;
	    }
	}
	
	//std::cout << i<<" "<<j<<" "<<k<<std::endl;
	assert(i<=k && k<=j);
	
	// don't recompute G here, since i is not changed
	trace_W(i,k-1);
    
	if (k!=j) {
	    trace_V(k,j);
	}
    }

    /** 
     * Trace from V entry
     * 
     * @param i row index
     * @param j column index
     * @param energy energy in V[i,j]
     */
    void
    trace_V(size_t i, size_t j) {
	assert (i+TURN+1<=j);
	
	structure_[i]='(';
	structure_[j]=')';

	if (exists_trace_arrow_from(i,j)) {
	    const TraceArrow &arrow = trace_arrow_from(i,j);
	    if(arrow.is_V()) {
		trace_V(arrow.k(i,j),arrow.l(i,j));
	    } else if (arrow.is_W()) {
		compute_W(arrow.k(i,j),arrow.l(i,j));
		trace_W(arrow.k(i,j),arrow.l(i,j));
	    }
	}
    }
public:

    /* recursion evaluation (forward, sparse) */
    void
    forward_evaluation() {
    	
	for (size_t i=n_; i>0; --i) {
	    for ( size_t j=i+TURN+1; j<=n_; j++ ) {
	    
		energy_t d = INF;
		
		for ( auto &x : CLW_[j] ) {
		    d = std::min( d, W_[x.first-1] + x.second );
		}
		
		int ptype_closing = pair_type(i,j);
		
		energy_t c = INF;
		size_t i_mod=i%MAXLOOP;

		if(ptype_closing>0) {

		    c = HairpinE(i,j);
		    V_(i_mod,j) = c;

		    size_t best_l=0;
		    size_t best_k=0;
		    energy_t best_e;
		    
		    // constraints for interior loops
		    // i<k; l<j            
		    // k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
		    //            ==> l >= k+j-i-MAXLOOP-2
		    // l-k>=TURN+1         ==> k <= j-TURN-2
		    //            ==> l >= k+TURN+1
		    // j-i>=TURN+3
		    //
		    size_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
		    for ( size_t k=i+1; k<=max_k; k++) {
			size_t k_mod=k%MAXLOOP;
			
			size_t min_l=std::max(k+TURN+1  +MAXLOOP+2, k+j-i)-MAXLOOP-2;
			
			for (size_t l=min_l; l<j; l++) {
			    
			    assert(k-i+j-l-2<=MAXLOOP);
			    
			    energy_t c_kl=V_(k_mod,l) + ILoopE(ptype_closing,i,j,k,l) + ILoopBonus_;
			    if ( c_kl < c ) {
				c = c_kl;
				best_l=l;
				best_k=k;
				best_e=V_(k_mod,l);
			    }
			}
		    }
	    
		    /*
		      std::cout << i << " " << j << " " << CL[j].size()
		      << " " << best_l << " " << best_k
		      << " " << c << " " << d << std::endl;
		    */
		
		    energy_t c_non_iloop = W_[j-1] + NonILoopPenalty_; // W_(i+1,j-1) is still in W_[j-1] !
		    //register trace arrows from C-entry (i,j); even for non-candidates 
		    if ( c_non_iloop < c ) {
			c = c_non_iloop;
			register_trace_arrow(i,j,'G',i+1,j-1,W_[j-1]);
		    } else {
			if (best_l!=0) { /* not hairpin */
			    assert(best_k<best_l);
			    if (!is_candidateW(best_k,best_l)) {
				//std::cout << "Reg TA "<<best_k<<" "<<best_l<<std::endl;
				register_trace_arrow(i,j,'C',best_k,best_l,best_e);
			    } else {
				//std::cout << "Avoid TA "<<best_k<<" "<<best_l<<std::endl;
			    }
			}
		    }
		    
		    if ( c < d ) {
			// got a candidate => register candidate
			CLW_[j].push_back( cand_entry_t(i, c) );
			//std::cout << "Register candidate "<<i<<" "<<j<<std::endl;
			inc_source_ref_count(i,j); // <- always keep arrows starting from candidates 
		    }
		}
		
		V_(i_mod,j) = c;
		W_[j] = std::min(c,d);
	    }
	    
	    // Clean up trace arrows in i+MAXLOOP
	    if ( i+MAXLOOP <= n_) {
		gc_row_of_trace_arrows( i + MAXLOOP );
	    }
	    
	    // Reallocate candidate lists in i
	    for ( auto &x: CLW_ ) {
		if (x.capacity() > 1.5*x.size()) {
		    cand_list_t vec(x.size());
		    copy(x.begin(),x.end(),vec.begin());
		    vec.swap(x);
		}
	    }
	    for ( auto &x: trace_arrow_ ) {
		if (x.capacity() > 1.2 * x.size()) {
		    x.reallocate();
		}
	    }
	    
	}
	
	std::cout << "MFE: \t"<<(W_[n_]/100.0)<<std::endl;
    }

    size_t
    num_of_candidates() const {
	size_t c=0;
	for ( auto &x: CLW_ ) {
	    c += x.size();
	}
	return c;
    }

    size_t
    capacity_of_candidates() const {
	size_t c=0;
	for ( auto &x: CLW_ ) {
	    c += x.capacity();
	}
	return c;
    }

    size_t
    num_of_tas() const {
	size_t c=0;
	for ( auto &x: trace_arrow_ ) {
	    c += x.size();
	}
	return c;
    }

    size_t
    capacity_of_tas() const {
	size_t c=0;
	for ( auto &x: trace_arrow_ ) {
	    c += x.capacity();
	}
	return c;
    }

    size_t ta_count() const {return ta_count_;}
    size_t ta_erase() const {return ta_erase_;}
    size_t ta_max() const {return ta_max_;}
    
};


int
main(int argc,char **argv) {
        
    std::string seq;
    if (argc>1) {
	seq=argv[1];
    } else {
	std::getline(std::cin,seq);
    }
    HalfZuker hz(seq);

    std::cout << seq << std::endl;
    std::cout << "Len:\t"<<seq.length()<<std::endl<<std::endl;
    
    hz.forward_evaluation();

    hz.trace_back();

    size_t n=seq.length();
    
    float factor=1024;
    const std::string unit=" kB";

    std::cout << "[Unsp:\t"<< 2* n*(n+1)/2 * sizeof(int)/factor<<unit<<"]"<<std::endl;
    std::cout << "W:\t"<<n*sizeof(HalfZuker::energy_t)/factor<<unit<<std::endl;
    std::cout << "V:\t"<<n*sizeof(HalfZuker::energy_t)*MAXLOOP/factor<<unit<<std::endl;
    std::cout << "Cands:\t"<<hz.num_of_candidates()*sizeof(HalfZuker::cand_entry_t)/factor<<unit<<std::endl;
    //std::cout << "TAs:\t"<<hz.ta_count()*sizeof(HalfZuker::TraceArrow)/factor<<unit<<"; size of ta="<<sizeof(HalfZuker::TraceArrow)<<std::endl;

    std::cout << "TAs+ov:\t"<<hz.ta_count()*(sizeof(HalfZuker::TraceArrow)+8)/factor<<unit<<"; size of ta="<<sizeof(HalfZuker::TraceArrow)<<std::endl;

    // trace arrows with overhead of maps
    // std::cout << "TAs+map ov:\t"<<(n*48+hz.ta_count()*(32+sizeof(size_t)+sizeof(HalfZuker::TraceArrow)))/factor<<unit<<"; size of ta="<<sizeof(HalfZuker::TraceArrow)<<std::endl;

    std::cout << "TA cnt:\t"<<hz.ta_count()<<std::endl;
    std::cout << "TA max:\t"<<hz.ta_max()<<std::endl;
    std::cout << "TA rm:\t"<<hz.ta_erase()<<std::endl;
    
    std::cout <<std::endl;
    std::cout << "Can num:\t"<<hz.num_of_candidates()<<std::endl;
    std::cout << "Can cap:\t"<<hz.capacity_of_candidates()<<std::endl;
    std::cout << "TAs num:\t"<<hz.num_of_tas()<<std::endl;
    std::cout << "TAs cap:\t"<<hz.capacity_of_tas()<<std::endl;

    return 0;
}
