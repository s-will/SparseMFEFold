/*
  Zuker with about half of the sugar.  Sparse variant of an
  RNA-folding algorithm "in-between Nussinov and Zuker" (interior loop
  energies, but no proper multi-loop energies)
  
  Demonstration of sparsification and sparse trace-back. For the
  purpose of trace-back, this algorithm is fundamentally more complex
  than the Nussinov algorithm.

  Recursions:
  G(i,j) = min { min_i<=k<=j  G(i,k-1) + C(k,j),
                 C(i,j) 
	       }
  C(i,j) = min { G(i+1,j-1) + NonILoopPenalty if pair_type(S[i],S[j])>0,
                 HairpinE(i,j),
		 min_kl C(i,j)+ILoopE(i,j,k,l)
               }
*/

#include <LocARNA/matrices.hh>

#include <limits>

#include <vector>
#include <list>
#include <map>

#include <cstring>

extern "C" {
#   include "ViennaRNA/data_structures.h"
#   include "ViennaRNA/utils.h"
#   include "ViennaRNA/fold_vars.h"
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

    
    LocARNA::Matrix<energy_t> C_;
    std::vector<energy_t> G_;
    
    typedef std::pair<size_t,energy_t> entry_t;
    typedef std::vector< entry_t > entry_list_t;
    std::vector< entry_list_t > CL_;

    std::string structure_;

    /**
     * @brief Trace arrow
     * 
     * Describes a trace arrow from i,j to k_,l_ of give type
     * type_. The sourcce (i,j) is not represented in the data
     * structure. However, each trace arrow is associated to exactly
     * one source.
     */
    struct trace_arrow_t {
	char type_; //!< type of the arrow (pointing to 'C' or 'G')
	size_t k_; //!< target row of arrow
	size_t l_; //!< target column of arrow
	size_t source_ref_count_; //!< counts how many trace arrows point to the source
	
	/**
	 * @brief construct by target coordinates
	 * @param type
	 * @param k
	 * @param l
	 */ 
	trace_arrow_t(char type,
	     size_t k,
	     size_t l
	     )
	    : type_(type),
	      k_(k),
	      l_(l),
	      source_ref_count_(0)
	{}
	
	/**
	 * @brief empty c'tor
	 */ 
	trace_arrow_t() {}
    };
    
    
    typedef std::map<size_t ,trace_arrow_t> trace_arrow_row_map_t;
    typedef std::map< size_t, trace_arrow_row_map_t > trace_arrow_map_t;

    trace_arrow_map_t trace_arrow_;
    
    /**
     * Get target of trace arrow by source
     *
     * @param i source row index
     * @param j source column index
     */
    const trace_arrow_t &
    trace_arrow_from(size_t i, size_t j) const {
	return trace_arrow_.find(i)->second.find(j)->second;
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
	trace_arrow_map_t::const_iterator row = trace_arrow_.find(i);
	return row != trace_arrow_.end()
	    &&
	    row->second.find(j) != row->second.end();
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
    register_trace_arrow(size_t i, size_t j,char type, size_t k, size_t l) {
	//std::cout << "register_trace_arrow "<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
	trace_arrow_[i][j] = trace_arrow_t(type,k,l);
	
	inc_source_ref_count(k,l);
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
	trace_arrow_map_t::iterator row = trace_arrow_.find(i);
	if (row == trace_arrow_.end()) return;

	trace_arrow_row_map_t::iterator it=row->second.find(j);
	if (it == row->second.end()) return;
	
	trace_arrow_t &ta=it->second;
	
	ta.source_ref_count_++;
    }

    /**
     * Garbage collect trace arrow
     *
     * if count = 0 then remove trace arrow; recursively decrement
     * targets and remove if count drops to 0
     */
    void gc_trace_arrow(trace_arrow_map_t::iterator row, trace_arrow_row_map_t::iterator it) {
	const trace_arrow_t &ta=it->second;
	
	if (ta.source_ref_count_ == 0) {
	    // get trace arrow from the target if the arrow exists
	    trace_arrow_map_t::iterator target_row = trace_arrow_.find(ta.k_);
	    if (target_row != trace_arrow_.end()) {
	    
		trace_arrow_row_map_t::iterator target_it=target_row->second.find(ta.l_);
		if (target_it != target_row->second.end()) {
	    
		    trace_arrow_t &target_ta=target_it->second;
		    assert(target_ta.source_ref_count_>0);
		    target_ta.source_ref_count_--;
		    
		    gc_trace_arrow(target_row,target_it);
		}
	    }
	    
	    row->second.erase(it);
	    if (row->second.size()==0) {
		trace_arrow_.erase(row);
	    }
	}
    }
    
    /** 
     * Garbage collect all trace arrows in given row
     * 
     * @param row_index index of the row 
     */
    void
    gc_row_of_trace_arrows( size_t row_index ) {
	trace_arrow_map_t::iterator row = trace_arrow_.find(row_index);
	if (row == trace_arrow_.end()) return;
	
	for (trace_arrow_row_map_t::iterator it=row->second.begin();
	     row->second.end() != it;
	     ++it) {
	    gc_trace_arrow(row,it);
	}
    }
    
public:
    HalfZuker(const std::string &seq)
	: NonILoopPenalty_(50),
	  ILoopBonus_(-300),
	  seq_(seq),
	  n_(seq.length()),
	  S_(encode_sequence(seq.c_str(),0)),
	  S1_(encode_sequence(seq.c_str(),1)),
	  params_(scale_parameters())
    {
	make_pair_matrix();
	
	C_.resize(MAXLOOP,n_+1);
	G_.resize(n_+1,0);
	
	/* init candidate lists */
	CL_.resize(n_+1);
	for (size_t j=n_; j>0; --j) {
	    CL_[j].push_back( entry_t(j,0) );
	}
    }


    void
    trace_back() {
	structure_.resize(n_+1,'.');

	/* Traceback */
	trace_G(1,n_);
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
    compute_G(size_t i, size_t max_j) {
	//std::cout << "Compute G " <<i<<" "<<max_j<<std::endl; 
	
	for ( size_t j=i-1; j<i+TURN+1; j++ ) { G_[j]=0; }
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
	
	    energy_t d = INF;
	    
	    for ( entry_list_t::iterator it = CL_[j].begin(); 
		  CL_[j].end()!=it && it->first>=i ; ++it ) {
		d = std::min( d, G_[it->first-1] + it->second );
	    }
	
	    G_[j] = d;
	}
	
	// std::cout << "G["<<i<<"]["<<i<<".."<<max_j<<"]: " ;
	// for ( size_t j=i; j<=max_j; j++ ) { 
	//     std::cout << G_[j] << " ";
	// }
	// std::cout << std::endl;
    }

    /** 
     * Trace from G entry
     * 
     * @param i row index
     * @param j column index
     * @param[out] structure result structure
     * pre: structure is string of size (n+1)
     * pre: G contains values of row i in interval i..j
     */
    void
    trace_G(size_t i, size_t j) {
	// std::cout << "Trace G "<<i<<" "<<j<<std::endl;
	if (i+TURN+1>=j) return;
	
	size_t k=j+1;
	
	for ( entry_list_t::iterator it = CL_[j].begin(); 
	      CL_[j].end()!=it && it->first>=i;
	      ++it ) {
	    energy_t d_it = G_[it->first-1] + it->second;
	    
	    if (G_[j] == d_it) {
		k = it->first;
		break;
	    }
	}
	
	//std::cout << i<<" "<<j<<" "<<k<<std::endl;
	assert(i<=k && k<=j);
	
	// don't recompute G here, since i is not changed
	trace_G(i,k-1);
    
	if (k!=j) {
	    trace_C(k,j);
	}
    }

    /** 
     * Trace from C entry
     * 
     * @param i row index
     * @param j column index
     * @param[out] structure result structure
     */
    void
    trace_C(size_t i, size_t j) {
	//std::cout << "Trace C "<<i<<" "<<j<<std::endl;
	assert (i+TURN+1<=j);
	
	structure_[i]='(';
	structure_[j]=')';
	
	if (exists_trace_arrow_from(i,j)) {
	    const trace_arrow_t &arrow = trace_arrow_from(i,j);
	    if(arrow.type_=='C') {
		trace_C(arrow.k_,arrow.l_);
	    } else if (arrow.type_=='G') {
		compute_G(arrow.k_,arrow.l_);
		trace_G(arrow.k_,arrow.l_);
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
		
		for ( entry_list_t::iterator it = CL_[j].begin(); CL_[j].end()!=it; ++it ) {
		    d = std::min( d, G_[it->first-1] + it->second );
		}
		

		int ptype_closing = pair_type(i,j);
		
		energy_t c = INF;
		size_t i_mod=i%MAXLOOP;

		if(ptype_closing>0) {

		    c = HairpinE(i,j);
		    C_(i_mod,j) = c;

		    size_t best_l=0;
		    size_t best_k=0;

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
			    
			    energy_t c_kl=C_(k_mod,l) + ILoopE(ptype_closing,i,j,k,l) + ILoopBonus_;
			    if ( c_kl < c ) {
				c = c_kl;
				best_l=l;
				best_k=k;
			    }
			}
		    }
	    
		    /*
		      std::cout << i << " " << j << " " << CL[j].size()
		      << " " << best_l << " " << best_k
		      << " " << c << " " << d << std::endl;
		    */
		
		    energy_t c_non_iloop = G_[j-1] + NonILoopPenalty_; // G_(i+1,j-1) is still in G_[j-1] !
		    //register trace arrows from C-entry (i,j); even for non-candidates 
		    if ( c_non_iloop < c ) {
			c = c_non_iloop;
			register_trace_arrow(i,j,'G',i+1,j-1);
		    } else {
			if (best_l!=0) { /* not hairpin */
			    register_trace_arrow(i,j,'C',best_k,best_l);
			}
		    }
		    
		    if ( c < d ) {
			// got a candidate => register candidate
			CL_[j].push_back( entry_t(i, c) );
			inc_source_ref_count(i,j); // <- always keep arrows starting from candidates 
		    }
		}
		
		C_(i_mod,j) = c;
		G_[j] = std::min(c,d);
	    }
	    
	    // Clean up trace arrows in i+MAXLOOP
	    gc_row_of_trace_arrows( i + MAXLOOP );
	    
	}
	
	std::cout << "MFE: \t"<<(G_[n_]/100.0)<<std::endl;
    }

    size_t
    num_of_candidates() const {
	size_t c=0;
	for(size_t i=1;i<=n_;i++) {
	    c += CL_[i].size();
	}
	return c;
    }
    size_t
    num_of_tas() const {
	size_t c=0;
	for(size_t i=1;i<=n_;i++) {
	    trace_arrow_map_t::const_iterator row=trace_arrow_.find(i); 
	    if (row!=trace_arrow_.end()) { 
		c += row->second.size();
	    }
	}
	return c;
    }
    
};

int
main(int,char**) {

    std::string seq;
    std::getline(std::cin,seq);
    
    HalfZuker hz(seq);

    std::cout << seq << std::endl;
    std::cout << "Len:\t"<<seq.length()<<std::endl<<std::endl;
    
    hz.forward_evaluation();

    hz.trace_back();

    std::cout << "Unsparse matrix entries: "<< 2* (seq.length()*(seq.length()+1)/2) <<std::endl;
    std::cout << "G:\t"<<seq.length()<<std::endl;
    std::cout << "C:\t"<<seq.length()*MAXLOOP<<std::endl;
    std::cout << "Cands:\t"<<hz.num_of_candidates()<<std::endl;
    std::cout << "TAs:\t"<<hz.num_of_tas()<<std::endl;

    return 0;
}