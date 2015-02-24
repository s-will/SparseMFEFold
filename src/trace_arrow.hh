#ifndef TRACE_ARROW_HH
#define TRACE_ARROW_HH

#include "base.hh"
#include "simple_map.hh"
#include <cassert>

/**
 * @brief Trace arrow
 * 
 * Describes a trace arrow from i,j to k_,l_. The source (X,i,j)
 * is not represented in the data structure. However, each trace arrow
 * is associated with exactly one source.  Source and target matrix
 * types are omitted, since we don't need them here, but in more
 * general scenarios, such information has to be maintained.
 */
class TraceArrow {
    unsigned char k_; //!< target row of arrow
    unsigned char l_; //!< target column of arrow
    energy_t energy_; //!<target energy
    uint count_; //!< counts how many trace arrows point to the source
public:
    /**
     * @brief construct by target coordinates
     * @param type
     * @param k
     * @param l
     */
    TraceArrow(size_t i,size_t j,size_t k,size_t l,energy_t e)
	: k_(k-i),l_(j-l),energy_(e),count_(0)
    {}

    /**
     * @brief empty c'tor
     */ 
    TraceArrow() {}
	
    size_t k(size_t i,size_t j) const {return k_+i;}
    size_t l(size_t i,size_t j) const {return j-l_;}
    energy_t target_energy() const {return energy_;}
    size_t source_ref_count() const {return count_;}
    
    void inc_src() {count_++;}
    void dec_src() {count_--;}

};

/** 
 * @brief Collection of trace arrows
 *
 * Stores trace arrows to be accessible by row and col index.  Access
 * by column index is logarithmic. TAs of one row are
 * traversable. Supports garbage collection of TAs. Keeps track of
 * several statistics on TAs.
 */
class TraceArrows {

public:
    typedef SimpleMap< size_t, TraceArrow >       trace_arrow_row_map_t;
    typedef std::vector< trace_arrow_row_map_t >  trace_arrow_map_t;

private:    
    trace_arrow_map_t trace_arrow_;

    size_t n_; // sequence length
    
    size_t ta_count_; // count all generated tas
    size_t ta_avoid_; // count all avoided tas (since they point to candidates)
    size_t ta_erase_; // count all erased tas (in gc)
    size_t ta_max_; // keep track of maximum number of tas, existing simultaneously
    
public:
    
    /**
     * @brief Empty constructor
     */
    TraceArrows(size_t n);
    
    void
    resize(size_t n);
    
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
     * @param srctype source matrix type
     * @param i source row
     * @param j source column
     * @param tgttype target matrix type
     * @param k target row
     * @param l target column
     */    
    void
    register_trace_arrow(size_t i, size_t j,
			 size_t k, size_t l,
			 energy_t e) {
	// std::cout << "register_trace_arrow "<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
	trace_arrow_[i].push_ascending( j, TraceArrow(i,j,k,l,e) );
	
	inc_source_ref_count(k,l);
	
	ta_count_++;
	ta_max_ = std::max(ta_max_,ta_count_);
    }

    /**
     * avoid one trace arrow (for statistics only)
     */
    void
    avoid_trace_arrow() {
	ta_avoid_++;
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

private:
    /**
     * Garbage collect trace arrow
     *
     * if count = 0 then remove trace arrow; recursively decrement
     * targets and remove if count drops to 0
     */
    void
    gc_trace_arrow(size_t i, size_t j);

public:    
    /** 
     * Garbage collect all trace arrows in given row
     * 
     * @param i index of the row 
     */
    void
    gc_row( size_t i );

    /**
     * @brief Compactify heap space
     */
    void
    compactify();

    /** @brief Number of trace arrows
     * @return number
     */
    size_t
    number() const;    
    
    /** @brief Capacity of trace arrows vectors
     * @return capacity
     */
    size_t
    capacity() const;

    size_t size() const {return ta_count_;}
    size_t erased() const {return ta_erase_;}
    size_t avoided() const {return ta_avoid_;}
    size_t max() const {return ta_max_;}
    
};



#endif // TRACE_ARROW_HH
