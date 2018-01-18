#include "trace_arrow.hh"


TraceArrows::TraceArrows(size_t n)
    : n_(n),
      ta_count_(0),
      ta_avoid_(0),
      ta_erase_(0),
      ta_max_(0)
{}

void
TraceArrows::resize(size_t n) {
    trace_arrow_.resize(n);
}

void
TraceArrows::gc_trace_arrow(size_t i, size_t j) {

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


void
TraceArrows::gc_row( size_t i ) {
    assert(i<=n_);

    for (size_t j=1; j<=n_ ; j++) {
	if (! trace_arrow_[i].exists(j)) continue;
	gc_trace_arrow(i,j);
    }
}

void
TraceArrows::compactify() {
    for ( auto &x: trace_arrow_ ) {
	if (x.capacity() > 1.2 * x.size()) {
	    x.reallocate();
	}
    }
}

size_t
TraceArrows::number() const {
    size_t c=0;
    for ( auto &x: trace_arrow_ ) {
	c += x.size();
    }
    return c;
}

size_t
TraceArrows::capacity() const {
    size_t c=0;
    for ( auto &x: trace_arrow_ ) {
	c += x.capacity();
    }
    return c;
}
