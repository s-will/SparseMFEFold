#ifndef BASE_HH
#define BASE_HH

#include <iostream>
#include <utility>

//! type of energy
typedef int energy_t;


template<class T1,class T2>
std::ostream &
operator << (std::ostream &out, const std::pair<T1,T2> &x) {
    out<<"("<<x.first<<","<<x.second<<")";
    return out;
}


#endif // BASE_HH
