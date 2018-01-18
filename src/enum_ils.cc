#include <iostream>

#include <cstring>
#include <cassert>

extern "C" {
#   include "ViennaRNA/pair_mat.h"
#   include "ViennaRNA/loop_energies.h"
}

int
main(int argc,char **argv) {
    paramT *params = scale_parameters();
    make_pair_matrix();

    short *S = encode_sequence("ACGU",0);
    short *S1 = encode_sequence("ACGU",1);

    for(int ptype_closing=1; ptype_closing<=6; ptype_closing++) {
	for(int ptype_enclosed=1; ptype_enclosed<=6; ptype_enclosed++) {

	    for (int l1=0; l1<=3; l1++) {
		for (int l2=0; l2<=3; l2++) {

		    //if (l1==0 && l2==0) {continue;}

		    for (int i=1; i<=4; i++) {

			int e = E_IntLoop(l1,l2,
					  ptype_closing,
					  ptype_enclosed,
					  S1[i],
					  S1[i],
					  S1[i],
					  S1[i],
					  const_cast<paramT *>(params));
			if (e<0) {
			    std::cout
				<< e << ":"
				<<  ptype_closing << " " << ptype_enclosed << " " << l1 << " " <<l2
				<< " " << i << std::endl;
			}
		    }
		}
	    }
	}
    }

}
