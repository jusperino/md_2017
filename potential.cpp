#include "potential.hpp"
#include <iostream>


real Potential::distance(Particle &p, Particle &q) {
	//calculate euclidean norm of distance vector by looping over the dimensions
	real sqrsum = 0;
	for (int i = 0; i<DIM; i++) {
		sqrsum += sqr(q.x[i] - p.x[i]);
	}
	return sqrt(sqrsum);
}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
