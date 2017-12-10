#include "potential.hpp"
#include "world.hpp"
#include <iostream>


real Potential::distance(Particle &p, Particle &q) {
	//calculate euclidean norm of distance vector by looping over the dimensions
	real sqrsum = 0;
	for (int i = 0; i<DIM; i++) {
		real dist = World::distance_DIM(p, q, i);
		sqrsum += sqr(dist);
	}
	return sqrt(sqrsum);
}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
