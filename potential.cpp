#include "potential.hpp"
#include "world.hpp"
#include <iostream>
#include <algorithm>


Potential::Potential(World& _W) : W(_W){
    // empty constructor
}

real Potential::distance_DIM(Particle &p, Particle &q, int dim){
	// distance respecting periodic lower border, initialized with maximal distance since smallest calculated distance will be returned
	real dist_lower = W.length[dim];
	// distance respecting periodic upper border, initialized with maximal distance since smallest calculated distance will be returned
	real dist_upper = W.length[dim];
	// distance without consideration of BorderType
	real dist = std::abs(p.x[dim] - q.x[dim]);

	// calculate dist_lower, dist_upper if BorderType is set "periodic"
	if (W.lower_border[dim] == periodic){
		dist_lower = std::min(p.x[dim], q.x[dim]) + W.length[dim] - std::max(p.x[dim], q.x[dim]);
	}

	if (W.upper_border[dim] == periodic){
		dist_upper = std::max(p.x[dim], q.x[dim]) + W.length[dim] - std::min(p.x[dim], q.x[dim]);
	}

	// return smallest calculated distance
	return std::min({dist_lower, dist_upper, dist});
}

real Potential::distance(Particle &p, Particle &q) {
	//calculate euclidean norm of distance vector by looping over the dimensions
	real sqrsum = 0;
	for (int i = 0; i<DIM; i++) {
		real dist = distance_DIM(p, q, i);
		sqrsum += sqr(dist);
	}
	return sqrt(sqrsum);
}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
