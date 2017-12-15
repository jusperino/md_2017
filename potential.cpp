#include "potential.hpp"
#include "world.hpp"
#include <iostream>
#include <algorithm>


Potential::Potential(World& _W) : W(_W){
    // empty constructor
}

real Potential::distance_DIM(Particle &p, Particle &q, int i){
	// if borders are not periodic, simply return regular distance
	if (W.lower_border[i] != periodic){
		return q.x[i] - p.x[i];
	} else {
		// if one of the borders is not periodic, make sure that the distance calculated does not make use of a non-periodic border
		// W.length is always greater than any other distance, so it will never win the following comparisons
        real dist_lower = q.x[i] - W.length[i] - p.x[i];
        real dist_upper = q.x[i] + W.length[i] - p.x[i];
		real dist = q.x[i] - p.x[i];

		//determine shortest of 3 options
        real shortest = W.length[i];
		if (std::abs(dist_lower) <= std::abs(dist_upper)){
            shortest = dist_lower;
		} else {
		    shortest = dist_upper;
		}

		if (std::abs(dist) < std::abs(shortest)) {
            shortest = dist;
        }
        return shortest;
	}
}

real Potential::distance_2(Particle &p, Particle &q) {
	//calculate square of euclidean norm of distance vector by looping over the dimensions
	real sqrsum = 0;
	for (int i = 0; i<DIM; i++) {
		real dist = distance_DIM(p, q, i);
		sqrsum += sqr(dist);
	}
	return sqrsum;
}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
