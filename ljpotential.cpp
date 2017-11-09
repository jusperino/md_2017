#include "ljpotential.hpp"
#include <iostream>


real LjPotential::force(Particle &p, Particle &q) {
	//calculate additional Force vector by looping over the dimensions
	real epsilon = 1;
	real signum = 1;
	
	real dist = distance(p, q);
    real potential = 4 * epsilon * pow(signum/dist,6) * ( pow(signum/dist,6)-1 );

    real factor=24 * epsilon * sqr(1/dist) * pow(signum/dist,6) * (1-2* pow(signum/dist,6));
    
	for (int i = 0; i<DIM; i++) {
		real add_force = factor * (q.x[i] - p.x[i]);
		p.F[i] += add_force;
	}
    return potential;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
