#include "gravitypotential.hpp"
#include <iostream>



real GravityPotential::force(Particle &p, Particle &q) {
	//calculate additional Force vector by looping over the dimensions
	real dist = distance(p, q);
    real potential = -1* p.m * q.m / dist;
	for (int i = 0; i<DIM; i++) {
		real add_force = p.m * q.m * (q.x[i] - p.x[i]) / pow(dist,3);
		p.F[i] += add_force;
	}
    return potential;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
