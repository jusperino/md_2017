#include "ljpotential.hpp"
#include <iostream>

template <class T>
T IntPow (T a, int b){
	T result = a;
	for(int i = 1; i < b; i++){
		result *= a;
	}
	return result;
}

LjPotential::LjPotential(World& _W, real& epsilon, real& sigma) : Potential(_W),epsilon(epsilon),sigma(sigma) {
    // empty constructor
}

real LjPotential::force(Particle &p, Particle &q) {
    //berechne Lennard-Jones-Potential und Distanz zwischen den Partikeln
	real dist = Potential::distance_2(p, q);
    real potential = 4 * epsilon * IntPow(sigma*sigma/dist,3) * (IntPow(sigma*sigma/dist,3)-1);

	//berechnet richtungsunabhängigen Vorfaktor
    real factor = 24 * epsilon * (1/dist) * IntPow(sigma*sigma/dist,3) * (1 - 2*IntPow(sigma*sigma/dist,3));

    //calculate additional Force vector by looping over the dimensions

	for (int i = 0; i<DIM; ++i) {
		real add_force = factor * (Potential::distance_DIM(p,q,i));
		p.F[i] += add_force;
	}
    return potential;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
