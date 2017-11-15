#include "velocityverlet.hpp"

VelocityVerlet::VelocityVerlet(World& _W, Potential& _Pot, Observer& _O) : TimeDiscretization(_W,_Pot,_O) {
    // empty constructor
}

VelocityVerlet::VelocityVerlet(World& _W, Potential* _Pot, Observer& _O) : TimeDiscretization(_W,(*_Pot),_O) {
    // empty constructor
}

void VelocityVerlet::simulate() {
    // while simulation end time not reached

	// initial calculation of forces
	comp_F();
	// remove particles initialised outside of sim area
	handle_borders();

    // simulate over set timeperiod as long as particles are left
    while (W.t < W.t_end && W.particles.size() != 0) {
        timestep(W.delta_t);
    }
}

void VelocityVerlet::timestep(real delta_t) {
    // test output
    //std::cout << "timestep()\t" << W.e_tot << std::endl;

    // calculate position
    update_X();
    // remove particles outside of simulation area
    handle_borders();
    // calculate forces
    comp_F();
    // calculate velocity
    update_V();

    // increase time
    W.t += delta_t;

    // set total energy
    W.e_tot = W.e_pot + W.e_kin;

    // notify observer
    O.notify();

}

void VelocityVerlet::comp_F() {

	// set forces and potential energy to 0 in respect of conservation of energy
	for (auto & p: W.particles){
		for(size_t i = 0; i<DIM; i++){
			p.F[i] = 0;
		}
	}
	W.e_pot = 0;

	// calculate forces of particles pairwise and sum up potential energy
    for (auto & it1: W.particles){
		for (auto & it2: W.particles){
			if (&it1 != &it2){
				W.e_pot += 0.5*Pot.force(it1, it2);
				// multiply potential by 0.5, as system potential energy is the sum of *unordered* pairwise potential
			}
		}
	}
}

void VelocityVerlet::update_V() {
	// set kinetic energy to 0 in respect of conservation of energy
	W.e_kin = 0;

	// update velocity v for each particle in each dimension and sum up kinetic energy
    for (auto & p: W.particles){
    	for(size_t i = 0; i<DIM; i++){
    		p.v[i] = p.v[i] + W.delta_t * 0.5 / p.m * (p.F[i] + p.F_old[i]);
    		W.e_kin += 0.5 * p.m * sqr(p.v[i]);
    	}
    }
}

void VelocityVerlet::update_X() {
	 // update position x for each particle in each dimension
	 for (auto & p: W.particles){
	    	for(size_t i = 0; i<DIM; i++){
	    		p.x[i] = p.x[i] + W.delta_t * (p.v[i] + 0.5 / p.m * p.F[i] * W.delta_t);
	    		p.F_old[i] = p.F[i];
	    	}
	   }
}

void VelocityVerlet::handle_borders() {
    // check if any particle has left the simulation area
    // if so, remove the particle entirely from the storage vector
    auto it = W.particles.begin();

    // loop over the particles with an iterator
    // when a particle is removed, the iterator is pushed to the one after the removed one

    while (it != W.particles.end()) {
        bool removed = false;
        auto &p = *it;
        for (int i = 0; i<DIM; ++i){
            if (p.x[i] > W.length[i]) {
                it = W.particles.erase(it);
                removed = true;
                break;
            }
        }
        if (!removed) {
            ++it;
        }
    }
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
