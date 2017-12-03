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

    // simulate over set timeperiod as long as particles are left
    while (W.t < W.t_end && W.particle_count != 0) {
        timestep(W.delta_t);
    }
}

void VelocityVerlet::timestep(real delta_t) {
    // test output
    //std::cout << "timestep()\t" << W.e_tot << std::endl;

    // calculate position
    update_X();
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

	// set potential energy to 0 in respect of conservation of energy
	W.e_pot = 0;

	// calculate forces of particles pairwise and sum up potential energy
    for (auto &cell: W.cells){
		for (auto &p: cell.particles){
            for (int i=0; i<DIM; ++i){
                p.F[i] = 0;
            }
            for (auto &n: cell.adj_cells){
                for (auto &q: W.cells[n].particles){
                    if (&p != &q){
                        W.e_pot += 0.5*Pot.force(p, q);
                        // multiply potential by 0.5, as system potential energy is the sum of *unordered* pairwise potential
                    }
                }
			}
		}
	}
}

void VelocityVerlet::update_V() {
	// set kinetic energy to 0 in respect of conservation of energy
	W.e_kin = 0;

	// update velocity v for each particle in each dimension and sum up kinetic energy
    for (auto &cell: W.cells){
        for (auto &p: cell.particles){
            for(size_t i = 0; i<DIM; i++){
                p.v[i] = p.v[i] + W.delta_t * 0.5 / p.m * (p.F[i] + p.F_old[i]);
                W.e_kin += 0.5 * p.m * sqr(p.v[i]);
            }
        }
    }
}

void VelocityVerlet::update_X() {
	 // update position x for each particle in each dimension
	 // after updating position,check if particle has left sim area
	 // if so, remove particle from storage vector
    for (auto &cell: W.cells){
        auto p = cell.particles.begin();
        while (p!= cell.particles.end()){
            for(size_t i = 0; i<DIM; i++){
                (*p).x[i] = (*p).x[i] + W.delta_t * ((*p).v[i] + 0.5 / (*p).m * (*p).F[i] * W.delta_t);
                (*p).F_old[i] = (*p).F[i];
            }
            if (W.check_if_outside(*p) == false) {
                ++p;
            } else {
                p = cell.particles.erase(p);
                W.particle_count--;
            }
        }
    }
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
