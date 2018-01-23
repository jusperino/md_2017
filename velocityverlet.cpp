#include "velocityverlet.hpp"
#include <string>
#include <sstream>


VelocityVerlet::VelocityVerlet(Subdomain& S, World& _W, Potential& _Pot, Observer& _O) : TimeDiscretization(S,_W,_Pot,_O) {
    // empty constructor
}

VelocityVerlet::VelocityVerlet(Subdomain& S, World& _W, Potential* _Pot, Observer& _O) : TimeDiscretization(S,_W,(*_Pot),_O) {
    // empty constructor
}

void VelocityVerlet::simulate() {
    // while simulation end time not reached

    // zaehler wie lange der letze output her ist
	// wenn output_interval erreicht ist wir er auf null zueruck gesetzt und ein output gemacht
	numberOfTimestepsSinceOutput=W.output_interval;

	t_count = 1;

	// initial calculation of forces
	comp_F();

    // simulate over set timeperiod as long as particles are left
    while (W.t < W.t_end) {
        timestep(W.delta_t);
        t_count++;
    }
}

void VelocityVerlet::timestep(real delta_t) {

    // calculate position
    update_X();
    // sort particles into correct cell according to new position
    update_Cells();
    // calculate forces
    comp_F();
    // calculate velocity
    update_V();

    // increase time
    W.t += delta_t;

    // set total energy
    real e_tot_loc = W.e_pot + W.e_kin;

    // communicate total system energy
    // MPI::COMM_WORLD.Allreduce(&e_tot_loc, &W.e_tot, 1, MPI_DOUBLE, MPI_SUM);

    // notify observer if output_interval is reached
    if (t_count%W.output_interval == 0) {
    	O.notify();
    }
}

void VelocityVerlet::comp_F() {

	// set potential energy to 0 in respect of conservation of energy
	W.e_pot = 0;

	W.clear_BorderCells();
	W.communicate_InsideBorder();

	// calculate forces of particles pairwise and sum up potential energy
    for (auto &c: S.cells){
		for (auto &p: W.cells[c].particles){
            for (int i=0; i<DIM; ++i){
                p.F[i] = 0;
            }
            for (auto &n: W.cells[c].adj_cells){
                for (auto &q: W.cells[n].particles){
                    if (&p != &q){
                        W.e_pot += 0.5*Pot.force(p, q);
                        // multiply potential by 0.5, as system potential energy is the sum of *unordered* pairwise potential
                    }
                }
			}
		}
	}
    W.clear_BorderCells();
}

void VelocityVerlet::update_V() {
	// set kinetic energy to 0 in respect of conservation of energy
	W.e_kin = 0;

	// update velocity v for each particle in each dimension and sum up kinetic energy
    for (auto &c: S.cells){
        for (auto &p: W.cells[c].particles){
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
	 // if so and BorderType is chosen "leaving", remove particle from storage vector
    for (auto &c: S.cells){
        auto p = W.cells[c].particles.begin();
        while (p!= W.cells[c].particles.end()){
            for(size_t i = 0; i<DIM; i++){
                (*p).x[i] = (*p).x[i] + W.delta_t * ((*p).v[i] + 0.5 / (*p).m * (*p).F[i] * W.delta_t);

                // check if BorderType is set "periodic" and if so adjust coordinates
                if(W.lower_border[i] == periodic && (*p).x[i] < 0){
                	(*p).x[i] += W.length[i];
                }

                if(W.upper_border[i] == periodic && (*p).x[i] >= W.length[i]){
                	(*p).x[i] -= W.length[i];
                }

                (*p).F_old[i] = (*p).F[i];
            }
            // check_if_outside also checks if BorderType is chosen "leaving"
            if (W.check_if_outside(*p) == false) {
                ++p;
            } else {
                p = W.cells[c].particles.erase(p);
                W.particle_count--;
            }
        }
    }
}

void VelocityVerlet::update_Cells() {

    for (auto &c: S.cells){
    	// use iterator to skip over moved/deleted particles
        auto p = W.cells[c].particles.begin();
        while (p!= W.cells[c].particles.end()){
            std::vector<int> new_cell_coord (W.determine_cell_coord(*p));
            int new_cell_index = W.get_cell_index(new_cell_coord);
            if (c == new_cell_index) {
            	// do nothing if particle is still in correct cell
                ++p;
            } else {
            	W.cells[new_cell_index].particles.push_back(*p);
            	p = W.cells[c].particles.erase(p);
            }
        }
    }

    W.communicate_OutsideBorder();
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
