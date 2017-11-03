#include "velocityverlet.hpp"

VelocityVerlet::VelocityVerlet(World& _W, Potential& _Pot, Observer& _O) : TimeDiscretization(_W,_Pot,_O) {
    // empty constructor
}

VelocityVerlet::VelocityVerlet(World& _W, Potential* _Pot, Observer& _O) : TimeDiscretization(_W,(*_Pot),_O) {
    // empty constructor
}

void VelocityVerlet::simulate() {
	//compute Forces for t=0
	comp_F();
    // while simulation end time not reached
    while (W.t < W.t_end)
    {
        timestep(W.delta_t);
    }
}

void VelocityVerlet::timestep(real delta_t) {
    // test output
    //std::cout << "timestep()" << std::endl;

	//compute position for t=n using forces F_(n) computed in previous loop
	update_X();
	//compute forces for t=n+1
	comp_F();
	//compute velocities for t=n+1 using F_(n+1) from this loop and F_(n) computed in previous loop
	update_V();

	//calculate total system energy
	W.e_tot = W.e_kin + W.e_pot;

	//output total energy
	//std::cout << W.t << ": " << W.e_tot << std::endl;

    // increase time
    W.t += delta_t;

    // notify observer
    O.notify();
}

void VelocityVerlet::comp_F() {
	W.e_pot = 0;
    for (int i = 0; i<W.particles.size(); i++){
		for (int j = i+1; j<W.particles.size(); j++){
            W.e_pot += Pot.force(W.particles[i], W.particles[j]);
		}
	}

}

void VelocityVerlet::update_V() {
	W.e_kin = 0;
    for (int k = 0; k<W.particles.size(); k++){
		real v_abs_squared = 0;
		for (int i = 0; i<DIM; i++){
			W.particles[k].v[i] += (W.particles[k].F[i] + W.particles[k].F_old[i])*W.delta_t/(2*W.particles[k].m);
			v_abs_squared += pow(W.particles[k].v[i], 2);
		}
		W.e_kin += W.particles[k].m * v_abs_squared;
	}
	W.e_kin *= 0.5;
}

void VelocityVerlet::update_X() {
    for (int k = 0; k<W.particles.size(); k++){
		for (int i = 0; i<DIM; i++){
			W.particles[k].x[i] += W.delta_t * W.particles[k].v[i] + (W.particles[k].F[i] * pow(W.delta_t, 2))/(2*W.particles[k].m);
			W.particles[k].F_old[i] =  W.particles[k].F[i];
		}
	}
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
