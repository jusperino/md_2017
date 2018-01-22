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

	exch_bord();

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
}

void VelocityVerlet::send_particle(Particle &p, int ip, int ic){
	// send particle p to process ip

	// string buffer for particle data
	std::stringstream strstr;

	// send cell index
	strstr << &ic;

	// send particle id and mass
	const char* id = p.id.c_str();
	strstr << id << p.m;

	for (int i=0; i<DIM; ++i){
		strstr << p.x[i] << p.v[i] << p.F[i] << p.F_old[i];
	}

	const std::string tmp = strstr.str();		//create temp string since strstr.str() gets deleted after call
	const char* msg = tmp.c_str();				//create char to send via MPI
	MPI::COMM_WORLD.Send(&msg, 1, MPI_BYTE, ip, t_count);
}

void VelocityVerlet::receive_particle(){
	// receive particle message from any process
	int ic;
	std::string msg;
	std::stringstream strstr;

	// receive message containing char of particle data
	MPI::COMM_WORLD.Recv(&msg, 1, MPI_BYTE, MPI_ANY_SOURCE, t_count);
	strstr << msg;
	strstr >> ic;

	// receive particle data
	Particle p;

	strstr >> p.id >> p.m;

	for(int i = 0; i < DIM; ++i){
		strstr >> p.x[i] >> p.v[i] >> p.F[i] >> p.F_old[i];
	}

	W.cells[ic].particles.push_back(p);
}

void VelocityVerlet::send_cell(int ic, int ip){
	// strstr stores the particle information as a buffer
	//std::stringstream strstr;

	// send cell index
	//strstr << &ic;
	MPI::COMM_WORLD.Send(&ic, 1, MPI_INT, ip, t_count);
	//MPI::COMM_WORLD.Recv();

	// send number of particles contained by cell
	int number = W.cells[ic].particles.size();
	//strstr << number;
	MPI::COMM_WORLD.Send(&number, 1, MPI_INT, ip, t_count);

	// send particle data
	for(auto &p: W.cells[ic].particles){
		const char* id = p.id.c_str();
		//strstr << id;
		MPI::COMM_WORLD.Send(&id, 1, MPI_BYTE, ip, t_count);

		real m = p.m;
		//strstr << m;
		MPI::COMM_WORLD.Send(&m, 1, MPI_DOUBLE, ip, t_count);

		for (int i = 0; i < DIM; i++){
			real x = p.x[i];
			//strstr << x;
			MPI::COMM_WORLD.Send(&x, 1, MPI_DOUBLE, ip, t_count);
			real v = p.v[i];
			//strstr << v;
			MPI::COMM_WORLD.Send(&v, 1, MPI_DOUBLE, ip, t_count);
			real F = p.F[i];
			//strstr << F;
			MPI::COMM_WORLD.Send(&F, 1, MPI_DOUBLE, ip, t_count);
			real F_old = p.F_old[i];
			//strstr << F_old;
			MPI::COMM_WORLD.Send(&F_old, 1, MPI_DOUBLE, ip, t_count);
		}

	}

	/*const std::string tmp = strstr.str();		//create temp string as strstr.str() gets deleted after call
	const char* msg = tmp.c_str();				//create char to send via MPI
	MPI::COMM_WORLD.Send(&msg, 1, MPI_BYTE, ip, 1);*/

}

void VelocityVerlet::recv_cell(int ip){
	int ic;
	int number;


	// receive cell index
	MPI::COMM_WORLD.Recv(&ic, 1, MPI_INT, ip, t_count);

	// clear all particles in cell
	W.cells[ic].particles.clear();

	// receive number of particles contained by cell
	MPI::COMM_WORLD.Recv(&number, 1, MPI_INT, ip, t_count);

	// receive particle data
	for(int k = 0; k < number; k++){
		Particle p;

		MPI::COMM_WORLD.Recv(&p.id, 1, MPI_BYTE, ip, t_count);
		MPI::COMM_WORLD.Recv(&p.m, 1, MPI_DOUBLE, ip, t_count);
		for(int i = 0; i < DIM; i++){
			MPI::COMM_WORLD.Recv(&p.x[i], 1, MPI_DOUBLE, ip, t_count);
			MPI::COMM_WORLD.Recv(&p.v[i], 1, MPI_DOUBLE, ip, t_count);
			MPI::COMM_WORLD.Recv(&p.F[i], 1, MPI_DOUBLE, ip, t_count);
			MPI::COMM_WORLD.Recv(&p.F_old[i], 1, MPI_DOUBLE, ip, t_count);
		}

		W.cells[ic].particles.push_back(p);
	}

}

void VelocityVerlet::exch_block(std::vector<int> I, std::vector<int> J, int ip){
	// block to be exchanged is spanned by cells I and J
	// they represent intervals [I[DIM],J[DIM]) in every dimension
	// loop over every cell in the block, send and receive cell data
	int MAX_NUMBERS = W.particle_count*(3+DIM*4);
	int n_particles = 0;

	while(I[DIM-1] < J[DIM-1]){
		while(I[DIM-2] < J[DIM-2]){
			while(DIM == 2 || I[DIM-3] < J[DIM-3]){
				int c = W.get_cell_index(I);
				n_particles += W.cells[c].particles.size();
				if(DIM == 3) I[DIM-3]++;
			}
		I[DIM-2]++;
		}
	I[DIM-1]++;
	}

	int n_send = n_particles*(3+DIM*4);
	real msg_send[n_send];

	int i = 0;
	while(I[DIM-1] < J[DIM-1]){
		while(I[DIM-2] < J[DIM-2]){
			while(DIM == 2 || I[DIM-3] < J[DIM-3]){
				int c = W.get_cell_index(I);
				for(auto &p: W.cells[c].particles){
					msg_send[i] = c;
					i++;
					const char* temp = p.id.c_str();
					msg_send[i] = atof(temp);
					i++;
					msg_send[i] = p.m;
					i++;
					for(int d=0; d<DIM; d++){
						msg_send[i] = p.x[d];
						i++;
						msg_send[i] = p.v[d];
						i++;
						msg_send[i] = p.F[d];
						i++;
						msg_send[i] = p.F_old[d];
						i++;
					}
				}
				if(DIM == 3) I[DIM-3]++;

			}
			I[DIM-2]++;
		}
		I[DIM-1]++;
	}

	real msg_recv[MAX_NUMBERS];
	int msg_count;
	MPI_Status status;
	if(S.myrank%2 == 0){
		MPI_Recv(&msg_recv, MAX_NUMBERS, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &msg_count);
		MPI_Send(&msg_send, n_send, MPI_DOUBLE, ip, 0, MPI_COMM_WORLD);
	}
	else{
		MPI_Send(&msg_send, n_send, MPI_DOUBLE, ip, 0, MPI_COMM_WORLD);
		MPI_Recv(&msg_recv, MAX_NUMBERS, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &msg_count);
	}

	int j = 0;
	std::vector<bool> clear (W.cells.size(), true);
	while(j<msg_count){
		Particle p;
		int cell_id = msg_recv[j];
		if(clear[cell_id]){
			W.cells[cell_id].particles.clear();
			clear[cell_id] = false;
		}
		j++;
		p.id = msg_recv[j];
		j++;
		p.m = msg_recv[j];
		j++;
		for(int d=0; d<DIM; d++){
			p.x[d] = msg_recv[j];
			j++;
			p.v[d] = msg_recv[j];
			j++;
			p.F[d] = msg_recv[j];
			j++;
			p.F_old[d] = msg_recv[j];
			j++;
		}
		W.cells[cell_id].particles.push_back(p);
	}
}

void VelocityVerlet::exch_bord(){
	// exchange lower and upper border margins in every dimension
	// for this purpose firstly calculate two cells spanning the border margin block (front, bottom, left and back, top, right)
	// then commit calculated coordinates and neighboring process' numbers to exch_block() in order to send and receive cell data

	std::vector<int> I (DIM);
	std::vector<int> J (DIM);

	if(DIM == 3){
		// lower, x3
		if(S.ip_lower[DIM-1] != -1){ 				// check if there is a neighboring process
			for(int i=0; i < DIM; i++){
				I[i] = S.ic_lower_global[i] + S.ic_start[i];
			}

			J[DIM-1] = S.ic_lower_global[DIM-1] + 2*S.ic_start[DIM-1];
			J[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_stop[DIM-2];
			J[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_stop[DIM-3];

			exch_block(I, J, S.ip_lower[DIM-1]);

		}
		// upper, x3
		if(S.ip_upper[DIM-1] != -1){
			I[DIM-1] = S.ic_lower_global[DIM-1] + S.ic_stop[DIM-1] - S.ic_start[DIM-1];
			I[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_start[DIM-2];
			I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_start[DIM-3];

			for(int i=0; i < DIM; i++){
				J[i] = S.ic_lower_global[i] + S.ic_stop[i];
			}

			exch_block(I, J, S.ip_upper[DIM-1]);
		}
	}

	// lower, x2
	if(S.ip_lower[DIM-2] != -1){
		I[DIM-1] = S.ic_lower_global[DIM-1];
		I[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_start[DIM-2];
		if(DIM == 3) I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_start[DIM-3];

		J[DIM-1] = S.ic_lower_global[DIM-1] + S.ic_number[DIM-1];
		J[DIM-2] = S.ic_lower_global[DIM-2] + 2*S.ic_start[DIM-2];
		if(DIM == 3) J[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_stop[DIM-3];

		exch_block(I, J, S.ip_lower[DIM-2]);
	}

	// upper, x2
	if(S.ip_upper[DIM-2] != -1){
		I[DIM-1] = S.ic_lower_global[DIM-1];
		I[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_stop[DIM-2] - S.ic_start[DIM-2];
		if(DIM == 3) I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_start[DIM-3];

		J[DIM-1] = S.ic_lower_global[DIM-1] + S.ic_number[DIM-1];
		J[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_stop[DIM-2];
		if(DIM == 3) J[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_stop[DIM-3];

		exch_block(I, J, S.ip_upper[DIM-2]);
	}

	// lower, x1
	if(S.ip_lower[DIM-3] != -1){
		I[DIM-1] = S.ic_lower_global[DIM-1];
		I[DIM-2] = S.ic_lower_global[DIM-2];
		I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_start[DIM-3];

		J[DIM-1] = S.ic_lower_global[DIM-1] + S.ic_number[DIM-1];
		J[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_number[DIM-2];
		J[DIM-3] = S.ic_lower_global[DIM-3] + 2*S.ic_start[DIM-3];

		exch_block(I, J, S.ip_lower[DIM-3]);
	}

	// upper, x1
	if(S.ip_upper[DIM-3] != -1){
		I[DIM-1] = S.ic_lower_global[DIM-1];
		I[DIM-2] = S.ic_lower_global[DIM-2];
		I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_stop[DIM-3] - S.ic_start[DIM-3];

		J[DIM-1] = S.ic_lower_global[DIM-1] + S.ic_number[DIM-1];
		J[DIM-2] = S.ic_lower_global[DIM-2] + S.ic_number[DIM-2];
		I[DIM-3] = S.ic_lower_global[DIM-3] + S.ic_stop[DIM-3];

		exch_block(I, J, S.ip_upper[DIM-3]);
	}
}




// vim:set et sts=4 ts=4 sw=4 ai ci cin:
