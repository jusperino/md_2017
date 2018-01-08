#include "world.hpp"
#include "cell.hpp"
#include <stdexcept>
#include <sstream>
#include <string>
#include <cmath>
#include "subdomain.hpp"

World::World(Subdomain &S) : name("unknown"),t(0),delta_t(0),t_end(0),e_kin(0),e_pot(0),e_tot(0),S(S) {
    // empty constructor
}

void World::read_Parameter(const std::string &filename) {

    // create input filestream
    std::ifstream parfile(filename.c_str());

    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error(
            "read_Parameter(): Can't open file '" + filename + "' for reading."
            );

    // helper strings
    std::string line, option, value;

    // read file till eof
    while (parfile.good()) {
        // read line from file
        getline(parfile,line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option and value from stringstream
        strstr >> option;
        // check options and read values
        if (option=="delta_t")
            strstr >> delta_t;
        if (option=="t_end")
            strstr >> t_end;
        if (option=="name")
            strstr >> name;

        if (option=="length"){
            for (int i=0;i<DIM;i++){
                strstr >> length[i];
            }
        }

        if (option=="upper_border"){
            for (int i=0;i<DIM;i++){
                strstr >> value;
                if (value=="unknown"){
                    upper_border[i]=unknown;
                }
                else if (value == "leaving"){
                    upper_border[i]=leaving;
                }
                else{
                	upper_border[i]=periodic;
                }
            }
        }

        if (option=="lower_border"){
            for (int i=0;i<DIM;i++){
                strstr >> value;
                if (value=="unknown"){
                	lower_border[i]=unknown;
                }
                else if (value == "leaving"){
                	lower_border[i]=leaving;
                }
                else{
                	lower_border[i]=periodic;
                }
            }
        }

        if (option=="cell_r_cut"){
            strstr >> cell_r_cut;
        }

        if (option=="num_procs"){
            for (int i=0; i<DIM; ++i){
                strstr >> S.N_p[i];
            }

        }

        if (option=="epsilon"){
			strstr >> epsilon;
		}

		if (option=="sigma"){
			strstr >> sigma;
		}

		if (option=="output_interval"){
			strstr >> output_interval;
		}

		

        option="";
    }

    // calculate cell division parameters from freshly read r_cut
    // obtain parameters for the world (global) and this subdomain (local)
    int global_cell_count = 1;
    int local_cell_count = 1;

    for (int i=0; i<DIM; ++i){
        cell_N[i] = std::max(floor(length[i]/cell_r_cut), 1.0);
        S.N_p[i] = cell_N[i]/S.N_p[i];
        local_cell_count *= S.N_p[i];
        global_cell_count *= cell_N[i];
        cell_length[i] = length[i]/cell_N[i];
    }

    // construct the world grid
    for (int i = 0; i<global_cell_count; ++i) {
        Cell c;
        cells.push_back(c);
    }

    generate_subdomain_cells();
    // close file
    parfile.close();
}

void World::generate_subdomain_cells() {
    // calculate global sequential indices of the cells in this subdomain
    // for each of those, calculate the global sequential indices of the adjacent cells

    // generate displacement modifiers dimension-wise
    std::vector<std::vector<int>> displacement;

    for (int j1 = -1; j1<2; ++j1) {
        for (int j2 = -1; j2<2; ++j2) {
            for (int j3 = -1; j3<2; ++j3) {
                    std::vector<int> comb(3);
                    comb[0] = j1;
                    comb[1] = j2;
                    comb[2] = j3;
                    displacement.push_back(comb);
            }
        }
    }

    std::vector<int> j = get_subd_dim_index(S.myrank);

    // while at it, fill the subdomain class parameters with their values
    std::copy(j.begin(), j.end(), S.ip);
    for (int i = 0; i<DIM; ++i) {
        S.ip_lower[i] = S.ip[i] - 1;
        S.ip_upper[i] = S.ip[i] + 1;
        S.ic_start[i] = 1;      // assume particle can not pass through an entire cell in one timestep 
        S.ic_stop[i] = S.ic_start[i] + S.N_p[i];
        S.ic_number[i] = S.N_p[i] + 2*S.ic_start[i];
        S.ic_lower_global[i] = (S.ip[i] + 1) * S.N_p[i]; // each subdomain is the same size, as are the cells
    }

    for (int k1 = j[0]*S.N_p[0]; k1 < (j[0]+1)*S.N_p[0]; ++k1) {
        for (int k2 = j[1]*S.N_p[1]; k2 < (j[1]+1)*S.N_p[1]; ++k2) {
            for (int k3 = j[2]*S.N_p[2]; k3 < (j[2]+1)*S.N_p[2]; ++k3) {

                std::vector<int> v {k1,k2,k3};
                int K = get_cell_index(v);
                S.cells.push_back(K);

                std::vector<int> nborcand(3,0);
                for (auto &disp: displacement) {
                    // calculate dim-indices of neighbour-candidates
                    nborcand[0] = k1+disp[0];
                    nborcand[1] = k2+disp[1];
                    nborcand[2] = k3+disp[2];
                    bool valid = true;
                    for (int i = 0; i<DIM; ++i){
                        if (nborcand[i] < 0 || nborcand[i] > cell_N[i]-1) {
                            valid = false;
                        }
                    }
                    if (valid) {
                        // add sequential cell index to list of neighbours of that particular cell
                        cells[K].adj_cells.push_back(get_cell_index(nborcand));
                    }
                }
            }
        }
    }
}

bool World::check_if_outside(Particle &p) {
    // given a particle, check if it is within the simulation area
    bool discard = false;
    for (int i = 0; i<DIM; ++i) {
        if (((p.x[i] >= length[i]) && (upper_border[i] == leaving)) || (((p.x[i] < 0 )) && (lower_border[i] == leaving))) {
            discard = true;
        }
    }
    return(discard);
}

int World::get_cell_index(std::vector<int> &j) {
    // calculate sequential cell index from dimension-wise indices
    int J = j[0];
	for(size_t i = 1; i < DIM; ++i){
		J *= cell_N[i];
		J += j[i];
	}
	return(J);
}

std::vector<int> World::get_subd_dim_index(int J) {
    std::vector<int> j;
    // calculate subdomain dimenion-wise index from sequential index
    if (DIM == 3) {
        j.push_back(J/(S.N_p[1]*S.N_p[2]));
        j.push_back((J - S.N_p[1]*S.N_p[2]*j[0])/S.N_p[2]);
        j.push_back(J - S.N_p[2]*(j[1]+S.N_p[1]*j[0]));
    } else if (DIM == 2) {
        j.push_back(J/S.N_p[1]);
        j.push_back(J - S.N_p[1]*j[0]);
    }
    return j;
}



int World::determine_corr_cell(const Particle &p) {
    // calculate cell coordinates
    int J = std::min(std::floor(p.x[0]/cell_length[0]), cell_N[0]/S.N_p[0]-1.0);
	for (size_t i = 1; i < DIM; ++i){
		int j = std::min(std::floor(p.x[i]/cell_length[i]), cell_N[i]-1.0);
		J *= cell_N[i];
		J += j;
	}

	return J;
}

void World::fill_Cell(const Particle &p){
	// calculate cell index
    int J = determine_corr_cell(p);

	// store particle in correct cell
	cells[J].particles.push_back(p);
}

void World::read_Particles(const std::string &filename) {
     // create input filestream
    std::ifstream parfile(filename.c_str());

    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error(
            "read_Parameter(): Can't open file '" + filename + "' for reading."
            );

    // helper strings
    std::string line, data;

    // set particle counter to zero
    particle_count = 0;

    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);

        // ensure no empty particle is created
        if (line != "") {
            // create a string stream
            std::stringstream strstr;
            // put line into string stream
            strstr << line;

            // create new particle
            Particle p;

            // read particle id and mass from current line
            strstr >> p.id >> p.m;

            // read particle position vector dimension-wise from current line
            for (int i = 0; i<DIM; ++i){
                strstr >> p.x[i];
            }

            // read particle velocity vector dimension-wise from current line
            for (int i = 0; i<DIM; ++i){
                strstr >> p.v[i];
            }


            // initialise particle force to zero dimension-wise
            for (int i = 0; i<DIM; ++i){
                p.F[i] = 0;
            }

            // store initiated particle in cell if it is not outside the sim area
            if (check_if_outside(p) == false) {
                fill_Cell(p);
                particle_count++;
            }
        }
    }
    // close file
    parfile.close();
}


std::ostream& operator << (std::ostream& os, World& W) {
    os << "t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end <<" Number of Particles=" << W.particle_count << " Number of Cells=" << W.cells.size() << std::endl;

    std::stringstream str_length, str_upper_border, str_lower_border;


    for (int i=0;i<DIM;i++) {
        str_length << "length["<<i<<"]=" << W.length[i]<< " ";
        str_upper_border << "upper_border["<<i<<"]=" << W.upper_border[i]<< " ";
        str_lower_border << "lower_border["<<i<<"]=" << W.lower_border[i]<< " ";
    }
    os << str_length.str() << std::endl << str_upper_border.str() << std::endl << str_lower_border.str() << std::endl;
    return os;
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
