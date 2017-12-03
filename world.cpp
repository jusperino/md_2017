#include "world.hpp"
#include "cell.hpp"
#include <stdexcept>
#include <sstream>
#include <string>
#include <cmath>

World::World() : name("unknown"),t(0),delta_t(0),t_end(0),e_kin(0),e_pot(0),e_tot(0)
{
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
                if (value=="leaving"){
                    upper_border[i]=leaving;
                }
                else{
                    upper_border[i]=unknown;
                }
            }
        }

        if (option=="lower_border"){
            for (int i=0;i<DIM;i++){
                strstr >> value;
                if (value=="leaving"){
                    lower_border[i]=leaving;
                }
                else{
                    lower_border[i]=unknown;
                }
            }
        }

        if (option=="cell_r_cut"){
            strstr >> cell_r_cut;
        }
        option="";
    }

    // calculate cell division parameters from freshly read r_cut

    int cell_count = 1;
    for (int i=0; i<DIM; ++i){
        cell_N[i] = floor(length[i]/cell_r_cut);
        cell_count *= cell_N[i];
        cell_length[i] = length[i]/cell_N[i];
    }

    // construct the grid by creating N1*N2(*N3) empty cells and adding them to the cells vector
    for (int i = 0; i<cell_count; ++i) {
        Cell c;
        cells.push_back(c);
    }

    generate_adj_cells()

    // close file
    parfile.close();
}

void World::generate_adj_cells(){
    // (hideously) generate the set of neighbour cells for each cell (thankfully only once)

    //start by constructing the dimensional indices of the neighbours relative to a specified cell (e.g. {-1,-1,0}, {1,0,-1} etc
    std::vector<std::vector<int>> displacement;

    for (int j1 = -1; j1<2; ++j1) {
        for (int j2 = -1; j2<2; ++j2) {
            if (DIM == 3){
                for (int j3 = -1; j3<2; ++j3) {
                    std::vector<int> comb(3);
                    comb[0] = j1;
                    comb[1] = j2;
                    comb[2] = j3;
                    displacement.push_back(comb);
                }
            } else {
                std::vector<int> comb(3);
                comb[0] = j1;
                comb[1] = j2;
                displacement.push_back(comb);
            }
        }
    }


    // using these displacement vectors, calculate the neighbours
    for (int j1 = 0; j1<cell_N[0]; ++j1) {
        for (int j2 = 0; j2<cell_N[1]; ++j2) {
            if (DIM == 3){
                for (int j3 = 0; j3<cell_N[2]; ++j3) {
                    // consider each cell using these loops
                    std::vector<int> v(3,0);
                            v[0] = j1;
                            v[1] = j2;
                            v[2] = j3;
                    int J = get_cell_index(v);

                    std::vector<int> nborcand(3,0);
                    for (auto &disp: displacement) {
                        // calculate dim-indices of neighbour-candidates
                        nborcand[0] = j1+disp[0];
                        nborcand[1] = j2+disp[1];
                        nborcand[2] = j3+disp[2];
                        bool valid = true;
                        for (int i = 0; i<DIM; ++i){
                            //check whether the neighbour is a valid cell (if in any dimension the index is greater than the no of cells in that direction, it cannot be valid)
                            if (nborcand[i] < 0 || nborcand[i] > cell_N[i]-1) {
                                valid = false;
                            }
                        }
                        if (valid) {
                            // add sequential cell index to list of neighbours of that particular cell
                            cells[J].adj_cells.push_back(get_cell_index(nborcand));
                        }
                    }
                }
            } else {
                for (auto &disp: displacement) {
                    std::vector<int> nborcand(2,0);
                    nborcand[0] = j1+disp[0];
                    nborcand[1] = j2+disp[1];
                    bool valid = true;
                    for (int i = 0; i<DIM; ++i){
                        if (nborcand[i] < 0 || nborcand[i] > cell_N[i]) {
                            valid = false;
                        }
                    }
                    if (valid) {
                        std::vector<int> v = {j1,j2};
                        cells[get_cell_index(v)].adj_cells.push_back(get_cell_index(nborcand));
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

void World::fill_Cell(const Particle p){
	// calculate cell coordinates
	std::vector<int> j(DIM,0);
	for (size_t i = 0; i < DIM; ++i){
		j[i] = std::floor(p.x[i]/cell_length[i]);
	}

    std::cout << j[0] << " " << j[1] << " " << j[2] <<std::endl;
	// calculate cell index
    int J = get_cell_index(j);
    std::cout << J << std::endl;

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
    os << "t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end <<" Number of Particles=" << W.particle_count<< std::endl;

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
