#ifndef _WORLD_HPP
#define _WORLD_HPP

#include "defines.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

///Border Type
enum BorderType{unknown = 0, leaving = 1};

/**
 * @brief the world class holds all information of the simulation environment
 */
class World {
public:
    World();

    /**
     * @brief read the world parameters from the given parameter file
     *
     * parameter file example
     * \code
     * delta_t 0.1
     * t_end 1.0
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    void read_Parameter(const std::string &filename);

    /**
     * @brief read the particles from the given data file
     *
     * @param filename filename of the particle data file
     */
    void read_Particles(const std::string &filename);

    /**
     * @brief calculate index of cell from dimension-wise indices
     *
     * @param vector of dimension-wise indices
     */
    int get_cell_index(const std::vector<int> j);

    /**
    * @brief write particle into cell according to coordinates
    *
    * @param p particle to be stored in cell
    */
    void fill_Cell(const Particle p);

    // data structures
    /// Name of the simulated world
    std::string name;
    /// Current time
    real t;
    /// Timestep
    real delta_t;
    /// End of simulation
    real t_end;
    /// world length
    real length[DIM];
    /// world upper border flags
    BorderType upper_border[DIM];
    /// lower world border flags
    BorderType lower_border[DIM];
    /// kinetic energy
    real e_kin;
    /// potential energy
    real e_pot;
    /// total energy
    real e_tot;
    /*
     *	/// Vector of particles
     *	std::vector<Particle> particles;
     **/
    /// cells
    std::vector<Cell> cells;
    /// Number of cells in every dimension
    int cell_N[DIM];
    /// length of cells
    real cell_length[DIM];
    /// r_cut used for calculation of the cell length
    real cell_r_cut;
    /// number of currently simulated particles
    int particle_count;

};

/**
 * @brief a ostream operator for the World class
 *
 * @param os stream object
 * @param W the world
 *
 * @return resulting stream object
 */
std::ostream& operator << (std::ostream& os, World& W);

#endif // _WORLD_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
