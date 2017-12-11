#ifndef _POTENTIAL_HPP
#define _POTENTIAL_HPP

#include "particle.hpp"
#include "world.hpp"
#include <cmath>

/**
 * @brief abstract Potential class
 */
class Potential {
public:
    /**
     * @brief constructor
     * @param _W World to make use of potential methods
     */
    Potential(World& _W);

    /**
    * @brief calculate distance dimension wise
    *
    * @param particles whose distance is to be calculated, dimension
    */
    real distance_DIM(Particle &p, Particle &q, int dim);

    /**
     * @brief calculate the absolute distance between the two particles
     *
     * @param p particle p
     * @param q particl q
     *
     * @return absolute distance
     */
	virtual real distance(Particle &p, Particle &q);

    /**
     * @brief calculate the force between the two particles and add it to p
     *
     * @param p particle p
     * @param q particl q
     *
     * @return potential energy
     */
    virtual real force(Particle &p, Particle &q) = 0;

protected:
    World &W;

private:
    Potential();
};

#endif // _POTENTIAL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
