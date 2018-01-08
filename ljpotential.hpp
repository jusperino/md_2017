#ifndef _LJPOTENTIAL_HPP
#define _LJPOTENTIAL_HPP

#include "potential.hpp"
#include <cmath>


class LjPotential : public Potential {
public:

    /**
     * @brief constructor
     * @param _W World to make use of potential methods
     */
    LjPotential(World& _W);

	/**
     * @brief calculate the force between the two particles and add it to p
     *
     * @param p particle p
     * @param q particl q
     *
     * @return potential energy
     */
	virtual real force(Particle &p, Particle &q);

	real epsilon;
	real sigma;

protected:
};

#endif // _LJPOTENTIAL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
