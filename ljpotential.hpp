#ifndef _LJPOTENTIAL_HPP
#define _LJPOTENTIAL_HPP

#include "potential.hpp"
#include <cmath>

/**
 * @brief TODO add the documentation
 */
class LjPotential : public Potential {
public:
	
	
	/**
     * @brief calculate the force between the two particles and add it to p
     *
     * @param p particle p
     * @param q particl q
     *
     * @return potential energy
     */
	virtual real force(Particle &p, Particle &q);
};

#endif // _LJPOTENTIAL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
