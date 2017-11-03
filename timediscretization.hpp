#ifndef _TIMEDISCRETIZATION_HPP
#define _TIMEDISCRETIZATION_HPP

#include "world.hpp"
#include "potential.hpp"
#include "observer.hpp"
#include <iostream>

/**
 * @brief Base Class for TimeDiscretization Algorithms
 */
class TimeDiscretization {
public:
    /**
     * @brief constructor
     *
     * @param _W world configuration
     * @param _Pot potential used for force calculation
     * @param _O Observer of the simulation
     */
    TimeDiscretization(World& _W, Potential& _Pot, Observer& _O);

    /**
     * @brief run a single timestep
     *
     * @param delta_t length of the timestep
     */
    virtual void timestep(real delta_t) = 0;

    /**
     * @brief run the simulation
     */
    virtual void simulate() = 0;

    /**
     * @brief calculates the forces affecting the particles at the current time
     */
    virtual void comp_F() = 0;

    /**
     * @brief calculates the new velocity of the particles
     */
    virtual void update_V() = 0;

    /**
     * @brief calculate the new position of all particles according to their velocity
     */
    virtual void update_X() = 0;

protected:
    // data structures
    /// the world where the particles live in
    World &W;
    /// the potential used for force calculation
    Potential& Pot;
    /// the observer of the simulation
    Observer &O;

private:
    TimeDiscretization();
};

#endif // _TIMEDISCRETIZATION_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
