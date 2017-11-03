#ifndef _OBSERVER_HPP
#define _OBSERVER_HPP

#include "world.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief an observer for the timediscretization
 */
class Observer {
public:
    /**
     * @brief constructor
     *
     * opens and creates the files written during observation
     *
     * @param _W
     */
    Observer(World& _W);

    /**
     * @brief destructor
     *
     * closes the files written during the obervation
     */
    ~Observer();

    /**
     * @brief notify the observer that the world has changed
     */
    void notify();

    /** 
     * @brief output statistics like kinetic, potential and total energy
     */
    void output_statistics();

    /** 
     * @brief output the coordinates of the particles
     */
    void output_coordinates();

protected:
    /// The world we are observing
    World &W;
    /// Statistics filestream
    std::ofstream statistics;
    /// coordiantes filestream
    std::ofstream coordinates;

private:
    /// Disabled Constructor
    Observer();
};

#endif // _OBSERVER_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
