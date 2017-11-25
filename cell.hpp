/*
 * cell.hpp
 *
 *  Created on: 25.11.2017
 *      Author: LSchiefelbein
 */

#ifndef CELL_HPP_
#define CELL_HPP_

#include "particle.hpp"
#include "world.hpp"
#include <vector>

class Cell{
public:

	/// The cell contains particles ...
	std::vector<Particle> particles;
};



#endif /* CELL_HPP_ */
