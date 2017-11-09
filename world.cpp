#include "world.hpp"
#include <stdexcept>
#include <sstream>
#include <string>

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
    while (parfile.good())
    {
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
	if (option=="length")
	{
	    for (int i=0;i<DIM;i++)
	    {
		strstr >> lenght[i];
	    }
	}
	if (option=="upper_border")
	{
	    for (int i=0;i<DIM;i++)
	    {
		strstr >> value;
		if (value=="leaving")
		{
		    upper_border[i]=leaving;
		}
		else
		{
		    upper_border[i]=unknown;
		}
	    }
	}
	if (option=="lower_border")
	{
	    for (int i=0;i<DIM;i++)
	    {
		strstr >> value;
		if (value=="leaving")
		{
		    lower_border[i]=leaving;
		}
		else
		{
		    lower_border[i]=unknown;
		}
	    }
	}
	option="";
    }
    // close file
    parfile.close();
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

    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
		
		// 
		std::string id;
		real masse, x0, x1, v0, v1;
		
		// read particle data from stringstream
        strstr >> id >> masse >> x0 >> x1 >> v0 >> v1;
		
		// create new particle and set particle params
		Particle p;
		p.id = id;
		p.m = masse;
		p.x[0] = x0;
		p.x[1] = x1;
		p.v[0] = v0;
		p.v[1] = v1;
		p.F[0] = 0;
		p.F[1] = 0;
		
		// store initiated particle
        particles.push_back(p);
    }
    // close file
    parfile.close();
}	

std::ostream& operator << (std::ostream& os, World& W) {
    os << "t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end <<" Number of Particles=" << W.particles.size()<< std::endl;

    std::stringstream str_lengh, str_upper_border, str_lower_border;

    
    for (int i=0;i<DIM;i++)
    {
	str_lengh << "lengh["<<i<<"]=" << W.lenght[i]<< " ";
	str_upper_border << "upper_border["<<i<<"]=" << W.upper_border[i]<< " ";
	str_lower_border << "lower_border["<<i<<"]=" << W.lower_border[i]<< " ";
    }
    os << str_lengh.str() << std::endl << str_upper_border.str() << std::endl << str_lower_border.str() << std::endl;
    return os;
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
