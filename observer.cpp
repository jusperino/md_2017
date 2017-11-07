#include "observer.hpp"

Observer::Observer(World &_W) : W(_W)
{
    // open statistics file
    std::string statistics_filename = W.name + ".statistics";
    // open file, overwrite existing files, take no prisioners
    statistics.open(statistics_filename.c_str());
    // and tell the world
    std::cout << "Opened " << statistics_filename << " for writing." << std::endl;

    // open coordinates file
    std::string coordinates_filename = W.name + ".csv";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordinates_filename.c_str());
    // and tell the world
    std::cout << "Opened " << coordinates_filename << " for writing." << std::endl;
}


Observer::~Observer()
{
    // close the statistics file
    if ( statistics.is_open() )
        statistics.close();
	// close the coordinates file
	if ( coordinates.is_open() )
        coordinates.close();
}

void Observer::output_statistics()
{
    // write statistics into the filestream, seperated with tabulars
    statistics
        << W.t << "\t"
        << W.e_pot << "\t"
        << W.e_kin << "\t"
		<< W.e_tot << "\t"
        << std::endl;
}

void Observer::output_coordinates()
{
    coordinates
		<< W.t << "\t";
		for (int i = 0; i<W.particles.size(); i++){
			coordinates << W.particles[i].x[0] << "\t"
						<< W.particles[i].x[1] << "\t"
                        /*<< W.particles[i].v[0] << "\t"
                        << W.particles[i].v[1] << "\t"
                        << W.particles[i].F[0] << "\t"
                        << W.particles[i].F[1] << "\t"
                        << " | " << "\t"*/;
		}
    coordinates << std::endl;
}

void Observer::notify()
{
    // debug output
    //std::cout << "notify()" << std::endl;

    // call output functions
    output_statistics();
    output_coordinates();
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
