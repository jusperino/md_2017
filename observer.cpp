#include "observer.hpp"
#include <math.h>

Observer::Observer(World &_W) : W(_W)
{
    // open statistics file
    std::string statistics_filename = W.name + ".statistics";
    // open file, overwrite existing files, take no prisioners
    statistics.open(statistics_filename.c_str());
    // and tell the world
    std::cout << "Opened " << statistics_filename << " for writing." << std::endl;

    /*
    // open coordinates file
    std::string coordinates_filename = W.name + ".csv";
    // open file, overwrite existing files, take no prisoners
    coordinates.open(coordinates_filename.c_str());
    // and tell the world
    std::cout << "Opened " << coordinates_filename << " for writing." << std::endl;
    */

    // open xyz file
    std::string xyz_filename = W.name + ".xyz";
    // open file, overwrite existing files, take no prisoners
    xyz.open(xyz_filename.c_str());
    // and tell the world
    std::cout << "Opened " << xyz_filename << " for writing." << std::endl;
}


Observer::~Observer()
{
    // close the statistics file
    if ( statistics.is_open() )
        statistics.close();
	// close the coordinates file
	if ( coordinates.is_open() )
        coordinates.close();
	// close the xyz file
    if ( xyz.is_open() )
        xyz.close();
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
    coordinates << W.t << "\t";
		for (auto &cell: W.cells){
            for (auto &p: cell.particles){
			coordinates << p.x[0] << "\t"
						<< p.x[1] << "\t";
            // check if there is a third coordinate to enter, else enter 0 in order to be applicable for both DIM 2 and DIM 3
            if (DIM == 3){
                coordinates << p.x[2] << "\t";
			}
			else coordinates << 0 << "\t";
            }
		}
    coordinates << std::endl;
}

void Observer::output_xyz()
{
    // write configuration of particles into the filestream, according to .xyz-format
    xyz
    // number of atoms
    << W.particle_count << "\n"
    // comment line
    << "Zeit: " << W.t << "\n";
    // x,y,z coordinates
    for (auto &cell: W.cells){
        for (auto &p: cell.particles){
            // element (here: "H" for every particle)
            xyz << "H" << "\t"
                << p.x[0] << "\t"
                << p.x[1] << "\t";
            // check if there is a third coordinate to enter, else enter 0 in order to be applicable for both DIM 2 and DIM 3
            if (DIM == 3){
                xyz << p.x[2]
                    << std::endl;
            }
            else xyz << 0 << std::endl;
        }
    }
}

void Observer::notify()
{
    float progress = round (10 * 100 * W.t / W.t_end)/10;
    std::cout << "Progress: " << progress << "%" <<"\r";
    std::cout.flush();

    // call output functions
    output_statistics();
    //output_coordinates();
    output_xyz();
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
