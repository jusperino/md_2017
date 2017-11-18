    // open xyz file
    std::string xyz_filename = W.name + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyz_filename.c_str());
    // and tell the world
    std::cout << "Opened " << xyz_filename << " for writing." << std::endl;
}
	// close the xyz file
		if ( xyz.is_open() )
	        xyz.close();
}
void Observer::output_xyz()
{
    // write configuration of particles into the filestream, according to .xyz-format
    xyz
		// number of atoms
        << W.particles.size() << "\n"
		// comment line
        << "Zeit: " << W.t << "\n";
    	// x,y,z coordinates
		for (size_t i = 0; i<W.particles.size(); i++){
			// element (here: "H" for every particle)
	        xyz << "H" << "\t"
				<< W.particles[i].x[0] << "\t"
				<< W.particles[i].x[1] << "\t";

			// check if there is a third coordinate to enter, else enter 0 in order to be applicable for both DIM 2 and DIM 3
			if (DIM == 3){
				xyz << W.particles[i].x[2]
					<< std::endl;
			}
			else xyz << 0 << std::endl;
		}
}
    output_xyz();
