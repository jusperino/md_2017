#include <cstdlib>
#include <iostream>
#include <cstdio>

#include "world.hpp"
#include "ljpotential.hpp"
#include "velocityverlet.hpp"
#include "observer.hpp"


void print_particle(Particle &p){
	std::cout << p.id << "\tMass: " << p.m << "\n";
	std::cout << "x0: " << p.x[0] << "\tx1: " << p.x[1] << "\tx2: " << p.x[2] << "\n";
	//std::cout << "v0: " << p.v[0] << "\tv1: " << p.v[1] << "\n";
	//std::cout << "F0: " << p.F[0] << "\tF1: " << p.F[1] << "\n";
	//std::cout << "F_old0: " << p.F_old[0] << "\tF_old1: " << p.F_old[1] << "\n\n";
}

int main(int argc, char *argv[]) {

    // check arguments
    if (argc < 2) {
    std::cerr << "error: missing arguments" << std::endl;
    std::cerr << "usage: " << std::endl
                << "\t" << argv[0] << " parameterfile particledatafile" << std::endl;
    return EXIT_FAILURE;
    }


    // instantiate Potential
    LjPotential Pot;

    // create World
    World W;

    // read Parameters
    W.read_Parameter(argv[1]);

    // read Particles
    W.read_Particles(argv[2]);

    // print World configuration
    std::cout << W << std::endl;


    // print_particle(W.particles[0]);

    //std::cout << W.particles[1].id << std::endl;


    // create the Observer
    Observer O(W);

    // instantiate timediscretization
    // remark: & is used to get the address of Pot
    VelocityVerlet Verlet(W, &Pot, O);

    // run the simulation
    Verlet.simulate();

    /*for(auto p: W.particles){
    print_particle(p);
    }*/


    return EXIT_SUCCESS;
}

// vim:set et sts=2 ts=2 sw=2 ai ci cin:
