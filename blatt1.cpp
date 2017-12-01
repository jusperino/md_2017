#include <cstdlib>
#include <iostream>
#include <cstdio>

#include "world.hpp"
#include "gravitypotential.hpp"
#include "velocityverlet.hpp"
#include "observer.hpp"


void debug(World &W, GravityPotential &Pot, int a, int b) {
	std::string p1 = W.particles[a].id;
	std::string p2 = W.particles[b].id;
	std::cout << "Distance between " << p1 <<  " and " << p2 << ": "  << Pot.distance(W.particles[a], W.particles[b]) << std::endl;
	std::cout << "Potential between " << p1 <<  " and " << p2 << ": "  << Pot.force(W.particles[a], W.particles[b]) << std::endl;
	std::cout << "Force on " << p1 << ": "  << W.particles[a].F[0] << ", " << W.particles[a].F[1] << std::endl;
}

void print_particle(Particle &p){
	std::cout << p.id << "\tMass: " << p.m << "\n";
	std::cout << "x0: " << p.x[0] << "\tx1: " << p.x[1] << "\n";
	std::cout << "v0: " << p.v[0] << "\tv1: " << p.v[1] << "\n";
	std::cout << "F0: " << p.F[0] << "\tF1: " << p.F[1] << "\n";
	std::cout << "F_old0: " << p.F_old[0] << "\tF_old1: " << p.F_old[1] << "\n\n";
}

int main(int argc, char *argv[]) {

  // check arguments
  if (argc < 2) {
    std::cerr << "error: missing arguments" << std::endl;
    std::cerr << "usage: " << std::endl
              << "\t" << argv[0] << " parameterfile particledatafile" << std::endl;
    return EXIT_FAILURE;
  }


  // instanciate Potential
  GravityPotential Pot;

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

  // instanciate timediscretization
  // remark: & is used to get the address of Pot
  VelocityVerlet Verlet(W, &Pot, O);

  // run the simulation
  Verlet.simulate();

  /*for(auto p: W.particles){
    print_particle(p);
  }*/

  debug(W, Pot, 1, 0);

  return EXIT_SUCCESS;
}

// vim:set et sts=2 ts=2 sw=2 ai ci cin:
