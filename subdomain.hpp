#ifndef _SUBDOMAIN_HPP
#define _SUBDOMAIN_HPP

#include "world.hpp"
#include "defines.hpp"

class SubDomain {
	public:
	real L[DIM]; 		// Kantenlängen des Gesamtgebiets
	int N_c[DIM]; 		// Zahl der Zellen im Gesamtgebiet
	
	// zusätzliche Parameter fur Parallelisierung ¨
	int myrank; 		// Prozessnummer des lokalen Prozesses
	int numprocs; 		// Anzahl der gestarteten Prozesse
	int ip[DIM]; 		// Position des Prozesses im Prozessgitter
	int N_p[DIM]; 		// Gröoße des Prozessgitters, bzw. Zahl der Teilgebiete
	int ip_lower[DIM]; 	// Prozessnummern der Nachbarprozessoren
	int ip_upper[DIM];
	
	int ic_start[DIM]; 	// Breite der Randbordure, entspricht erstem ¨
						// lokalen Index im Inneren des Teilgebiets
	int ic_stop[DIM]; 	// erster lokaler Index in der oberen Randbordure ¨
	int ic_number[DIM]; // Zahl der Zellen im Teilgebiet mit Randbordure ¨
	real cellh[DIM]; 	// Kantenlängen einer Zelle
	int ic_lower_global[DIM]; 	// globaler Index der ersten
								// Zelle des Teilgebiets
};