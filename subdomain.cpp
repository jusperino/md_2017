#include "subdomain.hpp"

Subdomain::Subdomain(const int &_numprocs, const int &_myrank): L(),N_c(),myrank(_myrank),numprocs(_numprocs),ip(),N_p(),Cell_N(),ip_lower(),ip_upper(),ic_start(),ic_stop(),ic_number(),cellh(),ic_lower_global() {
    // empty constructor
}