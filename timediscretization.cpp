#include "timediscretization.hpp"

TimeDiscretization::TimeDiscretization(Subdomain& S, World& _W, Potential& _Pot, Observer &_O) : S(S), W(_W), Pot(_Pot), O(_O)
{
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
