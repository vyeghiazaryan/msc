/***
 * millipede: QuadEqn.h
 * Added by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_QUADEQN
#define H_MILLIPEDE_QUADEQN

#include <vector>

namespace mp {

namespace QuadEqn {

std::vector<double> roots(double A, double B, double C);
double largest_root(double A, double B, double C);
double smallest_root(double A, double B, double C);

}

}

#endif
